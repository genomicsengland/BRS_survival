#!/usr/bin/env python

#### Christian J. Bouwens
#### evaluate lines of treatment and determine outcomes from
#### real world data in cancer participants.
#### last update: 2025.01.23

#########################
# these functions assume a clean input source only containing entries of drugs
# direclty affecting the tumour. So any steroids or co-administered drugs
# should have been removed. 
#########################

import pandas as pd

def identify_treatment_lines_with_mechanism_mapping(
		patient_data, 
		treatment_data, 
		drug_mechanism_dict, 
		gap_days=30):
	"""
	Identifies treatment lines for each participant and maps drugs to their 
		mechanism of action.

	Args:
		patient_data (pd.DataFrame): Contains participant_id and 
			sample_diagnosis_date.
		treatment_data (pd.DataFrame): Contains participant_id, 
			administration_date, and drug_group.
		drug_mechanism_dict (dict): Dictionary mapping drug mechanisms of 
			action to drug names.
		gap_days (int): Number of days to consider a gap significant enough 
			to indicate a new line.

	Returns:
		pd.DataFrame: Treatment lines for each participant with start date, 
			end date, and mapped mechanisms of action.
	"""
	# Create copies to preserve original data
	patient_data_copy = patient_data.copy()
	treatment_data_copy = treatment_data.copy()

	# Ensure correct data types
	patient_data_copy['participant_id'] = patient_data_copy['participant_id'].astype(str)
	patient_data_copy['sample_diagnosis_date'] = pd.to_datetime(
		patient_data_copy['sample_diagnosis_date'], errors='coerce'
		)

	treatment_data_copy['participant_id'] = treatment_data_copy['participant_id'].astype(str)
	treatment_data_copy['administration_date'] = pd.to_datetime(
		treatment_data_copy['administration_date'], errors='coerce'
		)
	treatment_data_copy['drug_group'] = treatment_data_copy['drug_group'].astype(str).str.lower()

	# Map drugs to their mechanisms of action
	drug_to_mechanism = {
		drug.lower(): mechanism
		for mechanism, drugs in drug_mechanism_dict.items()
		for drug in drugs
	}

	treatment_data_copy['mechanism'] = treatment_data_copy['drug_group'].map(
		drug_to_mechanism
		).fillna('unknown')

	# Merge tables
	merged = pd.merge(
		treatment_data_copy, 
		patient_data_copy, 
		on='participant_id', 
		how='inner'
		)

	# Filter treatments after diagnosis
	merged = merged[merged['administration_date'] >= merged['sample_diagnosis_date']]

	# Sort by participant_id and administration_date
	merged = merged.sort_values(
		by=['participant_id', 'administration_date']
		)

	result = []

	# Group by participant to identify treatment lines
	for participant, group in merged.groupby('participant_id'):
		group = group.sort_values('administration_date')
		treatment_lines = []
		current_line = set()
		start_date = None
        # build up lines of treatment based on the mechanism of each drug
		# and the time between administration dates.
		for _, row in group.iterrows():
			if not current_line:
				current_line.add(row['mechanism'])
				start_date = row['administration_date']
			else:
				# Check if mechanism fits in current treatment line's window
				# this step is crucial to merge different mechanisms that
				# are part of the same treatment line.
				if (row['administration_date'] - start_date).days <= gap_days:
					current_line.add(row['mechanism'])
				else:
					# Close the current line
					end_date = group.loc[
						group['administration_date'] 
						<= row['administration_date'] 
						- pd.Timedelta(days=gap_days),
						  'administration_date'].max()
					treatment_lines.append({
						'participant_id': participant,
						'start_date': start_date,
						'end_date': end_date,
						'mechanisms': tuple(sorted(current_line))
					})
					
					current_line = {row['mechanism']}
					start_date = row['administration_date']

		# Append the last line
		if current_line:
			end_date = group['administration_date'].max()
			treatment_lines.append({
				'participant_id': participant,
				'start_date': start_date,
				'end_date': end_date,
				'mechanisms': tuple(sorted(current_line))
			})

		result.extend(treatment_lines)

	result_df = pd.DataFrame(result)
	result_df['start_date'] = pd.to_datetime(result_df['start_date'], errors='coerce')
	result_df['end_date'] = pd.to_datetime(result_df['end_date'], errors='coerce')
	result_df['participant_id'] = result_df['participant_id'].astype(str)

	return result_df

# this is currently limited to the first line of treatment,
# but it doesn't have to be. could easily extend it to capture 
# every line of treatment.
def map_line_of_treatment(line_df, treatment_paths):
	"""
	Maps the line of treatment mechanisms to a predefined treatment path.

	Args:
		first_line_df (pd.DataFrame): Table containing participant_id, 
			start_date, end_date, and mechanisms.
		treatment_paths (dict): Dictionary defining treatment paths and their 
			corresponding mechanism combinations.

	Returns:
		pd.DataFrame: Input table with an additional column 
			'first_line_treatment' mapped from the mechanisms.
	"""
	# Normalize treatment paths for matching
	normalized_paths = {
		key: [tuple(sorted(mech)) for mech in value]
		for key, value in treatment_paths.items()
	}

	# Function to find matching treatment path
	def find_treatment_path(mechanisms):
		mechanisms_sorted = tuple(sorted(mechanisms))
		for path, paths in normalized_paths.items():
			if mechanisms_sorted in paths:
				return path
		return 'Unknown'

	# Apply the mapping function to each row
	line_df['treatment_line'] = line_df['mechanisms'].apply(
		find_treatment_path
		)

	return line_df


# Get first line of treatment
def get_first_line_of_treatment(treatment_lines):
	"""
	Selects the first line of treatment based on start_date for each participant.

	Args:
		treatment_lines (pd.DataFrame): Table containing participant_id, start_date, end_date, and mechanisms.

	Returns:
		pd.DataFrame: Table with the first line of treatment for each participant.
	"""
	return treatment_lines.sort_values(
		by=['participant_id', 'start_date']
		).groupby('participant_id').first().reset_index()


def get_second_line_of_treatment(treatment_lines):
	"""
	Selects the second line of treatment based on start_date for each 
		participant (if it exists).
	Args:
		treatment_lines (pd.DataFrame): Table containing participant_id, 
			start_date, end_date, and mechanisms.
	Returns:
		pd.DataFrame: Table with the second line of treatment for each participant.
	"""
	return treatment_lines.sort_values(
		by=['participant_id', 'start_date']
	).groupby('participant_id').nth(1).reset_index()


def get_first_and_second_line_info(treatment_lines):
	"""Combines first and second line of treatment information, 
	calculating the gap between them.
	Args:
		treatment_lines (pd.DataFrame): Table containing participant_id, 
			start_date, end_date, and mechanisms.
	Returns:
		pd.DataFrame: Table containing participant_id, first line start and 
			end date, if a second line occurs, and the time in days between 
			first and second line (if applicable).
	"""
	# Get first and second treatment lines
	first_lines = get_first_line_of_treatment(treatment_lines)
	second_lines = get_second_line_of_treatment(treatment_lines)
	# Merge to align first and second lines
	merged = pd.merge(
		first_lines, 
		second_lines[['participant_id', 'start_date']],
		on='participant_id', how='left', suffixes=('_first', '_second')
	)
	# Determine if a second line exists
	merged['second_line_exists'] = ~merged['start_date_second'].isna()
	# Calculate the gap between first and second line (if applicable)
	merged['gap_to_second_line'] = (merged['start_date_second'] - merged['end_date']).dt.days
	merged['gap_to_second_line'] = merged['gap_to_second_line'].where(merged['second_line_exists'])
	# Keep only relevant columns
	final_df = merged[['participant_id', 'start_date_first', 'end_date', 'second_line_exists', 'gap_to_second_line']]
	
	return final_df


import plotly.graph_objects as go
import kaleido

def create_sankey_diagram(final_first_line_df, final_second_line_df, lot_surv):
	"""
	Creates a Sankey diagram to visualize treatment paths from first-line treatment 
	to second-line treatment (if any) and finally to survival status.

	Args:
		final_first_line_df (pd.DataFrame): Contains participant_id, start_date, end_date, treatment_line.
		final_second_line_df (pd.DataFrame): Contains participant_id, start_date, end_date, treatment_line.
		lot_surv (pd.DataFrame): Contains participant_id and alive_at_threshold (boolean).
	Returns:
		go.Figure: A Plotly Sankey figure object.
	"""
	# Merge first-line with second-line data
	merged_df = pd.merge(
		final_first_line_df, 
		final_second_line_df[['participant_id', 'treatment_line']], 
		on='participant_id', 
		how='left', 
		suffixes=('_first', '_second')
	)

	# Replace NaN second-line treatments with "No second line"
	merged_df['treatment_line_second'] = merged_df['treatment_line_second'].fillna("No second line")

	# Differentiate second-line treatments from first-line ones
	merged_df['treatment_line_second'] = merged_df['treatment_line_second'].apply(
		lambda x: f"{x} (Second line)" if x != "No second line" else x
	)

	# Merge with survival data
	merged_df = pd.merge(merged_df, lot_surv, on='participant_id', how='left')

	# Convert survival status to "Alive" or "Dead"
	merged_df['survival_status'] = merged_df['alive_at_threshold'].map({True: "Alive", False: "Dead"})

	# Create unique categories (nodes)
	all_categories = pd.concat([
		merged_df['treatment_line_first'],
		merged_df['treatment_line_second'],
		merged_df['survival_status']
	]).unique()

	category_to_index = {category: i for i, category in enumerate(all_categories)}

	# Create source-target link data
	links = []

	# Transition from first-line to second-line treatment
	first_to_second_counts = merged_df.groupby(['treatment_line_first', 'treatment_line_second']).size()
	for (first, second), count in first_to_second_counts.items():
		links.append({
			'source': category_to_index[first],
			'target': category_to_index[second],
			'value': count
		})

	# Transition from second-line to survival status
	second_to_survival_counts = merged_df.groupby(['treatment_line_second', 'survival_status']).size()
	for (second, survival), count in second_to_survival_counts.items():
		links.append({
			'source': category_to_index[second],
			'target': category_to_index[survival],
			'value': count
		})

	# Create Plotly Sankey diagram
	sankey_fig = go.Figure(go.Sankey(
		node=dict(
			pad=20,
			thickness=20,
			label=list(category_to_index.keys())
		),
		link=dict(
			source=[link['source'] for link in links],
			target=[link['target'] for link in links],
			value=[link['value'] for link in links]
		)
	))
	sankey_fig.update_layout(title_text="Treatment Pathway Sankey Diagram", font_size=12)
	
	return sankey_fig