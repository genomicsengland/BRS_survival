#### Christian J. Bouwens
#### BRSC team
#### visualisations
#### last update: 2023.10.05
def simple_count_plt(
	table, 
	x, 
	ax, 
	colour_pal=None, 
	hue=None,
	mask=False,
	scale=False):
	import seaborn as sns
	import matplotlib.ticker as mtick
	from itertools import chain, repeat, cycle
	import numpy as np

	def adjust_counts(group,scale=scale):
			# masking hightest column when only 1 column < threshold
			# Adjust the original data for those hues and add an identifier
			# Add a column to track adjustments
			group['adjusted'] = False
			if group['level_0'].iloc[0] in hue_below_threshold:
				# Find the index of the row with the highest 'count'
				max_idx = group['count'].idxmax()
				# Subtract the threshold from the highest 'count'
				if scale:
					group.loc[max_idx, 'count'] -= scale_threshold[group['level_0'].iloc[0]]['threshold']
				else:
					group.loc[max_idx, 'count'] -= threshold
				# Mark the row as adjusted
				group.loc[max_idx, 'adjusted'] = True
			# add if else for none hue instances.
			return group
	
	def adjust_counts_without_hue(group):
			group['adjusted'] = False
			below_threshold_count = (group['count'] < threshold).sum()
			if below_threshold_count == 1:
				max_idx = group['count'].idxmax()
				group.loc[max_idx, 'count'] -= threshold
				group.loc[max_idx, 'adjusted'] = True
			return group
	
	threshold = 5
	# hue is for several cohorts / tables to be visualised at once.
	if hue:
		if scale:
			table_count = (table[hue]
				.groupby(table[x])
				.value_counts(normalize=True)
				.reset_index(drop=False, name='count')
				)
			# scale_threshold = [(5/y) for y in table[hue].groupby(table[x]).size()]
			scale_threshold = {i:{'threshold':
				(5/y)} for i,y in zip(
					table[hue].groupby(table[x]).size().index,
					table[hue].groupby(table[x]).size())
					}  
			# we'd have to change the downstream handling to deal with a dictionary though
			
			# is there a more global way of doing this?
			# we only use the threshold if mask = True
			# can we differentiate between scales 
			# threshold = scale_threshold
		else:
			table_count = (table
				.groupby([hue, x])
				.size()
				.reset_index(drop=False, name='count')
				)
		# to prevent identification of the number of masked samples ouot of the total
		# in the case where only one column has been adjusted
		# remove the threshold from the higest column
		# and mark the column for addition of a '+' string.
		## for multiple scale_thresholds (hue=True, scale=True)
		## we end up with two thesholds. 

		# adjust the hue's largest column and add capture which one was adjusted.
		if mask:
			if scale:
				# maybe do this with an iterrows?
				for k in  scale_threshold.keys():
					scale_threshold[k]['count'] = 0
				for _,row in table_count.iterrows():
					if row['count'] < scale_threshold[row['level_0']]['threshold']:
						scale_threshold[row['level_0']]['count'] =+1
				masked_col_count2 = pd.DataFrame(scale_threshold).transpose()
				hue_below_threshold = masked_col_count2[masked_col_count2['count'] == 1].index
				# table_count = table_count.groupby('level_0', group_keys=False).apply(adjust_counts)
			else:
				masked_col_count = table_count[
				table_count['count'] < threshold
				].groupby('level_0').size()
				masked_col_count = masked_col_count.reindex(
					table_count['level_0'].unique(), fill_value=0
					)
				# why are we only grabbing the cases where col count is one?
				# don't we need to mask if multiple columns are below the threshold?
				hue_below_threshold = masked_col_count[masked_col_count == 1].index
			
			table_count = table_count.groupby('level_0', group_keys=False).apply(adjust_counts)
		else:
			table_count['adjusted'] = False
	else:
		if scale:
			table_count = (table[x]
				.value_counts(normalize=True)
				.reset_index(drop=False, name='count')
				.rename({'index':x},axis=1))
			scale_threshold = 5/len(table)
			threshold = scale_threshold
		else:
			table_count = (table
				.value_counts(x)
				.reset_index(drop=False, name='count'))
		# this is where we apply the mask.
		if mask:
			table_count = adjust_counts_without_hue(table_count)
		else:
			table_count['adjusted'] = False
		# how does the mask = False work?
		# table count won't have an adjsuted = False.


	sns.set_theme(style='whitegrid')
	sns.barplot(
		data=table_count,
		x=x,
		y='count',
		palette=colour_pal,
		ax=ax,
		hue=hue,
		)
	# removing bars that need masking.
	# the number of patches do not match the threshold, sometimes
	# there are 2, sometimes there are 3. 
	# we have to repeat scale_threshold for the number of bars divided by 
	# the number of hues.

	# is this usefull? we somehow need to attach the information of the 
	# threshold to the bars.
	if hue:
		# Use the hue-based logic for matching
		hue_levels = ax.get_legend_handles_labels()[1]
		for container, hue_level in zip(
			ax.containers, 
			hue_levels, 
			):
			
			# instances where a hue_level only appears in one hue
			# leads to the downstream_iterrows() to mess up.
			# we need to expand the subset.
			subset = table_count[table_count[hue] == hue_level]
			if len(subset) != len(container):
				# Calculate the difference in lengths
				rep = len(container) - len(subset)
				if subset.empty:
					# If subset is empty, create a placeholder DataFrame with empty values
					expanded_subset = pd.DataFrame([{}] * len(container), columns=table_count.columns)
				else:
					# If subset is not empty, repeat rows to match the container length
					expanded_subset = pd.concat([subset] * (rep // len(subset) + 1), ignore_index=True).iloc[:len(container)]
			else:
				expanded_subset = subset

			if expanded_subset.empty:
				labels = [""] * len(container)
			else:
				labels = []
				for bar, (_, row) in zip(container, expanded_subset.iterrows()):
					height = bar.get_height()
					if scale & mask: 
						threshold = scale_threshold[row['level_0']]['threshold']
					if np.isnan(height):
						labels.append("")
					elif (height < threshold) & mask:
						labels.append("<5")
					else:
						if scale:
							label = f"{height:.2f}"
						else:
							label = f"{int(height)}"
						if row["adjusted"]:
							label += "+"
						labels.append(label)
			ax.bar_label(container, labels=labels, padding=-1, fontsize=8)
	else:
		# Non-hue context: simpler approach
		for i, container in enumerate(ax.containers):
			labels = []
			for bar, (_, row) in zip(container, table_count.iterrows()):
				height = bar.get_height()
				if np.isnan(height):
					labels.append("")
				elif (height < threshold) & (mask):
					labels.append("<5")
				else:
					if scale:
						label = f"{height:.2f}"
					else:
						label = f"{int(height)}"
					if row["adjusted"]:
						label += "+"
					labels.append(label)

			# Apply the labels to the container
			ax.bar_label(
				container,
				labels=labels,
				padding=-1,
				fontsize=8
			)

	if scale:
		ax.yaxis.set_major_formatter(mtick.PercentFormatter(1))
	ax.set_ylabel(None)
	ax.xaxis.set_label_position('top')
	ax.margins(.1,0.25)




def simple_hist_plt(
	table,
	x,
	ax,
	binwidth=5,
	hue=None,
	col_labels=True,
	mask=False,
	scale=False,
	multiple='layer'):
	
	import seaborn as sns
	import matplotlib.ticker as mtick

	sns.set_theme(style='whitegrid')
	if not scale:
		sns.histplot(
			data=table,
			x=x,
			binwidth=binwidth,
			element="bars",
			ax=ax,
			hue=hue,
			multiple=multiple
		)

	if scale:
		sns.histplot(
			data=table,
			x=x,
			binwidth=binwidth,
			element="bars",
			ax=ax,
			hue=hue,
			stat='percent',
			common_norm=False,
			multiple=multiple
		)
		datasize = len(table[x])
	threshold = 5

	# Count the number of masked bars
	masked_bar_count = 0
	for bar in ax.patches:
		if mask and bar.get_height() < threshold:
			masked_bar_count += 1

	# Adjust the tallest bar if only one bar is masked
	if mask and masked_bar_count == 1:
		# Find the tallest bar and adjust it
		tallest_bar_index = max(range(len(ax.patches)), key=lambda i: ax.patches[i].get_height())
		tallest_bar = ax.patches[tallest_bar_index]
		tallest_bar.set_height(tallest_bar.get_height() - 5)

	# Apply masking logic and add labels
	for bar in ax.patches:
		if mask:
			if bar.get_height() < threshold:
				bar.set_height(h=0)

	for i, c in enumerate(ax.containers):
		if mask:
			if scale:
				labels = [
					f'{(v/100):.1%}+' if masked_bar_count == 1 and i == tallest_bar_index else 
					f'{(v/100):.1%}' if v >= threshold or v == 0 else 
					"<5" for v in c.datavalues
				]
			else:
				labels = [
					f'{v}+' if masked_bar_count == 1 and i == tallest_bar_index else 
					v if v >= threshold or v == 0 else 
					"<5" for v in c.datavalues
				]
		else:
			if scale:
				labels = [f'{(v/100):.1%}' for v in c.datavalues]
			else:
				labels = [v for v in c.datavalues]
		if col_labels:
			ax.bar_label(
				c,
				labels=labels,
				padding=-1,
				fontsize=8)

	ax.set_ylabel(None)
	ax.xaxis.set_label_position('top') 
	ax.margins(0, 0.25)



# function to create visualisation for A or AB cohorts.
def vis_cohorts(
	cohorts,
	figure,
	coldict=None, 
	names=None, 
	title=True,
	mask=True, 
	show=False,
	scale=False,
	multiple='layer'
	):
	"""creates a figure of age, ancestry, mortality and sex in a grid pattern.
	the data in the figure can be masked (no counts <5, or not).

	Args:
		cohorts (list): a list of cohort class objects or a single cohort class
		containing all_age, all_ancestry, all_sex and all_mort - generated by
		the Cohort.concat_all() class function.
		figure (matplotlib.pyplot.figure): figure set up my matplotlib.
		coldict (Dictionary): colour mapping per feature of tables, must include
			'sex_col', 'vital_col', 'ancestry_col' as keys.
		names (str): Names to give to the cohorts, if None names will be extracted
			from the cohort class.
		title (str, bool, optional): string for the plot to be titled, 
			if bool a title will be generated from the cohort names.
		mask (bool, optional): masks counts lower than 5. Defaults to True.
		show (bool,optional): prints the plot using plt.show if true, returns
			the fig variable if False. 
		scale (bool,optional): scales the values in the plot to percentages. 
			Defaults to False.
		multiple (str, optional): should multi-cohort plots overlap the age
			histogram ('later') or should they be plot side-by-side ('dodge).
			other options: {“layer”, “dodge”, “stack”, “fill”}.
			Defaults to 'layer'.
	"""
	# cohorts=[group1, group2]
	# coldict=None
	# names=None
	# mask=True
	# scale=True
	# show=True

	import seaborn as sns
	import matplotlib.pyplot as plt
	import pandas as pd

	if not coldict:
		coldict = {
		'all_gel':{
			'dark blue':'#2b2f3b',
			'nhs blue':'#005eb8',
			'cyan':'#07c5f5',
			'green':'#26913d',
			'light green':'#71c52d',
			'yellow':'#ffb300',
			'magenta':'#df007d',
			'grey':'#c3c4bf',
			'white': '#FFFFFF'
			},
		'sex_col':{
			'Male':'#26913D', 
			'Female':'#2B2F3B',
			'Indeterminate':'#c3c4bf'},
		'ancestry_col':{
			'AFR':'#26913D',
			'AMR':'#005eb8',
			'EAS':'#2B2F3B',
			'EUR':'#DF007D',
			'SAS':'#07c5f5',
			'UNA':'#ffb300',
			'UNO':'#c3c4bf',
			},
		'vital_cols':{
			'Alive':'#07c5f5',
			'Deceased':'#c3c4bf'
			}
		}
	sex_cols = coldict['sex_col']
	ancestry_cols = coldict['ancestry_col']
	vital_cols = coldict['vital_cols']

	try:
		iterator = iter(cohorts)
	except:
		mult_cohorts = False
		conc_age = cohorts.all_age
		conc_anc = cohorts.all_ancestry
		conc_sex = cohorts.all_sex
		conc_mort = cohorts.all_mortality
		names = cohorts.name
	else:
		# could optimize memory by using iterator instead of cohorts.
		mult_cohorts=True
		# len_cohorts = len(cohorts)  # to judge if the histogram needs to layer or dodge.
		# if there are no names set for the cohort the script will fail.
		if not names:
			names = [df.name for df in cohorts]
		if len(names) != len(cohorts):
			raise ValueError('No names can be attributed to the cohorts.')
		conc_age = pd.concat([df.all_age for df in cohorts], keys=names).reset_index()
		conc_anc = pd.concat([df.all_ancestry for df in cohorts], keys=names).reset_index()
		conc_sex = pd.concat([df.all_sex for df in cohorts], keys=names).reset_index()
		conc_mort = pd.concat([df.all_mortality for df in cohorts], keys=names).reset_index()
	

	# plt.close()
	# plt.clf()
	# fig = plt.figure()
	# let people decide their own figure size.
	fig = figure
	gs = fig.add_gridspec(ncols=3,nrows=2)
	axs1 = fig.add_subplot(gs[0,:-1])
	axs2 = fig.add_subplot(gs[0,2])
	axs3 = fig.add_subplot(gs[1,:-1])
	axs4 = fig.add_subplot(gs[1,2])
	
	### sex ###
	if mult_cohorts:  # can we scale per level_0 instead of by total?
		# (conc_sex
		# 	.groupby(conc_sex['level_0'])
		# 	.value_counts(normalize=True)
		# 	.reset_index())
		# 	.pipe((simple_count_plt, "data"), x='level_0', y=y, hue=hue))
		simple_count_plt(
			table=conc_sex,
			hue='sex',
			ax=axs2,
			colour_pal=sex_cols,
			x='level_0',
			mask=mask,
			scale=scale,
			)
		sns.move_legend(
			axs2, "upper left",
			bbox_to_anchor=(1, 1), 
			ncol=1,  
			frameon=False,
			# columnspacing=0.8,
			handletextpad=0.2,
		)
	else:
		simple_count_plt(
			table=conc_sex,
			x='sex',
			ax=axs2,
			colour_pal=sex_cols,
			hue=None,
			mask=mask,
			scale=scale
			)

	axs2.set_ylabel(None)
	axs2.set_xlabel('Sex')


	### ancestry ###
	if mult_cohorts:
		simple_count_plt(
			table=conc_anc,
			hue='predicted_ancestry',
			ax=axs3,
			colour_pal=ancestry_cols,
			x='level_0',
			mask=mask,
			scale=scale,
			# estimator=lambda x: sum(x=='level_0')*100.0/len(x)
			)
		sns.move_legend(
			axs3, "lower center",
			bbox_to_anchor=(.5, -0.25), 
			ncol=7, 
			frameon=False,
			columnspacing=0.1,
			handletextpad=0.2,
			title=None
		)
	else:
		simple_count_plt(
			table=conc_anc,
			x='predicted_ancestry',
			ax=axs3,
			colour_pal=ancestry_cols,
			hue=None,
			mask=mask,
			scale=scale
			)
	axs3.set_ylabel(None)
	axs3.set_xlabel('Predicted Ancestry')



	### mortality ###
	if mult_cohorts:
		simple_count_plt(
			table=conc_mort,
			hue='status',
			ax=axs4,
			colour_pal=vital_cols,
			x='level_0',
			mask=mask,
			scale=scale
			# estimator=lambda x: sum(x=='level_0')*100.0/len(x)
			)
		sns.move_legend(
			axs4, "upper left",
			bbox_to_anchor=(1, 1), 
			ncol=1,  
			frameon=False,
			# columnspacing=0.8,
			handletextpad=0.2,
			title=None
		)
	else:
		simple_count_plt(
			table=conc_mort,
			x='status',
			ax=axs4,
			colour_pal=vital_cols,
			hue=None,
			mask=mask,
			scale=scale
			)
	axs4.set_ylabel(None)
	axs4.set_xlabel('Vital status')

	### age ###
	if mult_cohorts:

		simple_hist_plt(
			table=conc_age,
			x='age_at_consent',
			ax=axs1,
			binwidth=5,
			hue='level_0',
			mask=mask,
			scale=scale,
			multiple=multiple
			)
		sns.move_legend(
			axs1,
			"upper right",
			frameon=False,
			title=None,
			handletextpad=0.2,
		)
	else:
		simple_hist_plt(
			table=conc_age,
			x='age_at_consent',
			ax=axs1,
			hue=None,
			mask=mask,
			scale=scale
			)
	axs1.set_ylabel(None)
	axs1.set_xlabel('Age')
	
	# adding this otherwise if names is a string it will be counted as len > 1.
	from gelpack.gel_utils import force_list
	names_list = force_list(names)
	if title is not None:
		# if it is none the plot title will be left out.
		if isinstance(title, str):
			fig.suptitle(title, fontsize=16)
		elif title:
			# create a title from cohort names.
			if len(names_list) > 1:
				fig.suptitle(f'Summary of the {*names_list,} cohorts',fontsize=16)
			else:
				fig.suptitle(f'Summary of the {names} cohort',fontsize=16)

	if show:
		plt.show()
	elif not show:
		return fig

