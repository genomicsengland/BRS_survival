#### Christian J. Bouwens
#### Matthieu Vizuete-Forster
#### Miruna Carmen Barbu
#### Chris Odhams
#### BRSC team
#### generate a GEL cohort from HPO terms and ICD-10 codes
#### can be applied for cancer as well as rare disease cases.
#### last update: 2023.08.08


import pandas as pd
import numpy as np
import warnings
import labkey
from gelpack.gel_utils import lab_to_df


class Cohort(object):

	def __init__(self, featdict, version, name=None):

		# check if the featdict has correctly been generated.
		if not any(
			key in [
				'icd10', 
				'hpo', 
				'cancer_terms', 
				'terms'
				] for key in featdict.keys()
				):
			# if there are no features to create a cohort on - 
			# return an error.
			raise KeyError(
				'featdict does not contain any of the following keys:'
				'icd10, hpo, cancer_terms, terms'
				)


		if 'icd10' in featdict.keys():
			self.icd10s = featdict['icd10']
		if 'hpo' in featdict.keys():
			self.hpo = featdict['hpo']
		if 'terms' in featdict.keys():
			self.dterms = featdict['terms']
		if 'cancer_terms' in featdict.keys():
			self.cterms = featdict['cancer_terms']
		# we could add the cancer disease types here.
		# can we build a cohort based on morphology/histology codes?
		# or let people do that themselves and just import the data from pids?
		self.version = version
		self.pids = {}
		self.name = name

	################# functions building the cohort ######################
	def get_term_pids(self):
		"""get participant_ids associated with the given disease_terms from
		various labkey tables: 
			- rare_diseases_participant_disease
		This function does not look for cancer disease terms.

		Returns:
			pd.DataFrame: a dataframe with participant_ids and normalised
			disease terms.
		"""
		# some disease terms have apostrophes in them, here we make a feeble
		# attempt at correcting this atrocity. 
		## SQL will interpret double apostophe as single.
		term_clean = [i.replace("'", "''") for i in self.dterms]
		term_join = ', '.join(rf"'{i}'" for i in term_clean)
		rd_sql = (f"""
			SELECT
				DISTINCT participant_id,
				normalised_specific_disease
			FROM
				rare_diseases_participant_disease
			WHERE 
				normalised_specific_disease IN ({term_join}) 
			""")
		rare_disease = lab_to_df(
			sql_query=rd_sql,
			dr=self.version
			)
		self.dterm_table = rare_disease
		self.pids['dterm'] = rare_disease['participant_id']
	

	def get_hpo_pids(self):
		"""get participant_ids, normalised_hpo_id and terms associated with the 
		given hpo terms from the labkey table:
			- rare_diseases_participant_phenotype

		Returns:
			pd.DataFrame: a dataframe with participant_ids and normalised
			disease terms.
		"""
		## participant_id, and normalised hpo_id / hpo_term
		hpo_sql = (f'''
			SELECT 
				DISTINCT participant_id,
				normalised_hpo_id,
				normalised_hpo_term,
			FROM
				rare_diseases_participant_phenotype
			WHERE
				hpo_present = 'Yes' AND
				normalised_hpo_id 
					IN ({', '.join(f"'{i}'" for i in self.hpo)}) 
			''')
		hpo_terms = lab_to_df(
			sql_query=hpo_sql,
			dr=self.version
			)
		self.hpo_table = hpo_terms
		self.pids['hpo'] = hpo_terms['participant_id']


	def get_cterm_pids(self):
		"""get participant_ids, cancer_disease_types associated with the 
		given cancer_terms from the labkey table:
			- cancer_participant_disease
		Note: this is not an extensive search / confirmation of the various
		tumours our participants may have gotten. But rather returns
		participants with disease_types reported on enrolment to the 100K genomes
		programme.

		Returns:
			pd.DataFrame: a dataframe with participant_ids and cancer
			disease types.
		"""
		format_ca_types =  ', '.join(f"'{i}'" for i in self.cterms)
		cancer_sql = (f'''
			SELECT 
				DISTINCT participant_id,
				cancer_disease_type
			FROM 
				cancer_participant_disease
			WHERE 
				cancer_disease_type IN ({format_ca_types})
			''')
		cancer_disease = lab_to_df(
			sql_query=cancer_sql,
			dr=self.version
			)
		self.cterm_table = cancer_disease
		self.pids['cterm'] = cancer_disease['participant_id']


	def get_icd10_pids(self, limit=['all']):
		"""retrieve participant_ids and cleaned icd-10 codes associated with the 
		given icd-10 codes from the following labkey tables:
		hes:
			- hes-apc
			- hes-ae
			- hes-op
		mort:
			- mortality
		mental_health:
			- mhsds_medical_history_previous_diagnosis
			- mhsds_provisional_diagnosis
			- mhsds_primary_diagnosis
			- mhsds_secondary_diagnosis
			- mhsds_care_activity
			- mhsds_indirect_activity
		cancer:
			- cancer_invest_sample_pathology
			- cancer_participant_tumour
			- cancer_registry
			- rtds
			- sact
			- av_tumour
		
		the sources queried can be limited using the limit argument.
		'all', 'hes', 'mort', 'mental_health', 'cancer'.

		Returns:
			pd.DataFrame: a dataframe with participant_ids and clean
			icd-10 codes.
		"""

		if not all(
			item in [
				'all', 
				'hes', 
				'mort', 
				'mental_health', 
				'cancer'] for item in limit
				):
			raise ValueError(
				'limit contains values not relevant to get_icd10_pids()'
				)


		collect_icd10_tables=[]

		def clean_icd10(table, icd10s, dcol='diag_all'):
			p=r""+'|'.join(f"({i}[0-9]*)"for i in icd10s)+r""
			table[[*icd10s]] = (table[dcol]
				.reset_index(drop=True)
				.str.extractall(p)
				.groupby(level=0)
				.first()
				)

			table = (pd.melt(
				table,
				id_vars=['participant_id'],
				value_vars=[*icd10s]
				)
				.dropna(axis=0)
				.drop_duplicates(['participant_id','value'])
				.drop(['variable'], axis=1)
				.rename({'value':'code'}, axis=1)
				.reset_index(drop=True)
				)
			return table


		# gather all the sources of icd-10 codes, concatening and cleaning them up
		# in the end we only want the relevant icd-10 codes (and one per participant).
		# keeping in mind this top down approach is not perse relevant for cancer.
		# which prefers a bottom-up approach- where we start from the cancer samples
		#   we have, as opposed to starting from diagnosis codes.
		# not using if limit in ['all', 'hes'] to allow different combination of sources.
		if any(lim in ['all','hes'] for lim in limit):
			hes_tables=[
				('apc'),
				('ae'),
				('op')
				]

			# writing a query for the various HES tables.
			# the replace step is uneccessary - 
			# a carry over from the different table names for diagnosis fields.
			# mutliple ICD10 codes possible; therefore the join.
			query = (f'''
				SELECT 
					DISTINCT x.participant_id, 
					diag_all 
				FROM 
					hes_x as x 
				WHERE
					diag_all LIKE {
						" OR diag_all LIKE ".join(f"'%{i}%'" for i in self.icd10s)}''')  # nested f string awh yeah.

			q_apc, q_ae, q_op = [
				query.replace('x',i) for i in hes_tables
				]

			apc_icd10 = lab_to_df(
				sql_query=q_apc,
				dr=self.version
				)
			ae_icd10= lab_to_df(
				sql_query=q_ae,
				dr=self.version
				)
			op_icd10= lab_to_df(
				sql_query=q_op,
				dr=self.version
				)

			hes_df_con = pd.concat([apc_icd10, ae_icd10, op_icd10])

			if hes_df_con.size > 0:
				hes_icd10_clean = clean_icd10(hes_df_con, self.icd10s)
				collect_icd10_tables.append(hes_icd10_clean)
			else:
				warnings.warn("ICD-10 codes not found in HES tables.")

		if any(lim in ['all','mort'] for lim in limit):
			# extract codes from mortality
			# only way to rename the column to diag_all is by preceding the WHERE
			# clause with a sub select (forcing the SELECT to run before the WHERE)
			# the f-string join iterates the 'WHERE diag_all LIKE ICD-10' clause for
			# all our icd-10 codes.
			mortsql = (f'''
				WITH mort AS
				(
					SELECT
						DISTINCT participant_id,
						icd10_multiple_cause_all AS diag_all
					FROM
						mortality
				)
				SELECT
					participant_id,
					diag_all
				FROM
					mort
				WHERE
					diag_all LIKE {
						" OR diag_all LIKE ".join(f"'%{i}%'" for i in self.icd10s)}''')
			
			mortality = lab_to_df(
				sql_query=mortsql,
				dr=self.version
				)
			if mortality.size > 0:
				mort_icd10_clean = clean_icd10(mortality, self.icd10s)
				collect_icd10_tables.append(mort_icd10_clean)
			else:
				warnings.warn("ICD-10 codes not found in mortality tables.")


		if any(lim in ['all','mental_health'] for lim in limit):
			### mental health
			# same approach as above - binding rows of all the various tables.
			# f-strings do not allow backslash. which is neccessary to escape the '.' 
			# as the icd-10 codes in the mental health data can be formated as:
			# G22.2
			mh_sql = f'''
				WITH mhmd AS
				(
					SELECT
						DISTINCT participant_id,
						ic_eve_primarydiagnosis AS diag_mhmdprim,
						ic_eve_secondarydiagnosis AS diag_mhmdsec,
					FROM
						x
				)
				SELECT
					participant_id,
					diag_mhmdprim,
					diag_mhmdsec
				FROM
					mhmd
				WHERE
					REGEXP_REPLACE(diag_mhmdprim, 'b.', '') LIKE {
						" OR REGEXP_REPLACE(diag_mhmdprim, 'b.', '') LIKE "
							.join(f"'%{i}%'" for i in self.icd10s)
						}
					OR REGEXP_REPLACE(diag_mhmdsec, 'b.', '') LIKE {
						" OR REGEXP_REPLACE(diag_mhmdsec, 'b.', '') LIKE "
						.join(f"'%{i}%'" for i in self.icd10s)
					}
					'''.replace('b', '\\')
			# mental health tables differ slightly on column names.
			mhmd_sql, mhmdds_sql = [
				mh_sql.replace('x',i) for i in ['mhmd_v4_event', 'mhldds_event']
				]

			mhml = lab_to_df(
				sql_query=mhmd_sql,
				dr=self.version
				)
			
			mhldds = lab_to_df(
				sql_query=mhmdds_sql,
				dr=self.version
				)
			mh_con = pd.concat([mhml,mhldds])
			
			if mh_con.size > 0:
				mh = (mh_con
					.melt(
						id_vars='participant_id',
						value_vars=['diag_mhmdsec','diag_mhmdprim'])
					.dropna()
					.rename({'value':'diag'}, axis=1)
					.drop('variable', axis=1)
				)
			else:
				mh=pd.DataFrame()

			mhsds_sql = f'''
				WITH table AS
				(
					SELECT
						DISTINCT participant_id,
						x AS diag,
						z AS diagscheme
					FROM
						y
				)
				SELECT
					participant_id,
					diag,
				FROM
					table
				WHERE
					(diagscheme = 'k' OR diagscheme = '0k') AND
					(REGEXP_REPLACE(diag, '\\.', '') LIKE {
						" OR REGEXP_REPLACE(diag, ';.', '') LIKE "
						.join(f"'%{i}%'" for i in self.icd10s)
					})'''.replace(';', '\\')

			# these tuples match the correct tables, column names and diagschemes
			# to query the various mental health tables with labkey for icd-10
			# codes.
			mhsds_str = [
				(
					'mhsds_medical_history_previous_diagnosis',
					'prevdiag',
					'diagschemeinuse',
					'2'
					),
				(
					'mhsds_provisional_diagnosis',
					'provdiag',
					'diagschemeinuse',
					'2'
					),
				(
					'mhsds_primary_diagnosis',
					'primdiag',
					'diagschemeinuse',
					'2'
					),
				(
					'mhsds_secondary_diagnosis',
					'secdiag',
					'diagschemeinuse',
					'2'
					),
				(
					'mhsds_care_activity',
					'codefind',
					'findschemeinuse',
					'1'
					),
				(
					'mhsds_indirect_activity',
					'codefind',
					'findschemeinuse',
					'1'
				)
				]

			mhsds_sql_repl  = [mhsds_sql
				.replace('x', i[1])
				.replace('y', i[0])
				.replace('z', i[2])
				.replace('k', i[3]) for i in mhsds_str]


			mhsds_diag_list = [
				lab_to_df(
					sql_query=i,
					dr=self.version
					) for i in mhsds_sql_repl]
			if mh.size < 0:
				mhsds_diag_list.append(mh)  # add the 2014-2016 style tables.

			mhsds_diag = pd.concat(mhsds_diag_list).reset_index(drop=True)

			if mhsds_diag.size > 0:
				mhsds_diag['code'] = (mhsds_diag['diag']
					.str.replace('.','', regex=False)
					.str.extract(r'([A-Z][0-9]+)')
				)
				mhsds_diag = (mhsds_diag
					.drop_duplicates(['participant_id','code'])
					.drop(['diag'],axis=1)
					)
				collect_icd10_tables.append(mhsds_diag)
			else:
				warnings.warn("ICD-10 codes not found in mental health tables.")

		# concatenate relevant tables (depending on limit.)
		if any(lim in ['all','cancer'] for lim in limit):
			####### icd-10 cancer-specific
			# we should display some kind of warning this is just to gather all known
			# participants which may have had cancer. But not a proper approach
			# to identify patients for which we have tumour samples.
			# cancer invest sample pathology

			## the R script is filtering out participants who may have multiple cancers here.
			# we are keeping duplicate participant ids if they have different cancer codes.
			path_sql = f'''
				WITH table AS
				(
					SELECT
						participant_id,
						x as code_raw
						
					FROM
						z
				)
				SELECT 
					DISTINCT participant_id,
					code_raw
				FROM
					table
				WHERE 
					(REGEXP_REPLACE(code_raw, '\\.', '') LIKE {
						" OR REGEXP_REPLACE(code_raw, ';.', '') LIKE "
						.join(f"'%{i}%'" for i in self.icd10s)
					})'''.replace(';', '\\')

			c_table_str = [
				(
					'cancer_invest_sample_pathology',
					'primary_diagnosis_icd_code',
					),
				(
					'cancer_participant_tumour',
					'diagnosis_icd_code',
					),
				(
					'cancer_registry',
					'cancer_site'
					),
				(
					'rtds',
					'radiotherapydiagnosisicd'
				),
				(
					'sact',
					'primary_diagnosis'
				),
				(
					'av_tumour',
					'site_icd10_o2'
				)
				]


			cancer_sql_repl  = [path_sql
				.replace('x', i[1])
				.replace('z', i[0]) for i in c_table_str]


			cancer_code_list= [
				lab_to_df(
					sql_query=i,
					dr=self.version
					) for i in cancer_sql_repl]

			cancer_diag = pd.concat(cancer_code_list).reset_index(drop=True)
			if cancer_diag.size > 0:
				cancer_diag['code'] = (cancer_diag['code_raw']
					.str.replace('.','', regex=False)
					.str.extract(r'([A-Z][0-9]+)')
				)
				cancer_diag = (cancer_diag
					.drop_duplicates(['participant_id','code'])
					.drop(['code_raw'],axis=1)
					)
				collect_icd10_tables.append(cancer_diag)
			else:
				warnings.warn("ICD-10 codes not found in cancer tables.")
		# print('length of icd10_tables: ', len(collect_icd10_tables))
		# concat the various icd-10 sources.
		if len(collect_icd10_tables) > 1:
			self.icd10_table = pd.concat(collect_icd10_tables)
			self.pids['icd10'] = pd.concat(collect_icd10_tables)['participant_id']
		elif len(collect_icd10_tables) == 1: 
			self.icd10_table = collect_icd10_tables[0]
			self.pids['icd10'] = collect_icd10_tables[0]['participant_id']
		else:
			self.warnings.warn("No participants found with these ICD-10 codes.")


	################ adding features to the cohort ##################
	# we grab the unique participants here.
	def age(self):
		"""Get the age at consent and current age for each participant in the 
		cohort. This function does go over each source (cancer, icd10, hpo...)
		individually, which may lead to a few duplicate entries being calculated.
		But allows for comparison of the source cohorts downstream.

		Returns:
			age_table dictionary: The keys of which correspond to the keys of the sources,
			for each key a pd.dataframe with 'date_of_consent', 'participant_id',
			'year_of_birth', 'current_age' and 'age_at_consent'. In addition a
			QC flag has been added - this flags instances where the year_of_birth
			was set at a default value, the date of consent occured before 
			the start of the 100,000 Genomes Programme, and if participants
			were over 100 years old at the date of consent.
		"""
		#### Age ####
		self.age_table = {}
		for key, pid in self.pids.items():
			age_sql = (f'''
				SELECT
					participant_id,
					year_of_birth,
					date_of_consent,
				FROM
					participant
				WHERE
					participant_id IN {*pid,}
				''')
			# import datetime

			age_df = lab_to_df(
				sql_query = age_sql,
				dr = self.version
				)

			age_df['date_of_consent'] = (age_df['date_of_consent']
				.apply(
					pd.to_datetime,
					format='%Y-%m-%d'
					)
				.fillna(
					pd.to_datetime('2015-01-01',format='%Y-%m-%d')
				))
			age_df['year_of_birth'] = (age_df['year_of_birth']
				.apply(pd.to_datetime,
					format='%Y'
				))


			def calculate_age(born):
				from datetime import date
				today = date.today()
				return today.year - born.year - ((today.month, today.day) < (born.month, born.day))

			def calculate_age_diag(born, diag):
				return diag.year - born.year - ((diag.month, diag.day) < (born.month, born.day))


			age_df['current_age'] = age_df.apply(
				lambda x:
					calculate_age(x['year_of_birth']),
				axis=1
			)

			age_df['age_at_consent'] = age_df.apply(
				lambda x:
					calculate_age_diag(
						born = x['year_of_birth'],
						diag = x['date_of_consent']),
				axis=1
			)

			# age qc:
			## TODO double check this is correct.
			age_df['age_qc'] = np.where(
				(  
				(age_df['year_of_birth'] == pd.to_datetime(
					'1900-01-01',format='%Y-%m-%d'))
				| (age_df['date_of_consent'] < pd.to_datetime('2000-01-01',format='%Y-%m-%d'))
				| (age_df['age_at_consent'] > 100)
				), False, True)

			self.age_table[key] = age_df


	def ancestry(self):
		"""Get the age at consent and current age for each participant in the 
		cohort. like age(), this function goes over each source (cancer, icd10, hpo...)
		individually, which may lead to a few duplicate entries being calculated.
		But allows for comparison of the source cohorts downstream.

		Returns:
			ancestry_table (dictionary): The keys of which correspond to the keys of 
			the sources, for each key a pd.dataframe with 'participant_id' 
			and 'predicted_ancestry'. ancestries are set by the highest scoring
			predicted ancestry in aggregate_gvcf_sample_stats. If no single 
			score was higher than 0.8 the ancestry is set to unassigned (UNA). 
			If the participant is not part of this table the ancestry is set 
			to unknown (UNO).
		"""
		## ancestry  ##
		self.ancestry_table = {}
		for key, pid in self.pids.items():
			ancestry_sql = (f'''
				SELECT
					participant_id,
					pred_african_ancestries,
					pred_south_asian_ancestries,
					pred_east_asian_ancestries,
					pred_european_ancestries,
					pred_american_ancestries
				FROM
					aggregate_gvcf_sample_stats
				WHERE
					participant_id IN {*pid,}''')

			ancestry = lab_to_df(
				sql_query=ancestry_sql,
				dr=self.version
				)
			ancestry['predicted_ancestry'] = pd.NA
			# idxmax returns the column-name with the greates value for each row.
			ancestry['predicted_ancestry'].fillna(
				ancestry[[
					'pred_african_ancestries',
					'pred_south_asian_ancestries',
					'pred_east_asian_ancestries',
					'pred_european_ancestries',
					'pred_american_ancestries']].idxmax(axis=1),
				inplace=True)

			ancestry.loc[
				(
				(ancestry['pred_african_ancestries'] < 0.8)
				& (ancestry['pred_south_asian_ancestries'] < 0.8)
				& (ancestry['pred_east_asian_ancestries'] < 0.8)
				& (ancestry['pred_european_ancestries'] < 0.8)
				& (ancestry['pred_american_ancestries'] < 0.8)
				),
				'predicted_ancestry'
				] = 'pred_unassigned'

			ancestry['predicted_ancestry'].replace(
				{
				'pred_african_ancestries':'AFR',
				'pred_south_asian_ancestries':'SAS',
				'pred_east_asian_ancestries':'EAS',
				'pred_european_ancestries':'EUR',
				'pred_american_ancestries':'AMR',
				'pred_unassigned':'UNA'
				}, inplace=True
				)
			# capture which participants we are missing so we can assign them unknown
			miss_par = [
				p for p in pid if 
					(p not in list(ancestry['participant_id']))
				]
			if len(miss_par) > 0:
				ancestry = pd.concat([
					ancestry[['participant_id', 'predicted_ancestry']],
					pd.DataFrame({
						'participant_id' : miss_par, 
						'predicted_ancestry' : 'UNO'})
						]
					)
			self.ancestry_table[key] = ancestry


	def sex(self):
		"""query the labkey tables for phenotypic sex per participant_id.
		This data is the participant's stated sex by the clinician at the GMC.
		"""
		## Sex ##
		## TODO add kayotype and illumina ploidy as sex query options.
		self.sex_table = {}

		for key, pid in self.pids.items():
			sex_sql = (f'''
				SELECT
					participant_id,
					participant_phenotypic_sex as sex,
				FROM
					participant
				WHERE
					participant_id IN {*pid,}
				''')
			sex = lab_to_df(
				sql_query=sex_sql,
				dr=self.version
				)
			self.sex_table[key] = sex

	def mortality(self):
		"""Extract the death date from labkey tables for a set of participant_id.

		Args:
			pids (list): list of participant ids to include
			version (str): Data release version

		Returns:
			mortality_table (pd.DataFrame): date of death per participant_id
		"""

		# import date_of_death per participant_id from 'mortality'
		self.mortality_table = {}
		for key, pid in self.pids.items():

			query2 = (
				f'''
				SELECT 
					DISTINCT participant_id, death_date
				FROM
					death_details
				WHERE 
					participant_id IN {*pid,}
				'''
			)
			ons1 = lab_to_df(
				sql_query=query2,
				dr=self.version)
			query3 = (
				f'''
				SELECT 
					participant_id, date_of_death
				FROM 
					mortality
				WHERE 
					participant_id IN {*pid,}
			''')
			ons2 = lab_to_df(
				sql_query=query3,
				dr=self.version)

			ons1.rename(columns={'death_date':'date_of_death'}, inplace=True)
			# cohorts may have every / most participants alive
			# if one of the two tables is empty the merge will fail, hence:
			check_size = [ons1.size > 0, ons2.size > 0]
			if sum(check_size) == 0:  # no data in tables
				warnings.warn("No participants found in mortality tables.")
			elif sum(check_size) == 2: # both tables have data
				ons = (ons1.merge(
						ons2, 
						how='outer')
					.sort_values(
						by='date_of_death', 
						ascending=True)
					.drop_duplicates(
						subset='participant_id', 
						keep='first')
				)
				tmp_mort = pd.merge(pid, ons, on='participant_id', how='left')
				tmp_mort['status'] = np.where(
					tmp_mort['date_of_death'].isna(), 'Alive', 'Deceased')

				self.mortality_table[key] = tmp_mort
			elif sum(check_size) == 1:  # only one of the two lists has got data.
				res = [[ons1, ons2][i] for i, val in enumerate(check_size) if val]
				tmp_mort = pd.merge(pid, res[0], on='participant_id', how='left')
				tmp_mort['status'] = np.where(
					tmp_mort['date_of_death'].isna(), 'Alive', 'Deceased')
				self.mortality_table[key] = tmp_mort  # unlisting the datarame.

				

	################ summarysing the cohort ##################
	def concat_cohort(cls, dc):
		"""concatenate the dataframes in a dictionary, note; the dicitonarys
		should have the same structure / columns. 

		Args:
			dc (dictionary): Dictionary value's are pd.DataFrames, one per Key.

		Returns:
			pd.DataFrame: A concatenated pandas dataframe with unique values in
			the dictionary.
		"""
		concatdc = (pd.concat(dc.values())
			.reset_index(drop=True)
			.drop_duplicates()
			)

		return concatdc
	
	def summarize(cls, pid, anc, mort, age, sex):
		"""gather up cohort data and return summary stats. including:
		- number of unique participants
		- mean age at consent
		- standard deviation of mean age at consent
		- mean current age
		- fraction female
		- fraction diseased
		- fraction of each ancestry.

		Args:
			pid (np.Array): array of participant ids.
			anc (pd.DataFrame): 2 column table with participant_id and 
				predicted ancestry.
			mort (pd.DatFrame): 2 column table with participant_id and 
				date_of_death
			age (pd.DataFrame): Table with age_at_conset and current_age 
				per participant id.
			sex (pd.DataFrame): 2 column table with a sex per participant_id.

		Returns:
			pd.DataFrame: single row DataFrame with summary stats.
		"""
		ancestry_distribution = (anc.predicted_ancestry
			.value_counts(normalize=True)
			.reset_index(name='0', drop=False)
			.T)
		ancestry_distribution.columns = ancestry_distribution.iloc[0]
		ancestry_distribution.drop(ancestry_distribution.index[0], inplace=True)

		# changed the structure of the mortality table so no longer needed,
		# percent_mort = len(
		# 	[x for x in list(pid) if x in list(mort['participant_id'])]
		# 	) / len(pid)

		summary = pd.DataFrame(
			{
			'n_unique_participants':len(pid),
			'mean_consent_age':np.mean(age['age_at_consent']),
			'std_consent_age':np.std(age['age_at_consent']),
			'mean_current_age':np.mean(age['current_age']),
			'std_current_age':np.std(age['current_age']),
			'fraction_female':
				sex['sex'].value_counts(normalize=True)[0],
				  # alphabetically ordered.
			'fraction_deceased': 
				mort['status'].value_counts(normalize=True)[1]
			}, index=[0])
		
		summary_anc = pd.concat(
			[
				summary, 
				ancestry_distribution.reset_index(drop=True)
				], 
			axis=1)

		return summary_anc
	
	def concat_all(self):
		self.all_pids = self.concat_cohort(self.pids)
		self.all_age = self.concat_cohort(self.age_table)
		self.all_ancestry = self.concat_cohort(self.ancestry_table)
		self.all_sex = self.concat_cohort(self.sex_table)
		self.all_mort = self.concat_cohort(self.mortality_table)

	def summary_stats(self):
	
		# collapse data.
		all_pids = self.concat_cohort(self.pids)
		all_age = self.concat_cohort(self.age_table)
		all_ancestry = self.concat_cohort(self.ancestry_table)
		all_sex = self.concat_cohort(self.sex_table)
		all_mortality = self.concat_cohort(self.mortality_table)

		self.summary = self.summarize(
			pid=all_pids,
			anc=all_ancestry,
			mort=all_mortality,
			age=all_age,
			sex=all_sex
			)
		# TODO add missingness.

	def load_icd10_data(self):
		from importlib import resources
		with resources.path("gelpack",'coding19.tsv') as f:
			self.icd10_lookup = pd.read_csv(
				f, 
				sep='\t',
				usecols=['coding','meaning'],
				).rename({'coding':'code'}, axis=1)

# 	).rename({'coding':'code'}, axis=1))
	def ontology_stats(self):
		summ_dict = {}
		for key, pid in self.pids.items():
			summ_dict[key] = self.summarize(
				pid=pid.drop_duplicates(),
				anc=self.ancestry_table[key].drop_duplicates(),
				mort=self.mortality_table[key].drop_duplicates(),
				age=self.age_table[key].drop_duplicates(),
				sex=self.sex_table[key].drop_duplicates()
				)
			summ_dict[key]['source'] = key
		
		self.ont_summary = pd.concat(
			summ_dict.values(), 
			ignore_index=True, 
			names=[key for key in self.pids.keys()]
			)

		# participant counts per searched term / ontology
		# count_dict = {} is already part of summary_ont now (n_unique_participants).
		# vcount should summarize the counts for each resource.		
		self.ont_vcount = {}
		if all([k in self.pids.keys() for k in ['dterm','hpo','cterm','icd10']]):
			raise RuntimeError('No valid cohorts to summarize found')
		else:
			# check which keys have been included in the cohort.
			for key in ['dterm', 'hpo', 'cterm' , 'icd10']:
				if key in self.pids.keys():
					if key == 'dterm':
						self.ont_vcount[key] = (self.dterm_table 
							.drop_duplicates([
								'participant_id',
								'normalised_specific_disease'])
							.normalised_specific_disease
							.value_counts()
						)  
					elif key == 'hpo':
						hpo_uniq = (self.hpo_table
							.drop_duplicates([
								'participant_id',
								'normalised_hpo_id'
								])
								)
						hpo_uniq['term'] = (hpo_uniq['normalised_hpo_id'] 
							+ ' (' 
							+ hpo_uniq['normalised_hpo_term']+')'
							)
						self.ont_vcount[key] = hpo_uniq.term.value_counts()  # return to self.
					# related to cancer -> is this how we would go about creating a
					# cancer cohort?
					# TODO: add a check for icd-10 codes (confirm diagnosis)
					# TODO: add a search for study_abbreviation.
					# TODO: add a ICD-10 -> sample / option.
					elif key == 'cterm':  
						self.ont_vcount[key] = (self.cterm_table
							.drop_duplicates([
								'participant_id', 
								'cancer_disease_type'
								])
							.cancer_disease_type
							.value_counts()
						)
					elif key == 'icd10':
						self.load_icd10_data()
						# return simple codes
						un_icd10 = self.icd10_table.drop_duplicates([
							'participant_id',
							'code'
						])
						full_icd10 = un_icd10.merge(
							self.icd10_lookup, 
							on='code',
							how='left')
						self.ont_vcount[key+'_full'] = (
							full_icd10
								.meaning
								.value_counts(dropna=False)
							)
						simple_icd10 = un_icd10.copy()
						simple_icd10['simple_code'] = (un_icd10['code']
							.str.extract(r'([A-Z][0-9]{2})')
							)
						self.ont_vcount[key+'_simple'] = (
							simple_icd10
								.simple_code
								.value_counts(dropna=False)
								.reset_index()
								)
						# generate a matrix covering the overlap of ICD-10 codes
						# and its association with multiple participants.
						parts = un_icd10.participant_id.value_counts().reset_index()
						multparts = parts[parts['participant_id']>1].iloc[:,0].unique()
						self.icd10_overlap_matrix = (
							simple_icd10.loc[
								simple_icd10['participant_id'].isin(multparts)
								]
							.groupby(['participant_id','simple_code'])
							.size()
							).unstack()
					


## utilities for figures.
## apply those utilities to gather figures.
# if we set the table strings in a dictionary per version - which is loaded 
# upon setting the version -> easy to do version control.

# if we set the table strings in a dictionary per version - which is loaded 
# upon setting the version -> easy to do version control.