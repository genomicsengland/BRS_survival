#### Christian J. Bouwens
#### Matthieu Vizuete-Forster
#### Miruna Carmen Barbu
#### Chris Odhams
#### BRSC team
#### generate a GEL self from HPO terms and ICD-10 codes
#### can be applied for cancer as well as rare disease cases.
#### last update: 2023.09.12

import pandas as pd
import numpy as np
import warnings
import labkey
from gelpack.gel_utils import lab_to_df, force_list
from functools import reduce


class Cohort(object):

	def __init__(self,
		version, 
		featdict=None, 
		participants=None, 
		platekeys=None, 
		name=None):
		# check if the featdict has correctly been generated.
		# how to handel multiple? (featdict + platekeys?)
		if (featdict is None and participants is None and platekeys is None):
			raise ValueError(
				'either featdict, participant ids or platekeys must be given to create cohort.'
			)
		if featdict:
			if not any(
				key in [
					'icd10', 
					'hpo', 
					'cancer_terms',
					'cancer_abbr',
					'terms'
					] for key in featdict.keys()
					):
				# if there are no features to create a self on - 
				# return an error.
				raise KeyError(
					'featdict does not contain any of the following keys:'
					'icd10, hpo, cancer_terms, terms or cancer_abbr'
					)

			if 'icd10' in featdict.keys():
				self.icd10s = force_list(featdict['icd10'])
			if 'hpo' in featdict.keys():
				self.hpo = force_list(featdict['hpo'])
			if 'terms' in featdict.keys():
				self.dterms = force_list(featdict['terms'])
			if 'cancer_terms' in featdict.keys():
				self.cterms = force_list(featdict['cancer_terms'])
			if 'cancer_abbr' in featdict.keys():
				self.cabbr = force_list(featdict['cancer_abbr'])
			# we could add the cancer disease types here.
			# can we build a self based on morphology/histology codes?
			# or let people do that themselves and just import the data from pids?
		self.featdict=featdict
		self.version = version
		self.pids = {}
		self.feature_tables = {}
		self.sample_tables = {}
		self.name = name
		self.platekeys = None
		self.groups = {}  # to assign intra-cohort groups.

		if participants is not None:
			if participants == 'all_cancer':
				ca_pids = self.get_cancer_parts(dr=version).participant_id
				self.custom_pids(
					pids_lst_file=ca_pids,
					action='include'
				)
			else:
				self.custom_pids(
					pids_lst_file=participants,
					action='include'
					)
		
		if platekeys is not None:
			plate_pids = self.get_pids_per_platekey(
					platekeys,dr=version
					).participant_id
			self.custom_pids(
				pids_lst_file=plate_pids,
				action='include'
				)
			self.platekeys = platekeys

	################# functions building the self ######################
	def custom_pids(self, pids_lst_file, action='include'):
		"""
		Take a list of participant ids or a csv/tsv file with a 
		'participant_id' column and return the participant ids in a format
		we can include in self.pids. we are keeping the return ambiguous as
		we may use this function to exclude pids from the cohort as well.

		Args:
			pids_lst_file (str or list): a filepath to a tsv/file with a header
			'participant_id' and a column of participant ids. Or: a list of 
			participant ids.
			aciton (str): 'include' or 'exclude': shoudl the participants be 
			included in self.pids or should they be excluded from self.pids and
			the present features (sex, mortality, ancestry, age)?

		Raises:
			ValueError: if the csv/tsv has no header.
			ValueError: if the delimiters of the tsv/csv file are not ',' or
			'\t'.
			ValueError: if the csv/tsv has no header participant_id.

		Returns:
			_type_: _description_
		"""
		# allow custom participants ids to be added to the cohort.
		if isinstance(pids_lst_file, list) and not isinstance(pids_lst_file, str):
			ids = pd.Series(pids_lst_file, name='participant_id')
		elif isinstance(pids_lst_file, pd.Series):
			ids = pids_lst_file
			ids.name ='participant_id'
		elif pids_lst_file.endswith(('.tsv', ".csv")):
			import csv
			# check if the file has one column and a header:
			with open(pids_lst_file, 'r') as f:
				try:
					csv.Sniffer().has_header(f.read(1024))
				except csv.Error:
					raise ValueError(
						'The participants tsv/csv file has no header.'
					)
				dialect = csv.Sniffer().sniff(f.read(1024), [",","\t"]) 
				if dialect.delimiter not in [',','\t']:
					raise ValueError(
						'The delimiters of the participants tsv/csv file are not supported'
					)
				f.seek(0)
				reader = csv.DictReader(f, dialect=dialect)
				lines = [line for line in reader]
				if 'participant_id' not in lines[0].keys():
					raise ValueError(
						'The participants tsv/csv file has no header: participant_id'
					)

			part_table = pd.read_csv(pids_lst_file, sep=dialect.delimiter)
			ids = part_table['participant_id']
		else:
			warnings.warn('Custom participant_ids not in list or tsv/csv.')

		if action == 'include':
			self.pids['custom'] = ids
		elif action == 'exclude':
			for keys, table in self.pids.items():
				self.pids[keys] = pd.Series(
					[x for x in table if x not in list(ids)], 
					name='participant_id')
				if keys == 'dterm':
					self.dterm_table = self.dterm_table.loc[
						~self.dterm_table['participant_id'].isin(ids)
					]
				elif keys == 'cterm':
					self.cterm_table = self.cterm_table.loc[
						~self.cterm_table['participant_id'].isin(ids)
					]
				elif keys == 'icd10':
					self.icd10_table = self.icd10_table.loc[
						~self.icd10_table['participant_id'].isin(ids)
					]
				elif keys == 'hpo':
					self.hpo_table = self.hpo_table.loc[
						~self.hpo_table['participant_id'].isin(ids)
					]
				elif keys == 'cabbr':
					self.cabbr_table = self.cabbr_table.loc[
						~self.cabbr_table['participant_id'].isin(ids)
					]
			# remove the data from the feature teables too if they exist.
			# just using if self.age_table is not None: does not work, 
			# hence the existance of feature_tables , we can track what
			# features are in the class, and easily iterate over them.
			for name, feature in self.feature_tables.items():
				if name == 'age':
					for keys, table in self.age_table.items():
						table.drop(table[table['participant_id'].isin(ids)].index,
						inplace=True)
				if name == 'mortality':
					for keys,table in self.mortality_table.items():
						table.drop(table[table['participant_id'].isin(ids)].index,
						inplace=True)
				if name == 'sex':
					for keys,table in self.sex_table.items():
						table.drop(table[table['participant_id'].isin(ids)].index,
						inplace=True)
				if name == 'ancestry':
					for keys,table in self.ancestry_table.items():
						table.drop(table[table['participant_id'].isin(ids)].index,
						inplace=True)
				if name == 'omics_sample_data':
					for keys,table in self.omics_sample_data.items():
						table.drop(table[table['participant_id'].isin(ids)].index,
						inplace=True)
				if name == 'gmc_registration':
					for keys,table in self.gmc_registration.items():
						table.drop(table[table['participant_id'].isin(ids)].index,
						inplace=True)
					
						# also recalculate the sample counts/location
						

			# remove the participants from the concatenated table:
			# except -> when all_data wasn't created yet.
			try:
				self.all_data = self.all_data.loc[
					~self.all_data['participant_id'].isin(ids)
				]
			except AttributeError:
				pass
			except NameError:
				pass
	

	def add_to_group(self, id, group):
		if group not in self.groups.keys():
			self.groups[group] = []
		# confirm the id is in the cohort.
		if [x.eq(id) for x in self.pids.values()][0].any():
			self.groups[group].append(id)


	def get_cancer_parts(cls, dr):
		sqlstr = (
			'''
			SELECT
				participant_id,
			FROM
				cancer_analysis
			''')
		return lab_to_df(sql_query=sqlstr, dr=dr)


	##### TODO #####
	# if the cohort has been created on platekeys this function is redundant.
	# if we are going from participant_id -> platekeys we definitly need to 
	# make sure we grab only samples of the relevant cancer type.
	# and remove the rest.
	def limit_cohort_to_cancer_type(self, pids, dr):

		# filter the cohort to only include participants where we have recruited
		# the participant for a particular cancer type. (we have the tumour
		# Sample of the cancer of interest).

		# depends on featdict
		# confirm disease type. -> translate to either study_abb or GEL disease type.
		# get all cancer samples from cancer_analysis
		# ca = get_cancer_parts(dr=self.version)

		self.ca_sample_data()

		# for each sample in ca_sample_data, 
		# check if its a sample recruited for featdict.
		# this requires the featdict to be translated to diseaes_types / study_abbreviations.
		from gelpack.gel_utils import translateicd
		feat_icd10_trans = translateicd(self.featdict['icd10'])
		self.cancer_samples['icd10']

		# only
		return None

	
	def get_pids_per_platekey(cls, platekeys, dr):
		sqlstr = (
			f'''
			SELECT
				participant_id,
				plate_key
			FROM
				plated_sample
			WHERE
				plate_key in {*platekeys,}
			''')
		return lab_to_df(sql_query=sqlstr, dr=dr)


	### TODO ###
	def get_platekey_per_pid(cls, pid, dr):
		sqlstr =()


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
				cancer_disease_type AS disease_type
			FROM 
				cancer_participant_disease
			WHERE 
				cancer_disease_type IN ({format_ca_types})
			''')
		cancer_disease = lab_to_df(
			sql_query=cancer_sql,
			dr=self.version
			)
		
		cancer_analysis_sql =(f'''
			SELECT
				DISTINCT participant_id,
				disease_type
			FROM
				cancer_analysis
			WHERE
				disease_type IN ({format_ca_types})
		''')

		cancer_analysis = lab_to_df(
			sql_query=cancer_analysis_sql,
			dr=self.version
			)
		disease_table = pd.concat(
			[cancer_disease,cancer_analysis], axis=0
			).drop_duplicates(keep='first')

		self.cterm_table = disease_table
		self.pids['cterm'] = disease_table['participant_id']


	def get_cabbr_pids(self):
		"""Query cancer analysis for participants with a given study abbreviation
		This is limited to the GEL 100K cancer cohort.
		"""
		format_abbr_types = ', '.join(f"'{i}'" for i in self.cabbr)
	
		
		cancer_analysis_sql =(f'''
			SELECT
				DISTINCT participant_id,
				study_abbreviation
			FROM
				cancer_analysis
			WHERE
				study_abbreviation IN ({format_abbr_types})
		''')

		cancer_analysis = lab_to_df(
			sql_query=cancer_analysis_sql,
			dr=self.version
			)

		self.cabbr_table = cancer_analysis
		self.pids['cabbr'] = cancer_analysis['participant_id']


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
			- mhmd_v4_event
			- mhldds_event
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


	################ adding features to the self ##################
	def age(self):
		"""Get the age at consent and current age for each participant in the 
		self. This function does go over each source (cancer, icd10, hpo...)
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
			age_df['age_qc'] = np.where(
				(  
				(age_df['year_of_birth'] == pd.to_datetime(
					'1900-01-01',format='%Y-%m-%d'))
				| (age_df['date_of_consent'] < pd.to_datetime('2000-01-01',format='%Y-%m-%d'))
				| (age_df['age_at_consent'] > 100)
				), False, True)

			self.age_table[key] = age_df
			self.feature_tables['age'] = self.age_table


	def ancestry(self):
		"""Get the age at consent and current age for each participant in the 
		self. like age(), this function goes over each source (cancer, icd10, hpo...)
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
			self.feature_tables['ancestry']=self.ancestry_table


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
			self.feature_tables['sex'] = self.sex_table


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
				''')
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
			
			query4 = (
				f'''
				SELECT
					participant_id, date_of_death
				FROM
					rare_diseases_pedigree_member
				WHERE
					participant_id IN {*pid,}
			''')
			rd_ons = lab_to_df(
				sql_query=query4,
				dr=self.version
				)
			rd_ons = rd_ons.dropna()
			ons1.rename(columns={'death_date':'date_of_death'}, inplace=True)
			# cohorts may have every / most participants alive
			# if one of the two tables is empty the merge will fail, hence:
			check_size = [ons1.size > 0, ons2.size > 0, rd_ons.size>0]
			if sum(check_size) == 0:  # no data in tables
				warnings.warn("No participants found in mortality tables.")
			elif sum(check_size) == 3: # all tables have data
				ons_tmp = ons1.merge(
						ons2, 
						how='outer'
						)
				ons_tmp = ons_tmp.merge(
					rd_ons,
					how='outer'
					)
				ons= (ons_tmp
					.sort_values(
						by='date_of_death', 
						ascending=True,
						na_position='last')
					.drop_duplicates(
						subset='participant_id', 
						keep='first')
				)
				tmp_mort = pd.merge(pid, ons, on='participant_id', how='left')
				tmp_mort['status'] = np.where(
					tmp_mort['date_of_death'].isna(), 'Alive', 'Deceased')

				self.mortality_table[key] = tmp_mort
				self.feature_tables['mortality'] = self.mortality_table
			elif sum(check_size) == 2:  # only one of the two lists has got data.
				res = [[ons1, ons2, rd_ons][i] for i, val in enumerate(check_size) if val]
				res_merge = pd.concat(res)
				
				res_merge_clean = (res_merge
					.sort_values(
						by='date_of_death', 
						ascending=True,
						na_position='last')
					.drop_duplicates(
						subset='participant_id', 
						keep='first')
					)
				tmp_mort = pd.merge(
					pid, 
					res_merge_clean, 
					on='participant_id', 
					how='left')
				tmp_mort['status'] = np.where(
					tmp_mort['date_of_death'].isna(), 'Alive', 'Deceased')
				self.mortality_table[key] = tmp_mort 
				self.feature_tables['mortality'] = self.mortality_table
			
			elif sum(check_size) == 1:  # only one of the two lists has got data.
				res = [[ons1, ons2, rd_ons][i] for i, val in enumerate(check_size) if val]
				tmp_mort = pd.merge(pid, res[0], on='participant_id', how='left')
				tmp_mort['status'] = np.where(
					tmp_mort['date_of_death'].isna(), 'Alive', 'Deceased')
				self.mortality_table[key] = tmp_mort  # unlisting the datarame.
				self.feature_tables['mortality'] = self.mortality_table
	

	########### sample level data ###################
	def rd_sample_data(self):
		self.rd_samples = {}

		if self.platekeys is not None:

			for key, pid in self.pids.items():
				rd_samp_query = (f'''
					SELECT
						participant_id,
						plate_key,
						family_id,
						assembly,
						is_proband,
						normalised_specific_disease_proband,
						affection_status,
						platypus_vcf_path,
						case_solved_family,
						last_status
					FROM
						rare_disease_interpreted
					WHERE
						plate_key IN {*self.platekeys,}
					''')
				rd_samp = lab_to_df(
					sql_query=rd_samp_query,
					dr=self.version
				)
				self.rd_samples[key] = rd_samp	
			self.sample_tables['rare_disease_samples'] =  self.rd_samples
		else:
			# for each participants retrieve the samples from cancer_analysis
			for key, pid in self.pids.items():
				rd_samp_query = (f'''
				SELECT
						participant_id,
						plate_key,
						family_id,
						assembly,
						is_proband,
						normalised_specific_disease_proband,
						affection_status,
						platypus_vcf_path,
						case_solved_family,
						last_status
					FROM
						rare_disease_interpreted
					
					WHERE
						participant_id IN {*pid,}
					''')
				rd_samp = lab_to_df(
					sql_query=rd_samp_query,
					dr=self.version
				)
				self.rd_samples[key] = rd_samp
			# self.all_cancer_samples = self.concat_cohort(self.cancer_samples)
			# concat before adding to sample_tables?
			self.sample_tables['rare_disease_samples'] =  self.rd_samples
			self.platekeys=self.concat_cohort(self.rd_samples)['plate_key']


	def ca_sample_data(self):
		"""query labkey for cancer sample information. Data is retrieved from 
		cancer analysis (for nsv4 sample stats) and supplemented with dragen 3.2
		realigned samples.
		We make a distinction between cohorts generated from platekeys and
		cohorts generated from participant_ids / disease terms. As we want to limit
		the output to specific platekeys in the former, but may want to grab all 
		samples per participant_id for the latter.
		(there are multiple tumour samples available)

		Args:
			platekeys (pd.Array or list): list of tumour_sample_platekeys to include.
			self.pids (dict of list): list of tumour_sample_platekeys to include.
			version (str): Data release version

		Returns:
			mortality_table (pd.DataFrame): date of death per participant_id
		"""
		self.cancer_samples = {}

		if self.platekeys is not None:

			for key, pid in self.pids.items():
				ca_samp_query = (f'''
					SELECT
						ca.participant_id,
						ca.tumour_sample_platekey,
						ca.germline_sample_platekey,
						ca.tumour_type,
						ca.disease_type,
						ca.disease_sub_type,
						ca.study_abbreviation,
						ca.morphology_icd,
						ca.histology_coded,
						ca.tumour_delivery_date,
						ca.germline_delivery_date,
						ca.preparation_method,
						ca.tumour_purity,
						ca.coverage_homogeneity,
						ca.somatic_coding_variants_per_mb AS tmb,
						ca.signature_1,
						ca.signature_2,
						ca.signature_3,
						ca.signature_4,
						ca.signature_5,
						ca.signature_6,
						ca.signature_7,
						ca.signature_8,
						ca.signature_9,
						ca.signature_10,
						ca.signature_11,
						ca.signature_12,
						ca.signature_13,
						ca.signature_14,
						ca.signature_15,
						ca.signature_16,
						ca.signature_17,
						ca.signature_18,
						ca.signature_19,
						ca.signature_20,
						ca.signature_21,
						ca.signature_22,
						ca.signature_23,
						ca.signature_24,
						ca.signature_25,
						ca.signature_26,
						ca.signature_27,
						ca.signature_28,
						ca.signature_29,
						ca.signature_30,
						csc.component_tnm_t,
						csc.component_tnm_n,
						csc.component_tnm_m,
						csc.final_figo_stage,
						csc.stage_best,
						csc.grade,
						csc.sact_tumour_pseudo_id,
						csc.er_status,
						csc.pr_status,
						csc.her2_status,
						csc.npi,
						csc.gleason_primary,
						csc.gleason_combined,
						csc.behaviour_coded_desc,
						csc.histology_coded_desc,
						csc.diagnosisdatebest AS diagnosis_date_best,
						ca.somatic_small_variants_annotation_vcf AS nsv4_somatic_small_variants_annotation_vcf,
						da.somatic_small_variants_annotation_vcf AS dragen_somatic_small_variants_annotation_vcf,
						ca.tumour_sv_vcf AS nsv4_somatic_sv_vcf,
						da.somatic_sv_vcf AS dragen_somatic_sv_vcf,
						da.somatic_cnv_vcf AS dragen_somatic_cnv_vcf
					FROM
						cancer_analysis ca
					LEFT JOIN
						cancer_100K_genomes_realigned_on_pipeline_2 da
					ON
						ca.tumour_sample_platekey = da.tumour_sample_platekey
					LEFT JOIN
						cancer_staging_consolidated csc
					ON
						ca.tumour_sample_platekey = csc.tumour_sample_platekey
					WHERE
						ca.tumour_sample_platekey IN {*self.platekeys,}
					''')
				ca_samp = lab_to_df(
					sql_query=ca_samp_query,
					dr=self.version
				)
				self.cancer_samples[key] = ca_samp	
			# cancer_samples to concat in concat_all
			# self.all_cancer_samples = self.concat_cohort(self.cancer_samples)
			self.sample_tables['cancer_samples'] =  self.cancer_samples
		else:
			# for each participants retrieve the samples from cancer_analysis
			for key, pid in self.pids.items():
				ca_samp_query = (f'''
					SELECT
						ca.participant_id,
						ca.tumour_sample_platekey,
						ca.germline_sample_platekey,
						ca.tumour_type,
						ca.disease_type,
						ca.disease_sub_type,
						ca.study_abbreviation,
						ca.morphology_icd,
						ca.histology_coded,
						ca.tumour_delivery_date,
						ca.germline_delivery_date,
						ca.preparation_method,
						ca.tumour_purity,
						ca.coverage_homogeneity,
						ca.somatic_coding_variants_per_mb AS tmb,
						ca.signature_1,
						ca.signature_2,
						ca.signature_3,
						ca.signature_4,
						ca.signature_5,
						ca.signature_6,
						ca.signature_7,
						ca.signature_8,
						ca.signature_9,
						ca.signature_10,
						ca.signature_11,
						ca.signature_12,
						ca.signature_13,
						ca.signature_14,
						ca.signature_15,
						ca.signature_16,
						ca.signature_17,
						ca.signature_18,
						ca.signature_19,
						ca.signature_20,
						ca.signature_21,
						ca.signature_22,
						ca.signature_23,
						ca.signature_24,
						ca.signature_25,
						ca.signature_26,
						ca.signature_27,
						ca.signature_28,
						ca.signature_29,
						ca.signature_30,
						csc.component_tnm_t,
						csc.component_tnm_n,
						csc.component_tnm_m,
						csc.final_figo_stage,
						csc.stage_best,
						csc.grade,
						csc.sact_tumour_pseudo_id,
						csc.er_status,
						csc.pr_status,
						csc.her2_status,
						csc.npi,
						csc.gleason_primary,
						csc.gleason_combined,
						csc.behaviour_coded_desc,
						csc.histology_coded_desc,
						csc.diagnosisdatebest AS diagnosis_date_best,
						ca.somatic_small_variants_annotation_vcf AS nsv4_somatic_small_variants_annotation_vcf,
						da.somatic_small_variants_annotation_vcf AS dragen_somatic_small_variants_annotation_vcf,
						ca.tumour_sv_vcf AS nsv4_somatic_sv_vcf,
						da.somatic_sv_vcf AS dragen_somatic_sv_vcf,
						da.somatic_cnv_vcf AS dragen_somatic_cnv_vcf
					FROM
						cancer_analysis ca
					LEFT JOIN
						cancer_100K_genomes_realigned_on_pipeline_2 da
					ON
						ca.tumour_sample_platekey = da.tumour_sample_platekey
					LEFT JOIN
						cancer_staging_consolidated csc
					ON
						ca.tumour_sample_platekey = csc.tumour_sample_platekey
					WHERE
						ca.participant_id IN {*pid,}
					''')
				ca_samp = lab_to_df(
					sql_query=ca_samp_query,
					dr=self.version
				)
				self.cancer_samples[key] = ca_samp
			# self.all_cancer_samples = self.concat_cohort(self.cancer_samples)
			# concat before adding to sample_tables?
			self.sample_tables['cancer_samples'] =  self.cancer_samples
			self.platekeys=self.concat_cohort(self.cancer_samples)['tumour_sample_platekey']


	def omics_sample_metadata(self):
		"""Extract sample metadata from labkey. Specifically looking at omics
		sample availability and location. Please keep in mind the sample counts
		are subject to change.
		
		Returns: 
			omics_sample_data (dictionary): The keys of which correspond to the
			different cohorts in self.pids. This table contains info on each sample
			collected.
			omics_sample_counts: the number of samples per sample_type 
			(RNA)
			) 
		"""
		self.omics_sample_data = {}
		self.omics_sample_counts = {}
		self.omics_sample_location = {}
		self.gmc_registration = {}
		for key, pid in self.pids.items():
			meta_query = (f'''
				SELECT
					samp.participant_id,
					samp.laboratory_sample_id,
					samp.laboratory_sample_gmc_ods_code,
					samp.laboratory_sample_gmc_trust,
					samp.sample_type,
					samp.sample_preparation_method,
					samp.tumour_germline_omics,
					samp.sample_received_at_biorepository,
					omics.aliquots
				FROM
					laboratory_sample samp
				LEFT JOIN 
					laboratory_sample_omics_availability omics
				ON 
					samp.laboratory_sample_id = omics.laboratory_sample_id
				WHERE
					samp.participant_id IN {*pid,}
				''')
			sample_meta = lab_to_df(
				sql_query=meta_query,
				dr=self.version
				)
			# create a histogram of the number of aliquots per sample type
			counts = sample_meta.groupby(
				['sample_type', 'aliquots']
				).size()

			location = sample_meta.groupby(
				[
					'sample_type',
					'laboratory_sample_gmc_trust',
					'laboratory_sample_gmc_ods_code'
					]
				).size()
			location.reset_index(drop=False, name='count')
			
			self.omics_sample_counts[key] = counts
			self.omics_sample_location[key] = location
			self.omics_sample_data[key] = sample_meta
		
		# should this be a seperate function/call? its not only sample related.
		# retrieve GMC registration site per participant
		for key, pid in self.pids.items():
			gmc_query = (f'''
				SELECT
					participant_id, 
					registered_at_gmc_trust,
					registered_at_gmc_ods_code
				FROM
					participant
				WHERE
					participant_id IN {*pid,} '''
					)
			registration_loc = lab_to_df(
				sql_query=gmc_query,
				dr=self.version
				)
			self.gmc_registration[key] = registration_loc

		self.feature_tables['gmc_registration'] = self.gmc_registration
		self.sample_tables['omics_sample_data'] = self.omics_sample_data

		
			
	################ summarysing the participants ##################
	def concat_cohort(cls, dc):
		"""concatenate the dataframes in a dictionary, note; the dicitonarys
		should have the same structure / columns. 

		Args:
			dc (dictionary): Dictionary value's are pd.DataFrames, one per Key.

		Returns:
			pd.DataFrame: A concatenated pandas dataframe with unique values in
			the dictionary.
		"""
		if len(dc.keys()) > 1:
			concatdc = (pd.concat(dc.values())
				.reset_index(drop=True)
				.drop_duplicates()
				)
		else:
			concatdc = next(iter(
				dc.values()
				)).drop_duplicates()

		return concatdc
	

	def summarize(cls, pid=None, anc=None, mort=None, age=None, sex=None):
		"""gather up self data and return summary stats. including:
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
		# not all features have to be calculated/present.
		# only include those that are.
		# instead of building a dataframe in one go we'll concat each
		# column that we actually have.
		summary_series = []  # add each feature we want to include as named series.
		if pid is not None:
			summary_series.append(
				pd.Series(len(pid), name='n_unique_participants')
			)
		if age is not None:
			summary_series.append(
				pd.Series(
					np.mean(age['age_at_consent']), 
					name='mean_consent_age'
					)
					)
			summary_series.append(
				pd.Series(
					np.std(age['age_at_consent']), 
					name='std_consent_age'
					)
					)
			summary_series.append(
				pd.Series(
					np.mean(age['current_age']), 
					name='mean_current_age'
					)
					)
			summary_series.append(
				pd.Series(
					np.std(age['current_age']), 
					name='std_current_age'
					)
					)
		if sex is not None:
			summary_series.append(
				pd.Series(
					sex['sex'].value_counts(normalize=True)[0], 
					name='fraction_female'
					)
					)
		if mort is not None:
			summary_series.append(
				pd.Series(
					mort['status'].value_counts(normalize=True)[1], 
					name='fraction_deceased'
					)
					)

		summary = pd.concat(summary_series, axis=1)

		# ancestry_distribution is already 'wide' - various columns.	
		if anc is not None:
			ancestry_distribution = (anc.predicted_ancestry
				.value_counts(normalize=True)
				.reset_index(name='0', drop=False)
				.T)
			ancestry_distribution.columns = ancestry_distribution.iloc[0]
			ancestry_distribution.drop(ancestry_distribution.index[0], inplace=True)
			

		# summary = pd.DataFrame(
		# 	{
		# 	'n_unique_participants':len(pid),
		# 	'mean_consent_age':np.mean(age['age_at_consent']),
		# 	'std_consent_age':np.std(age['age_at_consent']),
		# 	'mean_current_age':np.mean(age['current_age']),
		# 	'std_current_age':np.std(age['current_age']),
		# 	'fraction_female':
		# 		sex['sex'].value_counts(normalize=True)[0],
		# 		  # alphabetically ordered.
		# 	'fraction_deceased': 
		# 		mort['status'].value_counts(normalize=True)[1]
		# 	}, index=[0])
		
		summary_anc = pd.concat(
			[
				summary, 
				ancestry_distribution.reset_index(drop=True)
				], 
			axis=1)

		return summary_anc
	

	def concat_all(self):
		"""concatenates all the features and sources in the cohort to a long
		type table, matching each patient/diagnosis combination with its features
		this is usefull for filtering, especially when you want to exclude
		participants of a certain feature (e.g. age) and diagnosis - 
		but don't want to blanket remove based on features. (e.g. 
		filter participants age 30 and under for diagnosis A, but include all ages
		for diagnosis B: a participant who is 25 yo and was diagnosed for both
		A and B will be included in the cohort)

		Raises:
			RuntimeError: When no valid cohort is present
			RuntimeError: when no valid features are present
		"""
		# this function requires some restructuring, 
		# right now it checks if there are participant_ids,
		# if there are features, and if there are samples.
		# if either one of those are not present it will return a runtime error
		# probably better to:
		# 
		# 1. only return a runtime error if there are no participant ids
		diag_keys = [
			'icd10',
			'dterm',
			'cterm',
			'cabbr',
			'hpo',
			'custom'
		]
		
		if all([not k in self.pids.keys() for k in diag_keys]):
			raise RuntimeError('No valid cohorts to concatenate found')
		else:
			# 2. concat the different participantids in self.all_pids.
			self.all_pids = self.concat_cohort(self.pids)
			tmp_merge = [self.all_pids]
			# 3. concat the diagnosis fields.
			conc_tables = []
			for key in diag_keys:
				if key in self.pids.keys():
					if key == 'dterm':
						conc_tables.append(
							self.dterm_table.rename(
								{'normalised_specific_disease':'diag'}, 
								axis=1)
								)
					if key == 'icd10':
						conc_tables.append(
							self.icd10_table.rename(
								{'code':'diag'}, 
								axis=1)
								)	
					if key == 'cterm':
						conc_tables.append(
							self.cterm_table.rename(
								{'disease_type':'diag'},
								axis=1)
								)
					if key == 'cabbr':
						conc_tables.append(
							self.cabbr_table.rename(
								{'study_abbreviation':'diag'},
								axis=1)
								)
					if key == 'hpo':
						conc_tables.append(
							self.hpo_table.rename(
								{'normalised_hpo_id':'diag'},
								axis=1
							)
						)
					if key == 'custom':
						tmpframe = pd.DataFrame({
							'participant_id':self.pids['custom']
							})
						tmpframe['diag'] = 'custom'
						conc_tables.append(tmpframe)

			ntable = sum([k in self.pids.keys() for k in diag_keys])
			if ntable > 1:
				diag_table = pd.concat(conc_tables)
			else:
				diag_table = conc_tables[0]


			# 4. see if there are features (age,sex,mortality),
			# if there are, include them in the merge.

			# feature_tables is a method for tracking what
			# features we've calculated in the cohort
			# feature_tables = {
			# 	'age':self.age_table,
			# 	'ancestry':self.ancestry_table,
			# 	'sex':self.sex_table,
			# 	'mortality':self.mortality_table
			# }
			# slightly unusual double negative here,
			# if none of the feature tables exist in the class
			# we want to send a warning, otherwise continue.
			# alternatively we could:
			# if any ([k for k in feature tables:]) 
			# and then call the error at the end of the for loop.
			if all([not feature for name, feature in self.feature_tables.items()]):
				warnings.warn('Currently no features to concatenate found.')
			else:
				# we are concatenating all the potential diagnosis, hpo, icd10, disease terms
				# as the researcher may want to apply filtering to specific parts of the self.
				# e.g. filter out participants under age 30 with ICD-10 code A. But keep 
				# all participants (irrelevant of age) for a particular hpo term. If a participant
				# is included for the icd-10 code AND the hpo term we don't want to completely
				# remove the participant_id from self.pids.
				for name, feature in self.feature_tables.items():
					if feature:
						if name == 'age':
							self.all_age = self.concat_cohort(feature)
							tmp_merge.append(self.all_age)
						elif name == 'ancestry':
							self.all_ancestry = self.concat_cohort(feature)
							tmp_merge.append(self.all_ancestry)
						elif name == 'sex':
							self.all_sex = self.concat_cohort(feature)
							tmp_merge.append(self.all_sex)
						elif name == 'mortality':
							self.all_mortality = self.concat_cohort(feature)
							tmp_merge.append(self.all_mortality)
						elif name == 'gmc_registration':
							self.all_gmc_registration = self.concat_cohort(feature)
							tmp_merge.append(self.all_gmc_registration)

			# 5. see if there are samples (ca, rd, other),
			# include them in the merge.
			if all([not samples for name, samples in self.sample_tables.items()]):
				warnings.warn('Currently no sample data to concatenate found.')
			else:
				for name, samples in self.sample_tables.items():
					if name == 'omics_sample_data':
						self.all_omics_metadata = self.concat_cohort(samples)
						tmp_merge.append(self.all_omics_metadata)
					elif name == 'cancer_samples':
						self.all_cancer_samples = self.concat_cohort(samples)
						tmp_merge.append(self.all_cancer_samples)
					elif name == 'rare_disease_samples':
						self.all_rd_samples = self.concat_cohort(samples)
						tmp_merge.append(self.all_rd_samples)
						# since the key here are the samples and not the 
						# participants we shouldn't merge them to all_data.
						# where the key is diag+participant_id

		# determine which of the features are available in the class
		# then perform stepwise merges for those
		from functools import reduce
		df_merged = reduce(lambda left,right: 
			pd.merge(
				left, right,
				on=['participant_id'],
				how='left'), tmp_merge)

		self.all_data = pd.merge(
			diag_table,
			df_merged,
			on='participant_id',
			how='left'
			)

		
	def summary_stats(self):
		"""calculate feature stats of the cohort, this should ideally be all
		features (age, sex, ancestry, mortality). depends on the summarize() 
		function. First collapse all the participant ids in the featuers (as they
		are generated per diag/ontology, e.g. self.age has got a table for
		icd10, hpo part of the cohort.)

		Returns:
			self.summary: a summary table of the entire cohort. not per ontology.
		"""

		# setting tables to none to allow summary() to only use those which we
		# have calculated.
		all_ancestry = all_mortality = all_age = all_sex = all_omics = all_reg = None
		# collapse all the features we have identified
		all_pids = self.concat_cohort(self.pids)
		for name, feature in self.feature_tables.items():
			if name == 'age':
				all_age = self.concat_cohort(feature)
			elif name == 'ancestry':
				all_ancestry = self.concat_cohort(feature)
			elif name == 'sex':
				all_sex = self.concat_cohort(feature)
			elif name == 'mortality':
				all_mortality = self.concat_cohort(feature)
			elif name == 'omics_sample_data':
				all_omics = self.concat_cohort(feature)
			elif name == 'gmc_registration':
				all_reg = self.concat_cohort(feature)

			else:
				return KeyError('Feature in feature table not recognised.')


		self.summary = self.summarize(
			pid=all_pids,
			anc=all_ancestry,
			mort=all_mortality,
			age=all_age,
			sex=all_sex,
			)


	def load_icd10_data(self):
		"""gelpack contains a translation file for icd-10 codes, this function
		loads in the relevant file - used for generating summary statitics per
		icd-10 code.
		"""
		from importlib import resources
		with resources.path("gelpack",'coding19.tsv') as f:
			self.icd10_lookup = pd.read_csv(
				f, 
				sep='\t',
				usecols=['coding','meaning'],
				).rename({'coding':'code'}, axis=1)


	def ontology_stats(self):
		"""calculate cohort summary per ontology/source (hpo, icd10, cterm,
		dterm, custom).

		Raises:
			RuntimeError: Returns an error when there are no ontologies to
			summarize.

		Returns:
			self.ont_vcount and self.icd10_overlap_matrix
		"""
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
		# count_dict = {} is already part of summary_ont now 
		# (n_unique_participants).
		# vcount should summarize the counts for each resource.		
		self.ont_vcount = {}
		if all([not k in self.pids.keys() for k in [
			'dterm',
			'hpo',
			'cterm',
			'cabbr',
			'icd10',
			'custom']]):   # TODO add custom vcounts.
			raise RuntimeError('No valid cohorts to summarize found')
		else:
			# check which keys have been included in the self.
			for key in ['dterm', 'hpo', 'cterm', 'cabbr', 'icd10', 'custom']:
				if key in self.pids.keys():
					if key == 'custom':
						tmpframe = pd.DataFrame({
							'participant_id':self.pids['custom']
							})
						tmpframe['diag'] = 'custom'
						self.ont_vcount[key] = (tmpframe 
							.drop_duplicates([
								'participant_id',
								'diag'])
							.diag
							.value_counts()
						)  
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
					# cancer self?
					# TODO: add a check for icd-10 codes (confirm diagnosis)
					# TODO: add a ICD-10 -> sample / option.
					elif key == 'cterm':  
						self.ont_vcount[key] = (self.cterm_table
							.drop_duplicates([
								'participant_id', 
								'disease_type'
								])
							.disease_type
							.value_counts()
						)
					elif key == 'cabbr':  
						self.ont_vcount[key] = (self.cabbr_table
							.drop_duplicates([
								'participant_id', 
								'study_abbreviation'
								])
							.study_abbreviation
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


	############### filtering the cohort ####################
	def feature_filter(self, query_string, action='exclude'):
		try:
			filter_cohort = self.all_data.query(query_string, engine='python')
		except AttributeError:
			self.concat_all()
			filter_cohort = self.all_data.query(query_string, engine='python')
		
		# note that we implicitely do not remove participants when they may
		# remain in the cohort under different diagnosis codes.
		# TODO: this may cause confusion and needs to be reviewed.
		if action =='exclude':
			# if we drop all these we are inadvertantly going to drop participants who should be included.
			# we could drop the from a copy of the table and check for each participant id
			# if they remain in another query.
			main = pd.concat([self.all_data.copy(), filter_cohort])
			main.drop_duplicates(keep=False)
			# now check for filter_cohort which ones are not in main.participant_ids
			filter_pids = [x for x in filter_cohort['participant_id'] if x not in main.participant_id]
			
		elif action =='include':
			# not this 'include' is different from the costuom_pids() action == include
			# custom_pids include just adds the participant_ids to self.pids['custom']
			# if we want to filter our data to only `include` the selected participants
			# we have to exclude the inverse.
			filter_pids = [
				x for x in self.all_data['participant_id'] 
					if x not in filter_cohort['participant_id']
					]
		else:
			raise NameError('no valid argument for `action` included.')

		# subsequently remove the participant ids from the cohort with
		self.custom_pids(
			pids_lst_file=filter_pids,
			action='exclude'
			)
		# refresh trhe all_feature.
		self.concat_all()


	def select_single_ca_sample(self):
		"""This function attempts to select a single cancer sample from
		participants with multiple samples. In those cases it first removes
		samples that are not in cancer_analysis. If there are still remaining
		samples it will remove non-Primary tumour samples. Finally, if there 
		are still multiple primary tumour samples it will select one with the 
		highest tumour_purity, and highest coverage_homogeneity in case of ties.
		"""
		self.concat_all()
		# identify participants with more than 1 sample:
		subcounts = (self
			.all_cancer_samples
			.participant_id
			.value_counts()
			.reset_index(name='count')
			)
		dups = subcounts.loc[subcounts['count']>1, 'index']
		to_drop = []
		for pid in dups:
			tmp_drop = []
			samps = (self
				.all_cancer_samples
				.loc[self.all_cancer_samples['participant_id'] == pid]
				)

			# select a platekey, (first remove those not in cancer_analysis.)
			no_interp = samps.loc[
					samps['nsv4_somatic_small_variants_annotation_vcf'].isna(),
					'tumour_sample_platekey'
					].tolist()

			tmp_drop += no_interp
			
			if (len(tmp_drop) < len(samps)-1): # this means we keep non-primary tumours if no other samples are available.
				# remove non-primary
				non_primary = samps.loc[samps['tumour_type']!='PRIMARY',
					'tumour_sample_platekey'].tolist()
				tmp_drop += non_primary

			if (len(tmp_drop) < len(samps)-1):
				# then grab highest tumour_purity and coverage_homogeneity.
				pur_cov_drop = (samps
					.sort_values([
						'tumour_purity',
						'coverage_homogeneity'
						], 
						axis=0,
						ascending=[True,True]) 
					).tumour_sample_platekey.iloc[:-1].tolist()
				tmp_drop += pur_cov_drop
			
			to_drop += tmp_drop

		# this should leave one platekey per participant_id.
		# the other platekey should be removed from:
		# cohort.platekeys:
		self.platekeys = self.platekeys[~self.platekeys.isin(to_drop)]
		# cohort.sample_tables()
		for key,table in self.cancer_samples.items():
			if key == 'custom':
				self.cancer_samples['custom'] = table.loc[
					~table['tumour_sample_platekey'].isin(to_drop)
					]
			elif key == 'icd10':
				self.cancer_samples['icd10'] = table.loc[
					~table['tumour_sample_platekey'].isin(to_drop)
					]
			elif key == 'cterm':
				self.cancer_samples['cterm'] = table.loc[
					~table['tumour_sample_platekey'].isin(to_drop)
					]
			elif key == 'dterm':
				self.cancer_samples['dterm'] = table.loc[
					~table['tumour_sample_platekey'].isin(to_drop)
					]
			elif key == 'cabbr':
				self.cancer_samples['cabbr'] = table.loc[
					~table['tumour_sample_platekey'].isin(to_drop)
					]
			elif key == 'hpo':
				self.cancer_samples['hpo'] = table.loc[
					~table['tumour_sample_platekey'].isin(to_drop)
					]
		self.sample_tables['cancer_samples'] = self.cancer_samples
		# if the sample_tables have been filtered
		# all_cancer_samples gets reset with concat_all():
		self.concat_all()		
