#!/usr/bin/env python

#### Christian J. Bouwens
#### BRS team
#### generate Kaplann-meier survival curves on domain variants and
#### disease types.
#### based on Mariana Pareira's survival.R
#### last update: 2022.10.18

#########################
'''TODO
- create confluence page with documentation
- Check if part 2 of the script works. 
with a number of participant ids, return the survival time and / or plot
- check if part 3 works: with a list of survival times / censoring, plot.
- add coxPH.
'''
import argparse
import numpy as np
import pandas as pd
import warnings
import labkey
import re
from functools import reduce


# this version of pandas is still suffering from the _str_len bug.
pd.set_option("display.max_columns", None)
#################
# functions
#################

def argparser():
	"""Parse argument from bash script to python functions.

	Returns:
		options: arguments per flag in the format options.flag
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument('-g','--genes',
		dest='genes',
		nargs='+',
		default=None,
		required=True,
		help='<Required> Array of Domain variants to stratify'
			' survival analysis on.')
	parser.add_argument('-s', "--stratification",
		dest='strata',
		default='and',
		help='how should multi-genes stratification occur.'
			'["and", "or", "full"]: '
			'and = all genes require to be mutated to be in group.'
			'or = any gene requires to be mutated to be in group.'
			'full = all possible combinations of muts/wt are assigned groups.'
			'Default = "and".')
	parser.add_argument('-i', "--imputate",
		dest='imputate_flag',
		action='store_true',
		default=True,
		help='Flag: should missing diagnosis dates be imputated?'
			'Default=True')
	parser.add_argument('-v', "--version",
		dest='version',
		default="main-programme/main-programme_v16_2022-10-13",
		help='GEL data release version',)
	parser.add_argument('-o', "--out",
		dest='out',
		default=None,
		required=True,
		help='Full path of output directory.')
	parser.add_argument('-d', "--disease",
		dest='disease_type',
		nargs='+',
		default=None,
		help=('Cancer disease_types, or clinical indications to include, '
			'if none given all GEL disease_types will be included'))
	parser.add_argument('-t', "--title",
		dest='plttitle',
		default='Survival analysis',
		required=False,
		help='Title for Kaplan-Meier plot.')

	options=parser.parse_args()
	
	return options


def lab_to_df(sql_query, dr):
	"""generate an pandas dataframe from labkey sql query

	Args:
		sql_query (str): an sql query as string.
		dr (str): GEL datarelease version
	"""
	server_context = labkey.utils.create_server_context(
		domain= "labkey-embassy.gel.zone",
		container_path = dr,
		context_path = "labkey",
		use_ssl = True
	)
	
	results =  labkey.query.execute_sql(
		server_context,
		schema_name="lists",
		sql=sql_query,
		max_rows=5000000
	)
	return(pd.DataFrame(results['rows']))


ic_lookup = [  # the ic_lookup regex doesn't function the same in R and python.
	('ADULT_GLIOMA',r"C71[0-9]{0,1}|D43[0-4]{0,1}|D32[0-4]{0,1}|D33[0-4]{0,1}"),
	('BLADDER',r"C67[0-9]{0,1}|D090|D414"),
	('BREAST',r"C50[0-9]{0,1}|D05[0-9]{0,1}"),
	('CARCINOMA_OF_UNKNOWN_PRIMARY',r"C80[0-9]{0,1}|D489"),
	# ('OTHER',r"C80[0-9]{0,1}|D489"),  # putting unknown in OTHER, like childhood.
	('COLORECTAL',r"C18[0-9]{0,1}|C20|C19|D010|D012|D011|C17[0-9]{0,1}"
		r"|C21[0-9]{0,1}|C78$|C784|C785|C788"),  # changed from 
	('ENDOCRINE',r"C73|C74[0-9]{0,1}|C75[0-9]{0,1}|D093|D44[0-9]{0,1}"),
	('ENDOMETRIAL_CARCINOMA',r"C53[0-9]{0,1}|C54[0-9]{0,1}|C55[0-9]{0,1}|D070|D069"),
	('HAEMONC',r"C81[0-9]{0,1}|C82[0-9]{0,1}|C83[0-9]{0,1}|C84[0-9]{0,1}"
		r"|C85[0-9]{0,1}|C86[0-9]{0,1}|C88[0-9]{0,1}|C90[0-9]{0,1}|C91[0-9]{0,1}"
		r"|C92[0-9]{0,1}|C93[0-9]{0,1}|C94[0-9]{0,1}|C95[0-9]{0,1}|C96[0-9]{0,1}|"
		r"D47[0-9]{0,1}"),
	('HEPATOPANCREATOBILIARY',r"C22[0-9]{0,1}|C23|C24[0-9]{0,1}"
		r"|C25[0-9]{0,1}|C787|D015|D376"),
	('LUNG',r"C34[0-9]{0,1}|D022|C450|C459"),
	('MALIGNANT_MELANOMA',r"C43[0-9]{0,1}|D03[0-9]{0,1}"),
	('NASOPHARYNGEAL',r"C11[0-9]{0,1}"),
	('ORAL_OROPHARYNGEAL',r"C01[0-9]{0,1}|C02[0-9]{0,1}|C03[0-9]{0,1}|C04[0-9]"
		r"{0,1}|C05[0-9]{0,1}|C06[0-9]{0,1}|C07[0-9]{0,1}|C08[0-9]{0,1}"
		r"|C09[0-9]{0,1}|C10[0-9]{0,1}|C00[0-9]{0,1}|C14[0-9]{0,1}|D000|D370"),
	('OVARIAN',r"C56|C796|D391|C57[0-9]{0,1}"),
	('PROSTATE',r"C61|D075"),
	('RENAL',r"C64|C65|C66|C68[0-9]{0,1}|C790|D410|D411|D412"),
	('SARCOMA',r"C41[0-9]{0,1}|C45\b|C45[1-8]{1}|C46[0-9]{0,1}|C48[0-9]{0,1}"
		r"|C49[0-9]{0,1}|C40[0-9]{0,1}|C786|D481|D483|D484"),
	('SINONASAL',r"C12|C30[0-9]{0,1}|C31[0-9]{0,1}"),
	('TESTICULAR_GERM_CELL_TUMOURS',r"C62[0-9]{0,1}"),
	('UPPER_GASTROINTESTINAL',r"C15[0-9]{0,1}|D001|C16[0-9]{0,1}"
		r"|D002|C26[0-9]{0,1}|D017|D019|D379")
]


def translateicd(icd_vec, lookups=ic_lookup):
	"""As we are importing ICD-10 codes from different RWE sources their 
	formatting has to be unified. this function cleans the codes and returns
	a matching disease type. 
	

	Args:
		icd_vec (list): list of str, ICD-10 codes from different sources.
		icd_dictionary (DataFrame): reference dictionary of disease types
			and regex-form ICD-10 codes
	"""
	import re

	def transicd(icd, lookups=ic_lookup):
		for value, pattern in lookups:
			if re.search(pattern, icd):
				return value
		return 'OTHER'
	# remove trailing 'X' and any '.'
	icd_vec_clean = [re.sub('X$|\[.\]|[.]' , '', x) for x in icd_vec]
	return list(map(transicd, icd_vec_clean))


def assign_groups(dataframe, vars, type='and'):
	"""assigns groups/labels for survival analysis based on a list of variables.


	Args:
		dataframe (pd.DataFrame): A pandas dataframe with samples / participants 
			their survival time, censoring and variables to stratify on.
		vars (list): list of variable names in string format, corresponding to 
			columns of boolean series (if a var is True or False in a sample)
			in dataframe. 
		type (str, optional): determines how the groups are assigned:
			'and': binary groups (1, 2) where var1 AND var2 AND..var* are True.
			'or': binary groups where var1 OR var2 OR var* are True.
			'full': groups on each possible iteration of vars. 
			Defaults to 'and'.
	Returns:
		pd.series of grouping assigned as group_1, group_2, group_n
	"""	
	# check if each var is indeed in dataframe
	dataframe=dataframe.copy()
	for var in vars:
		if not var in dataframe:
			raise ValueError(
				f'{var} should only include columns in dataframe.'
				)
		if dataframe[var].isna().any():
			raise ValueError(
				f'{var} contains NA values, only True/False allowed'
				)

	if type=='and':
		dataframe['group'] = (dataframe[vars]
			.all(axis=1)
			.replace({True: 'group_1', False: 'group_2'})
		)
		# return dataframe('group')
	elif type=='or':
		dataframe['group'] = (dataframe[vars]
			.any(axis=1)
			.replace({True: 'group_1', False: 'group_2'})
		)
		# return dataframe['group']

	elif type=='full':
		dataframe['combination'] = dataframe[vars].agg(tuple, axis=1)
		dataframe['group'] = dataframe['combination'].factorize()[0]
		# dataframe.groupby(['snv1','snv2'], sort=False).ngroup()
		mapping = dataframe[['combination','group']].drop_duplicates()
		print(mapping)
		return dataframe, mapping
	else:
		raise ValueError(f'Unrecognized type: {type}.')

	return dataframe	


class Survdat(object):
	
	def __init__(self, df, pids, version, impute):
		import pandas as pd
		import warnings
		import labkey
		import re
		from functools import reduce
		
		self.ca = df
		self.pids = pids
		self.version = version
		self.impute = impute
	

	def quer_ons(self):
		"""Extract the death date from labkey tables for a set of participant_id.

		Args:
			pids (list): list of participant ids to include
			version (str): Data release version

		Returns:
			ons (pd.DataFrame): date of death per participant_id
		"""

		# note this f string SQL query is not infinitely scalable.
		# we'll run into the chacacter limit, but for 100K pids it is fine.
		###
		# determine dates_of_death 
		# functionize this, and let it be done on particular participant ids.
		# that way someone can identify participant ids, and merge it with their grouping.
		###
		# import date_of_death per participant_id from 'mortality'
		query2 = (
			f'''
			SELECT 
				DISTINCT participant_id, death_date
			FROM
				death_details
			WHERE participant_id IN {*self.pids,}
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
			WHERE participant_id IN {*self.pids,}
		''')
		ons2 = lab_to_df(
			sql_query=query3,
			dr=self.version)

		ons1.rename(columns={'death_date':'date_of_death'}, inplace=True)
		ons = ons1.merge(ons2, how='outer')
		ons.sort_values(by='date_of_death', ascending=True, inplace=True)
		ons.drop_duplicates(subset='participant_id', keep='first', inplace=True)

		self.ons = ons


	def quer_hes(self):
		"""Extract the date of last follow up from hospital episode statitsics
			labkey tables for a set of participant_id.

		Args:
			pids (list): list of participant ids to include
			version (str): Data release version

		Returns:
			hes (pd.DataFrame): date of last interaction per participant_id.
		"""
		###
		# determine date of last follow up (last time seem)
		###
		hes_tables=[
			('apc','disdate'),
			('ae','arrivaldate'),
			('op','apptdate'),
			('cc','ccdisdate')
			]
		# WHERE x.participant_id IN {*self.pids,}
		query = (f'''
		SELECT 
			x.participant_id, MAX(x.y) AS LastSeen 
		FROM 
			hes_x as x 
		GROUP BY x.participant_id 
		''')
		# replacing x and y with matched hes data and column name in.
		q_apc, q_ae, q_op, q_cc = [
			query.replace('x',i[0]).replace('y' ,i[1]) for i in hes_tables
			]


		apc_lastseen = lab_to_df(
			sql_query=q_apc,
			dr=self.version
			)

		ae_lastseen= lab_to_df(
			sql_query=q_ae,
			dr=self.version
			)

		op_lastseen= lab_to_df(
			sql_query=q_op,
			dr=self.version
			)

		cc_lastseen= lab_to_df(
			sql_query=q_cc,
			dr=self.version
			)

		dfs = [apc_lastseen, ae_lastseen, op_lastseen, cc_lastseen]

		# merge multiple dataframes to find the last entry in all.
		merged_lastseen = reduce(
			lambda left, right: pd.merge(
				left, right, on=['participant_id'],
				how = 'outer'),
			dfs)
		merged_lastseen.columns = ['participant_id', 'apc', 'ae', 'op', 'cc']

		# not all are datetimes, which we need in order to select last value.
		for i in [ 'apc', 'ae', 'op', 'cc']:
			merged_lastseen[i] = pd.to_datetime(merged_lastseen[i])
		merged_lastseen = merged_lastseen.loc[
			merged_lastseen['participant_id'].isin(self.pids)
			]

		hes=pd.DataFrame({
			'participant_id':merged_lastseen['participant_id'],
			'lastseen':merged_lastseen[['apc','ae','op','cc']].max(axis=1)})
				
		# hes=pd.DataFrame({
		# 	'participant_id':merged_lastseen['participant_id'],
		# 	'lastseen':merged_lastseen[['apc','ae','cc']].max(axis=1)})
		
		self.hes = hes


	def quer_dod(self):
		###
		# determine date of diagnosis
		##
		query1 =(f'''
			SELECT 
				DISTINCT participant_id, diagnosis_date, diagnosis_icd_code
			FROM
				cancer_participant_tumour
			WHERE 
				participant_id IN {*self.pids,}
		''')
		dod_participant = lab_to_df(
			sql_query=query1,
			dr=self.version
			)
		dod_participant['disease_type'] = translateicd(
			dod_participant['diagnosis_icd_code'],
			lookups=ic_lookup
			)
		dod_participant.sort_values(
			['diagnosis_date'],
			ascending=True,
			inplace=True)
		dod_participant.drop_duplicates(
			['participant_id', 'disease_type'],
			keep='first',
			inplace=True)

		query2=(f'''
			SELECT 
				DISTINCT participant_id, diagnosisdatebest, site_icd10_o2
			FROM
				av_tumour
			WHERE
				participant_id IN {*self.pids,}
		''')
		av_dod_participant = lab_to_df(
			sql_query=query2,
			dr=self.version
			)
		av_dod_participant['disease_type'] = translateicd(
			av_dod_participant['site_icd10_o2'],
			lookups=ic_lookup
			)
		av_dod_participant = (av_dod_participant
			.sort_values(
				['diagnosisdatebest'],
				ascending=True
				)
			.drop_duplicates(
				['participant_id', 'disease_type'],
				keep='first'
				)
			.rename(
				columns={'diagnosisdatebest':'diagnosis_date'}
				)
			)

		av_dod_participant_merged = (
			av_dod_participant[[
			'participant_id', 
			'diagnosis_date',
			'disease_type'
			]]
			.merge(dod_participant[[
				'disease_type',
				'participant_id',
				'diagnosis_date']],
				how='outer')
			.sort_values(['participant_id', 'diagnosis_date'], ascending=True)
			.drop_duplicates(['participant_id', 'disease_type'], keep='first')
			)

		self.dod = av_dod_participant_merged


	def merge_dod(self):
		###
		# merge diagnosis date with pID
		###
		### Merging diagnosis date with pID 
		### (by participant_id and then manually matching disease_type)
		### disease_type either matching
		# OR CHILDHOOD / CUP accept any TYPE 
		# OR if OTHER in translation accepted
		# this is neccessary as multiple diagnosis may have taken place
		# and we've grouped childhood/Carcinoma of unknown origin under other.
		pid_diag = pd.merge(
			self.ca,
			self.dod,
			how='left',
			on=['participant_id'])

		pid_diag['match'] = np.where(
			pid_diag['disease_type_x'] == pid_diag['disease_type_y'],
			1,0
			)
		self.pid_diag = (pid_diag
			.sort_values(
				['participant_id', 'disease_type_x', 'match'],
				ascending=[True, True, False]
				)
			.drop_duplicates(['participant_id', 'disease_type_x'])
			.rename(columns={'disease_type_x':'disease_type'})
			.loc[
				(pid_diag['disease_type_x'] == pid_diag['disease_type_y'])
				| (pid_diag['disease_type_x'].isin(
					['CHILDHOOD',
					'CARCINOMA_OF_UNKNOWN_PRIMARY'])
					)
				| (pid_diag['disease_type_y'] == 'OTHER')
			]
			.drop(['disease_type_y', 'match'],axis=1)
		)

		# which samples have no date of diagnosis?
		outer_join = self.ca.merge(
			self.pid_diag,
			how='left',
			on=['participant_id', 'disease_type'],
			indicator=True
			)
		no_diag = outer_join[~(outer_join._merge == 'both')]
		self.no_diag = no_diag.drop(['_merge'], axis=1)


	def dod_impute(self):
		"""fill in diagnosis_date for cancer participants whose diagnosis date is 
		unknown. Based on the average time between diagnosis and enrolment of 
		each cancer type.

		Args:
			miss_df (pd.DataFrame): Dataframe with 3 columns: 
			['disease_type', 'participant_id', 'diagnosis_date']
			Where 'diagnosis_date' is NA.

		Returns:
			pd.DataFrame: Dataframe with imputed diagnosis date filled in.
		"""

		query = (f'''
			SELECT
				DISTINCT participant_id, date_of_consent
			FROM 
				participant
			WHERE 
				CAST(participant_id AS Char) LIKE '2%'
			''')
		par = lab_to_df(query, dr=self.version)
		par = (par
			.sort_values(['participant_id', 'date_of_consent'], ascending=True)
			.drop_duplicates(['participant_id'], keep='first')
			)

		query = ('''
			SELECT	
				participant_id, disease_type
			FROM
				lists.cancer_analysis'''
		)
		fullca = lab_to_df(query, dr=self.version)
		
		################################################################
		# load in full av_dod / 
		# In quer_dod() we've limited it to particiants of interest.
		# now we want to generate an average time between enrolment and
		# diagnosis based on ALL participants
		################################################################
		query1 =(f'''
		SELECT 
			DISTINCT participant_id, diagnosis_date, diagnosis_icd_code
		FROM
			cancer_participant_tumour
		''')
		dod_participant = lab_to_df(
			sql_query=query1,
			dr=self.version
			)
		dod_participant['disease_type'] = translateicd(
			dod_participant['diagnosis_icd_code'],
			lookups=ic_lookup
			)
		dod_participant.sort_values(
			['diagnosis_date'],
			ascending=True,
			inplace=True)
		dod_participant.drop_duplicates(
			['participant_id', 'disease_type'],
			keep='first',
			inplace=True)

		query2=(f'''
		SELECT 
			DISTINCT participant_id, diagnosisdatebest, site_icd10_o2
		FROM
			av_tumour
		''')
		av_dod_participant = lab_to_df(
			sql_query=query2,
			dr=self.version
			)
		av_dod_participant['disease_type'] = translateicd(
			av_dod_participant['site_icd10_o2'],
			lookups=ic_lookup
			)
		av_dod_participant = (av_dod_participant
			.sort_values(
				['diagnosisdatebest'],
				ascending=True
				)
			.drop_duplicates(
				['participant_id', 'disease_type'],
				keep='first'
				)
			.rename(
				columns={'diagnosisdatebest':'diagnosis_date'}
				)
			)

		av_dod_full = (
			av_dod_participant[[
			'participant_id', 
			'diagnosis_date',
			'disease_type'
			]]
			.merge(dod_participant[[
				'disease_type',
				'participant_id',
				'diagnosis_date']],
				how='outer')
			.sort_values(['participant_id', 'diagnosis_date'], ascending=True)
			.drop_duplicates(['participant_id', 'disease_type'], keep='first')
			)


		imp_diag = pd.merge(
			fullca,
			av_dod_full, 
			how='left')
		imp_diag = imp_diag.merge(par, how='left')

		imp_diag['date_of_consent'] = imp_diag['date_of_consent'].apply(
			pd.to_datetime,
			format='%Y-%m-%d'
			)
		imp_diag['diagnosis_date'] = imp_diag['diagnosis_date'].apply(
			pd.to_datetime,
			format='%Y-%m-%d'
			)
		imp_diag['diff'] = (
			imp_diag['date_of_consent'] - imp_diag['diagnosis_date']
			).dt.days
		# impute DF is the key ref table. with average time between
		# diagnosis and enrollment for each disease type.
		impute_df = (imp_diag
			.loc[(imp_diag['diff'] >= 0), ['disease_type', 'diff']]
			.groupby(['disease_type'], as_index=False)
			.mean()
			.round(0)
			.rename(columns={'diff':'meandiff'})
			)

		imp_diag = pd.merge(
			imp_diag,
			impute_df,
			on='disease_type',
			how='left')
		imp_diag=(imp_diag
			.drop_duplicates(['participant_id','disease_type'], keep='first')
		)
		imp_diag['diagnosis_date'] = (imp_diag['diagnosis_date']
			.fillna(
				imp_diag['date_of_consent'] 
				- pd.to_timedelta(imp_diag['meandiff'], unit='d')
			)
		)

		filled_df = pd.merge(self.no_diag[['disease_type','participant_id']],
			imp_diag[['participant_id', 'disease_type', 'diagnosis_date']],
			how='left',
			on=['participant_id','disease_type'])

		self.fill_diag = filled_df

		if self.impute:
			self.full_diag = pd.concat(
				[self.pid_diag, self.fill_diag]
				).drop_duplicates(
					['participant_id','disease_type']
				)
			# this dropna could also occur outside the class.
			# self.pid_diag = self.pid_diag.dropna(
			# 	axis=0,
			# 	subset=['diagnosis_date'],
			# 	inplace=True)

	def surv_time(self):
		
		if None in (self.pid_diag, self.ons, self.hes):
			raise ValueError('[ERROR] in calculating survival time.'
			'Have quer_ons, quer_dod, quer_hes and merge_dod been run?')

		date_cutoff = pd.to_datetime( 
			max(c.ons['date_of_death']),
			format='%Y-%m-%d')

		surv_dat = pd.merge(
			self.pid_diag,
			self.ons, 
			how='left', 
			on='participant_id')
		surv_dat = pd.merge(
			surv_dat, 
			self.hes, 
			how='left', 
			on='participant_id')

		for x in ['lastseen', 'date_of_death', 'diagnosis_date']:
			surv_dat[x] = surv_dat[x].apply(
				pd.to_datetime,
					format='%Y-%m-%d'
				)
		# set the last date based on death date, lastseen or the maximum
		# date seen in the HES data.
		surv_dat.loc[
			~surv_dat['date_of_death'].isna(), 'last_date'
			] = surv_dat['date_of_death']
		surv_dat.loc[
			(surv_dat['date_of_death'].isna())
			& (~surv_dat['lastseen'].isna()), 'last_date'
			] = surv_dat['lastseen']
		surv_dat.loc[
			(surv_dat['date_of_death'].isna()) 
			& (surv_dat['lastseen'].isna()), 'last_date'
			] = date_cutoff

		surv_dat['last_date'] = surv_dat['last_date'].apply(
			pd.to_datetime,
				format='%Y-%m-%d'
			)
		# add live / dead flag.
		surv_dat['status'] = [
			0 if pd.isnull(x) else 1 for x in surv_dat['date_of_death']
			]
		surv_dat['survival'] = (
			surv_dat['last_date'] - surv_dat['diagnosis_date']
			)
		# filter out those with a last date before the diagnosis date.
		surv_dat = surv_dat.loc[
			~(surv_dat['survival'].dt.days <= 0)
			& ~(surv_dat['survival'].isna())]

		surv_dat = surv_dat.drop_duplicates([
			'participant_id', 
			'disease_type'
			])

		surv_data = surv_dat[[
			'participant_id',
			'disease_type',
			'diagnosis_date',
			'last_date',
			'survival',
			'status']]

		self.surv_dat = surv_data

def query_ctd(
	df, 
	version, 
	genes, 
	clinsig=['(likely)pathogenic','LoF','path_LoF']):
	"""Returns TRUE/FALSE for samples in dataframe by querying labkey,

	Args:
		df (pd.DataFrame): a pandas dataframe with at a column of participant_id
		version (str): version of labkey to query, DR 15 or later.
		genes (list): list of gene names in string format.
		clinsig (list): list of clinical significance to include. default:
		'(likely)pathogenic','LoF','path_LoF', excluding 'other'
	"""
	sqlstr = (f'''
	SELECT 
		participant_id,
		gene,
		relevance
	FROM
		lists.cancer_tier_and_domain_variants
	WHERE
		GENE IN {*genes,} AND relevance IN {*clinsig,}
	''')
	snvdb = lab_to_df(sqlstr, dr=version)

	for strat in genes:
		df[strat] = np.where(
			surv_dat['participant_id'].isin(snvdb['participant_id']),
			True,
			False
			)
	return df


def kmsurvival(data, strata, output,  plt_title, plotting=True, table=True):
	"""Calculate and plot Kaplan-Meier median survival time using the Lifelines
	Python package. 

	Args:
		data (pd.DataFrame): Dataframe with at 3 columns: survival(int), 
	    	status (1/0 if the event of interest occured), group 
			(some indicator of grouping)
		strata (list): Which groups to be compared, present in the data
		output (str): path to output folder.
		plt_title (str): Title of kaplain-meier plot.
		plotting (bool, optional): Should the function output a plot. 
			Defaults to True.
		table (bool, optional): should the function save a .csv table. 
			Defaults to True.

	Returns:
		pd.DataFrame: a Dataframe with each group a row, size of the group,
			number of events, median survival time and 95% confidence interval.
	"""
	
	from lifelines import KaplanMeierFitter
	from lifelines.utils import median_survival_times
	from matplotlib import pyplot as plt
	kmf = KaplanMeierFitter()
	# ngroup = range(0,len(pd.unique(data['group'])))
	out_d = []
	if plotting:
		ax = plt.subplot(111)
		plt.title(plt_title)
	for g in strata:
		s = (data['group'] == g)  
		fit = kmf.fit(
			data['survival'].dt.days[s],
			data['status'][s],
			label=g
			)
		out_d.append(
			{
				'group' : g,
				'n' : fit.event_observed.size,
				'events' : sum(surv_dat['status'][s]),
				'median' : fit.median_survival_time_,
				'upper_ci' : median_survival_times(
					fit.confidence_interval_
					).iloc[:,1].item(),
				'lower_ci' : median_survival_times(
					fit.confidence_interval_
					).iloc[:,0].item()
			}
		)
		if plotting:
			ax = kmf.plot_survival_function().plot(ax=ax)
	if plotting:
		plt.savefig(output+'surv.png',bbox_inches='tight', dpi=300)
		plt.close()
		plt.clf()
	outdf = pd.DataFrame(out_d)
	if table:
		outdf.to_csv(output+'surv.csv')
	# note lifelines uses the greenwood formulation for confidence intervals
	# these can be be larger than the bounds of survival and may differ
	# from the confidence intervals calculated by the R survival package.
	print((
		f'''strata\t| n {'':<6} | events {'':<6} | median {'':<6} |'''
		f''' 0.95UCL{'':<6}| 0.95LCL{'':<6}\n'''
		'''--------------------------------------------------------'''
		'''-----------------------\n'''),
		end=''
	)
	outdf.apply(
		lambda x:
			print(
				(f'''{x['group']}\t| {x['n']:<8} | {x['events']:<13} |'''
				f'''{x['median']:<14} | {x['upper_ci']:<13}| '''
				f'''{x['lower_ci']:<2}\n'''),
				end=''
			),
			axis=1
	)
	if not plotting:
		return kmf.plot_survival_function(), outdf
	return outdf


######################
# main
######################
if __name__ == '__main__':

	options = argparser()
	if options.out[-1] != '/':
		options.out=options.out+'/'

	# grab all cancer participants to find those with unknown diagnosis dates
	# or dates of death.
	if options.disease_type:
		options.disease_type = [x.upper() for x in options.disease_type]
		ca_query = (
			f'''
			SELECT 
				participant_id, disease_type
			FROM
				lists.cancer_analysis
			WHERE 
				tumour_type = 'PRIMARY'
				AND
				disease_type IN {*options.disease_type,}
			''')
	else:
		ca_query = (
			f'''
			SELECT 
				participant_id, disease_type
			FROM
				lists.cancer_analysis
			WHERE 
				tumour_type = 'PRIMARY'
			''')

	ca_analysis = lab_to_df(
		sql_query =ca_query,
		dr=options.version)

	c = Survdat(
		ca_analysis,
		ca_analysis['participant_id'], 
		options.version, 
		True)
	c.quer_ons()  # get survival data : c.ons
	c.quer_hes()  # query HES data for date of last follow up : c.hes
	c.quer_dod()  # get date of diagnosis :c.dod
	c.merge_dod()  # match date of diagnosis with cohort: c.pid_diag, c.no_diag
	c.dod_impute()  # impute date of diagnosis from average per disease type c.full_diag
	c.surv_time()  # use ons, hes, dod and pid_diag for survival data.


	############################################
	# bringing together the time to event data
	############################################
	date_cutoff = pd.to_datetime( 
		max(c.ons['date_of_death']),
		format='%Y-%m-%d')

	surv_dat = pd.merge(c.pid_diag, c.ons, how='left', on='participant_id')
	surv_dat = pd.merge(surv_dat, c.hes, how='left', on='participant_id')

	for x in ['lastseen', 'date_of_death', 'diagnosis_date']:
		surv_dat[x] = surv_dat[x].apply(
			pd.to_datetime,
				format='%Y-%m-%d'
			)
	# set the last date based on death date, lastseen or the maximum
	# date seen in the HES data.
	surv_dat.loc[
		~surv_dat['date_of_death'].isna(), 'last_date'
		] = surv_dat['date_of_death']
	surv_dat.loc[
		(surv_dat['date_of_death'].isna())
		& (~surv_dat['lastseen'].isna()), 'last_date'
		] = surv_dat['lastseen']
	surv_dat.loc[
		(surv_dat['date_of_death'].isna()) 
		& (surv_dat['lastseen'].isna()), 'last_date'
		] = date_cutoff

	surv_dat['last_date'] = surv_dat['last_date'].apply(
		pd.to_datetime,
			format='%Y-%m-%d'
		)
	# add live / dead flag.
	surv_dat['status'] = [
		0 if pd.isnull(x) else 1 for x in surv_dat['date_of_death']
		]
	surv_dat['survival'] = (surv_dat['last_date'] - surv_dat['diagnosis_date'])
	# filter out those with a last date before the diagnosis date.
	surv_dat = surv_dat.loc[
		~(surv_dat['survival'].dt.days <= 0)
		& ~(surv_dat['survival'].isna())]

	surv_dat = surv_dat.drop_duplicates([
		'participant_id', 
		'disease_type'
		])

	surv_dat = surv_dat[[
		'participant_id',
		'disease_type',
		'diagnosis_date',
		'last_date',
		'survival',
		'status']]

	
	####################
	# include SNV data
	####################
	# we create strata for multiple genes with assign_groups; full:
	# BRACA mut, PIK3CA mut
	# BRACA mut, PIK3CA wt
	# BRACA wt, PIK3CA mut
	# BRACA wt, PIK3CA wt
	surv_dat = query_ctd(  # loading in snvdb can be slow.
		df = surv_dat,
		version='main-programme/main-programme_v16_2022-10-13',
		genes=options.genes
		)

	# TODO, return a mapping for the and/ or too to later adjust labels with.
	if not options.strata == 'full':
		grouped_dat = assign_groups(
			dataframe=surv_dat, 
			vars=options.genes, 
			type=options.strata
			)
	else:
		grouped_dat, mapping = assign_groups(
			dataframe=surv_dat, 
			vars=options.genes, 
			type=options.strata
			)

	###
	# calculate KM survival
	###
	# TODO
	# check if And, Or, full all work.
	# How do we give names/labels to the groups? - 
	# 	we would have to keep track of their mapping.



	dat = kmsurvival(
		data=grouped_dat,
		strata=pd.unique(grouped_dat['group']),
		output=options.out,
		plt_title=options.plttitle,
		)
