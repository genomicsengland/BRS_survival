#!/usr/bin/env python

#### Christian J. Bouwens
#### BRS team
#### generate Kaplann-meier survival curves on domain variants and
#### disease types.
#### last update: 2024.12.04

#########################
import argparse
import numpy as np
from pathlib import Path
import pandas as pd
import warnings
import labkey
import re
from functools import reduce
from gelpack.gel_utils import lab_to_df, translateicd, assign_groups, create_name_map
import lifelines





# this version of pandas is still suffering from the _str_len bug.
# pd.set_option("display.max_columns", None)
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
		default=False,
		help='Flag: should missing diagnosis dates be imputated?'
			'Default=False')
	parser.add_argument('-v', "--version",
		dest='version',
		default="main-programme/main-programme_v17_2023-03-30",
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
	# for testing purposes:
	# argv = '--genes KRAS -s full -o ./test'.split()
	# options = parser.parse_args(argv)
	options=parser.parse_args()
	
	return options


class Survdat():
	
	def __init__(self, df, pids, version, impute):

		self.pids = {}
		self.platekeys = None
		self.ca = df
		self.pids['custom'] = pids
		self.version = version
		self.impute = impute
		self.feature_tables = {}
		self.sample_tables = {}
	 

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
			WHERE participant_id IN {*self.pids['custom'],}
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
			WHERE participant_id IN {*self.pids['custom'],}
		''')
		ons2 = lab_to_df(
			sql_query=query3,
			dr=self.version)

		ons1.rename(columns={'death_date':'date_of_death'}, inplace=True)
		ons = ons1.merge(ons2, how='outer')
		ons.sort_values(
			by='date_of_death',
			ascending=True, 
			inplace=True,
			na_position='last')
		ons.drop_duplicates(
			subset='participant_id', 
			keep='first', 
			inplace=True)

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
			('apc','epistart'),
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
		suffix = ('_apc', '_ae', '_op', '_cc')
		for i,df in enumerate(dfs):
			df.rename({'LastSeen':'LastSeen'+suffix[i]}, inplace=True, axis=1)

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
			merged_lastseen['participant_id'].isin(self.pids['custom'])
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
				participant_id IN {*self.pids['custom'],}
		''')
		dod_participant = lab_to_df(
			sql_query=query1,
			dr=self.version
			)
		dod_participant['disease_type'] = translateicd(
			dod_participant['diagnosis_icd_code']
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
				participant_id IN {*self.pids['custom'],}
		''')
		av_dod_participant = lab_to_df(
			sql_query=query2,
			dr=self.version
			)
		av_dod_participant['disease_type'] = translateicd(
			av_dod_participant['site_icd10_o2']
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

		# add in cancer registry data here:
		# Append this based on the data release.
		query3=(f'''
			SELECT
				DISTINCT participant_id, event_date, cancer_site
			FROM
				cancer_registry
			WHERE participant_id IN {*self.pids['custom'],}
			''')
		nhsd_dod = lab_to_df(
			sql_query=query3,
			dr=self.version
			)
		# some of the cancer sites are integers - leading to errors in
		# translateicd
		nhsd_dod['cancer_site'] = nhsd_dod['cancer_site'].map(str)
		nhsd_dod.rename(
			{'event_date':'diagnosis_date'}, 
			axis=1, 
			inplace=True
			)
		nhsd_dod['disease_type'] = translateicd(
			nhsd_dod['cancer_site']
			)

		av_dod_participant_merged = av_dod_participant[[
			'participant_id', 
			'diagnosis_date',
			'disease_type'
			]].merge(
				dod_participant[[
					'disease_type',
					'participant_id',
					'diagnosis_date']],
				how='outer'
				)

		av_nhsd_dod_participant_merged = (
			av_dod_participant_merged
			.merge(nhsd_dod[[
				'disease_type',
				'participant_id',
				'diagnosis_date']],
				how='outer'
				)
			.sort_values(['participant_id', 'diagnosis_date'], ascending=True)
			.drop_duplicates(['participant_id', 'disease_type'], keep='first')
			)

		self.dod = av_nhsd_dod_participant_merged


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
			dod_participant['diagnosis_icd_code']
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
			av_dod_participant['site_icd10_o2']
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
		# add in cancer registry data here:
		query3=(f'''
		SELECT
			DISTINCT participant_id, event_date, cancer_site
		FROM
			cancer_registry
		''')
		nhsd_dod = lab_to_df(
			sql_query=query3,
			dr=self.version
			)
		# some of the cancer sites are integers - leading to errors in
		# translateicd
		nhsd_dod['cancer_site']= nhsd_dod['cancer_site'].map(str)
		nhsd_dod.rename(
			{'event_date':'diagnosis_date'}, 
			axis=1, 
			inplace=True
			)
		nhsd_dod['disease_type'] = translateicd(
			nhsd_dod['cancer_site']
			)
		# TODO: increaes ICD10 code regex coverage.
		# nhsd_dod.loc[
		# 	nhsd_dod['disease_type']=='OTHER',['cancer_site']
		# 	].value_counts()

		av_dod_half = av_dod_participant[[
			'participant_id', 
			'diagnosis_date',
			'disease_type'
			]].merge(dod_participant[[
				'disease_type',
				'participant_id',
				'diagnosis_date']],
			how='outer'
			)
			
		av_dod_full = (
			av_dod_half
			.merge(nhsd_dod[[
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

		date_cutoff = pd.to_datetime( 
			max(self.ons['date_of_death']),
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
		# surv_dat = surv_dat.loc[
		# 	~(surv_dat['survival'].dt.days <= 0)
		# 	& ~(surv_dat['survival'].isna())]

		surv_dat = surv_dat.drop_duplicates([
			'participant_id', 
			'disease_type'
			])

		surv_data = surv_dat[[
			'participant_id',
			'disease_type',
			'diagnosis_date',
			'date_of_death',
			'last_date',
			'survival',
			'status']]

		self.surv_dat = surv_data


# force list                
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
	from . import force_list
	genes_list = force_list(genes)

	if len(genes_list) > 1:
		sqlstr = (f'''
		SELECT 
			participant_id,
			gene,
			relevance
		FROM
			lists.cancer_tier_and_domain_variants
		WHERE
			gene IN {*genes_list,} AND relevance IN {*clinsig,}
		''')
	else:
		sqlstr = (f'''
		SELECT 
			participant_id,
			gene,
			relevance
		FROM
			lists.cancer_tier_and_domain_variants
		WHERE
			gene IN ('{genes_list[0]}') AND relevance IN {*clinsig,}
		''')

	snvdb = lab_to_df(sqlstr, dr=version)

	for strat in genes_list:
		df[strat] = np.where(
			df['participant_id'].isin(
				snvdb.loc[snvdb['gene'] == strat, 'participant_id']
				),
			True,
			False
			)
	return df


def cohort_surv(cohorts, ca_df=None, group_pids=None, group_names=None):
	"""perform KM survival calculations on one or multiple Cohorts. If a list of
	Cohorts is given they should be named. If multiple groups are part of the same
	Cohort they should be assigned a group in Cohort.groups. Finaly, the function
	can also take in a nested lids of group_pids and associated group_names.

	Args:
		cohorts (Cohort class): one or multiple Cohorts presented as list.
		ca_df (pd.Dataframe, optional): A dataframe of GEL cancer data, in particular
			The dataframe should contain participant ids and disease_types.
			Defaults to None.
		group_pids (list, optional):  A nested list of participant_ids, each group 
			to be compared in its own list. Defaults to None.
		group_names (list, optional): a list of names for each group. Defaults to None.


	Returns:
		class.Survdat: custom class of diagnosis and survival data per participant.
		pd.DataFrame: A dataframe indiciating how many participants the survival
			workflow was unable to find accurate diagnosis dates for.
		pd.DataFrame: The survival data and groupings used in KM-survival analysis.
		pd.Dataframe: A mapping of group names.
	"""
	from itertools import chain
	from functools import reduce
	from gelpack.gel_utils import lab_to_df
	import numpy as np

	# from gelpack.survival import Survdat
	# from gelpack.gel_utils import assign_groups
	# check if the cohort is one or multiple Cohorts.
	
	try:
		iterator = iter(cohorts)
	except:
		mult_cohorts = False
		# larger than 1, otherwise there would be no comparison possible
		if len(cohorts.groups) > 1:
			group_pids = list(cohorts.groups.values())
			len(group_pids)
			group_names = list(cohorts.groups.keys())
		else:
			if len(cohorts.groups) == 1:
				warnings.warn('Only one group found in the Cohort.')
			if not group_pids:
				raise ValueError('Please provide lists of participant ids to compare.')
			if not group_names:
				raise ValueError('Please provide a corresponding group name for each list of participant ids.')
			if (len(group_pids) != len(group_names)):
				raise ValueError('not every group of participant ids has an associated name.')

		all_pids = list(chain(*group_pids))
		version = cohorts.version
	else:
		mult_cohorts=True
		# if we have multiple cohorts we assume each cohort is compared to the other cohorts.
		# in which case we don't need group_names, as we already have them.
		group_pids = [list(x.all_pids) for x in cohorts]
		group_names = [x.name for x in cohorts]
		if (len(group_pids) != len(group_names)):
				raise ValueError('not every cohort has an associated name.')

		all_pids = pd.concat([x.all_pids for x in cohorts])
		
		# concat_featdicts(featdict_llist):
		# 	from functools import reduce
		# 	for dict in featdict:
		# 		for df_merged = reduce(lambda left,right: 
		# 			pd.concat(
		# 				left, right,
		# 				axis=1), featdict)

		version = cohorts[1].version 


	# double check there are no overlapping participants in the groups:
	# TODO: this only checks if there is an overlap in the first and subsequent groups.
	overlap = reduce(np.intersect1d, (group_pids))
	# overlap = len(set(group_pids[0]).intersection(*group_pids[1:]))
	if len(overlap) > 0:
		print(overlap)
		raise ValueError('Overlap found in groups. ' +
		'Make sure each group is unique to avoid confounding.')

	if ca_df is None:
		ca_df = lab_to_df(
			sql_query=
				'''
				SELECT
					participant_id,
					disease_type
				FROM
					cancer_analysis
				''',
			dr=version)	

	simp_cohort = Survdat(
		version=version, 
		pids=all_pids,
		df=ca_df,
		impute=False
		)

	simp_cohort.quer_ons()  # get survival data : c.ons
	simp_cohort.quer_hes()  # query HES data for date of last follow up : c.hes
	simp_cohort.quer_dod()  # get date of diagnosis :c.dod
	simp_cohort.merge_dod() 
	simp_cohort.surv_time()  # use ons, hes, dod and pid_diag for survival data.
	
	simp_cohort.age()
	simp_cohort.ancestry()
	simp_cohort.sex()

	# simp_cohort.ca_sample_data()
	
	simp_cohort.concat_all()

	
	# merge the survival data with its features:
	all_surv = pd.merge(
		simp_cohort.surv_dat,
		simp_cohort.all_data,
		how='left',
		on='participant_id'
	)
	simp_cohort.group_features = pd.DataFrame({
		'participant_id':simp_cohort.pids['custom']
		})
	for x in group_names:
		simp_cohort.group_features[x] = pd.NA
	
	for feat, source in zip(group_names,group_pids):
		simp_cohort.group_features[feat] = np.where(
			simp_cohort.group_features['participant_id'].isin(source),
			True,
			False
		)
	groups, mapping = assign_groups(
		simp_cohort.group_features,
		vars=group_names,
		type='full')
	map_dict = create_name_map(mapping, group_names)
	

	survival_data = pd.merge(groups, all_surv, on = 'participant_id', how='left')
	# survival_data.dropna()
	missing_pgroup = pd.concat([
		survival_data.loc[survival_data['diagnosis_date'].isna(),'group'].value_counts(),
		survival_data['group'].value_counts()
		], axis=1)
	missing_pgroup.columns = ['missing','all']
	missing_pgroup['freq'] = missing_pgroup['missing'] / missing_pgroup['all'] 
	# sdat = survival_data.dropna()
	# TODO add a little boxplot that shows how much NA we have for each group

	return simp_cohort, missing_pgroup, survival_data, map_dict

def add_at_risk_counts(
    *fitters,
    labels= None,
    rows_to_show=None,
    ypos=-0.6,
    xticks=None,
    ax=None,
    at_risk_count_from_start_of_period=False,
    mask=False,
    **kwargs,
	):
    """
	--- Fork from Lifelines ---
    Add counts showing how many individuals were at risk, censored, and observed, at each time point in
    survival/hazard plots.

    Tip: you probably want to call ``plt.tight_layout()`` afterwards.

    Parameters
    ----------
    fitters:
      One or several fitters, for example KaplanMeierFitter, WeibullFitter,
      NelsonAalenFitter, etc...
    labels:
        provide labels for the fitters, default is to use the provided fitter label. Set to
        False for no labels.
    rows_to_show: list
        a sub-list of ['At risk', 'Censored', 'Events']. Default to show all.
    ypos:
        make more positive to move the table up.
    xticks: list
        specify the time periods (as a list) you want to evaluate the counts at.
    at_risk_count_from_start_of_period: bool, default False.
        By default, we use the at-risk count from the end of the period. This is what other packages, and KMunicate suggests, but
        the same issue keeps coming up with users. #1383, #1316 and discussion #1229. This makes the adjustment.
    ax:
        a matplotlib axes
    mask:
        mask counts lower than 5 with '<5'

    Returns
    --------
      ax:
        The axes which was used.

    Examples
    --------
    .. code:: python

        # First train some fitters and plot them
        fig = plt.figure()
        ax = plt.subplot(111)

        f1 = KaplanMeierFitter()
        f1.fit(data)
        f1.plot(ax=ax)

        f2 = KaplanMeierFitter()
        f2.fit(data)
        f2.plot(ax=ax)

        # These calls below are equivalent
        add_at_risk_counts(f1, f2)
        add_at_risk_counts(f1, f2, ax=ax, fig=fig)
        plt.tight_layout()

        # This overrides the labels
        add_at_risk_counts(f1, f2, labels=['fitter one', 'fitter two'])
        plt.tight_layout()

        # This hides the labels
        add_at_risk_counts(f1, f2, labels=False)
        plt.tight_layout()

        # Only show at-risk:
        add_at_risk_counts(f1, f2, rows_to_show=['At risk'])
        plt.tight_layout()

    References
    -----------
     Morris TP, Jarvis CI, Cragg W, et al. Proposals on Kaplan–Meier plots in medical research and a survey of stakeholder views: KMunicate. BMJ Open 2019;9:e030215. doi:10.1136/bmjopen-2019-030215

    """
    
    from matplotlib import pyplot as plt
    from lifelines.plotting import move_spines, remove_spines, remove_ticks, is_latex_enabled
    if ax is None:
        ax = plt.gca()
    fig = kwargs.pop("fig", None)
    if fig is None:
        fig = plt.gcf()
    if labels is None:
        labels = [f._label for f in fitters]
    elif labels is False:
        labels = [None] * len(fitters)
    if rows_to_show is None:
        rows_to_show = ["At risk", "Censored", "Events"]
    else:
        assert all(
            row in ["At risk", "Censored", "Events"] for row in rows_to_show
        ), 'must be one of ["At risk", "Censored", "Events"]'
    n_rows = len(rows_to_show)

    # Create another axes where we can put size ticks
    ax2 = plt.twiny(ax=ax)
    # Move the ticks below existing axes
    # Appropriate length scaled for 6 inches. Adjust for figure size.
    ax_height = (
        ax.get_position().y1 - ax.get_position().y0
    ) * fig.get_figheight()  # axis height
    ax2_ypos = ypos / ax_height

    move_spines(ax2, ["bottom"], [ax2_ypos])
    # Hide all fluff
    remove_spines(ax2, ["top", "right", "bottom", "left"])
    # Set ticks and labels on bottom
    ax2.xaxis.tick_bottom()
    # Set limit
    min_time, max_time = ax.get_xlim()
    ax2.set_xlim(min_time, max_time)
    # Set ticks to kwarg or visible ticks
    if xticks is None:
        xticks = [xtick for xtick in ax.get_xticks() if min_time <= xtick <= max_time]
    ax2.set_xticks(xticks)
    # Remove ticks, need to do this AFTER moving the ticks
    remove_ticks(ax2, x=True, y=True)
    
    ticklabels = []

    for tick in ax2.get_xticks():
        lbl = ""

        # Get counts at tick
        counts = []
        
        for f in fitters:
            # this is a messy:
            # a) to align with R (and intuition), we do a subtraction off the at_risk column
            # b) we group by the tick intervals
            # c) we want to start at 0, so we give it it's own interval
            if at_risk_count_from_start_of_period:
                event_table_slice = f.event_table.assign(at_risk=lambda x: x.at_risk)
            else:
                event_table_slice = f.event_table.assign(
                    at_risk=lambda x: x.at_risk - x.removed
                )
            if not event_table_slice.loc[:tick].empty:
                event_table_slice = (
                    event_table_slice.loc[:tick, ["at_risk", "censored", "observed"]]
                    .agg(
                        {
                            "at_risk": lambda x: x.tail(1).values,
                            "censored": "sum",
                            "observed": "sum",
                        }
                    )  # see #1385
                    .rename(
                        {
                            "at_risk": "At risk",
                            "censored": "Censored",
                            "observed": "Events",
                        }
                    )
                    .fillna(0)
                )
                
                counts.extend([int(c) for c in event_table_slice.loc[rows_to_show]])
            else:
                counts.extend([0 for _ in range(n_rows)])
        if n_rows > 1:
            if tick == ax2.get_xticks()[0]:
                max_length = len(str(max(counts)))
                for i, c in enumerate(counts):
                    if i % n_rows == 0:
                        if is_latex_enabled():
                            lbl += (
                                ("\n" if i > 0 else "")
                                + r"\textbf{%s}" % labels[int(i / n_rows)]
                                + "\n"
                            )
                        else:
                            lbl += (
                                ("\n" if i > 0 else "")
                                + r"%s" % labels[int(i / n_rows)]
                                + "\n"
                            )
                    l = rows_to_show[i % n_rows]
                    s = (
                        "{}".format(l.rjust(10, " "))
                        + (" " * (max_length - len(str(c)) + 3))
                        + "{:>{}}\n".format(
                            str(c) if not mask or c==0 else (str(c) if int(c) >= 5 else '<5'), max_length)
                    )
                    lbl += s.format(c)
            else:
                # Create tick label
                lbl += ""
                for i, c in enumerate(counts):
                    if i % n_rows == 0 and i > 0:
                        lbl += "\n\n"
                    s = "\n{}"
                    if not mask:
                        lbl += s.format(c)
                    else:
                        lbl +=s.format(c if c >= 5 or c==0 else '<5')
        else:
            # if only one row to show, show in "condensed" version
            if tick == ax2.get_xticks()[0]:
                max_length = len(str(max(counts)))

                lbl += rows_to_show[0] + "\n"

                for i, c in enumerate(counts):
                    s = (
                        "{}".format(labels[i].rjust(10, " "))
                        + (" " * (max_length - len(str(c)) + 3))
                        + "{:>{}}\n".format(
                            str(c) if not mask or c==0 else (str(c) if int(c) >= 5 else '<5'), max_length)
                    )
                    lbl += s

                    # if not mask:
                    #     lbl += s.format(c)
                    # else:
                    #     lbl +=s.format(c if c >= 5 else '<5')
            else:
                # Create tick label
                lbl += ""
                for i, c in enumerate(counts):
                    s = "\n{}"
                    if not mask:
                        lbl += s.format(c)
                    else:
                        lbl +=s.format(c if c >= 5 or c==0 else '<5')
        ticklabels.append(lbl)
    # Align labels to the right so numbers can be compared easily
    ax2.set_xticklabels(ticklabels, ha="right", **kwargs)

    return ax

def logrank(data, mapping, mult_test='multivariate', weightings=None):
	"""Calculate a logrank p-value with the lifelines package. distinguishes between
	univariate models (2 groups) and multivariate models (<2 groups.)

	Args:
		data (pd.DataFrame): pandas dataframe with at least 3 columns: 
		group, survival, status
		mapping (dicionary): a dictionary mapping the groups to names, 
		essential is that the key corresponds to the groups in data
		mult_test (str, optional): If there are more than 2 groups what type of 
		statistical test should be applied? options are: 'pair', and 'multivariate'
		corresponding to lifelines pairwise_logrank_test() and 
		multivariate_logrank_test(). Defaults to 'pair'.
		weightings (str, optional): should any weightings be applied?
		refer to lifelines.. Defaults to None.

	Returns:
		pd.DataFrame: a dataframe with a p-value an stat data.
	"""
	groups=list(mapping.keys())
	# check if its a single test.
	if len(groups) == 2:
		from lifelines.statistics import logrank_test
		results = logrank_test(
			durations_A=data.loc[data['group']==groups[0],'survival'].dt.days,
			durations_B=data.loc[data['group']==groups[1],'survival'].dt.days,
			event_observed_A=data.loc[data['group']==groups[0],'status'],
			event_observed_B=data.loc[data['group']==groups[1],'status'],
			weightings=weightings)
	elif len(groups) > 2:
		if mult_test =='pair':
			from lifelines.statistics import pairwise_logrank_test
			results = pairwise_logrank_test(
				event_durations=data['survival'].dt.days,
				event_observed=data['status'],
				groups=data['group'],
				weightings=weightings)
		elif mult_test =='multivariate':
			from lifelines.statistics import multivariate_logrank_test
			results=multivariate_logrank_test(
				event_durations=data['survival'].dt.days,
				event_observed=data['status'],
				groups=data['group'],
				weightings=weightings
			)
	else:
		raise ValueError('not enough groups to test.')
	
	return results


def km_survival(
	data,
	strata,
	output,
	plt_title,
	map_dict,
	col=None,
	plotting='show',
	table=True,
	interval='days',
	at_risk_counts=False,
	weightings=None,
	mask=True):
	"""Calculate and plot Kaplan-Meier median survival time using the Lifelines
	Python package. 

	Args:
		data (pd.DataFrame): Dataframe with at 3 columns: survival(int), 
			status (1/0 if the event of interest occured), group 
			(some indicator of grouping)
		strata (list): Which groups to be compared, present in the data
		output (str): path to save output as.
		plt_title (str): Title of kaplain-meier plot.
		map_dict (dictionary, optional): Dictionary to rename the groups, reflected
			in both table and plot output. Note: when included strata should 
			match the renamed groups.
		plotting (str, optional): Should the function create a plot: 
			['show', 'save'], or return the ax: 'None'. Defaults to show.
		table (bool, optional): should the function save a .csv table. 
			Defaults to True.
		interval (str, optional): what time interval should the plot be on?
			options =['days','months']. Defaults to 'days'.
		at_risk_counts (bool, optional): to add summary tables. Defaults to False.
		weightings (str, optional): Should the logrank p-value be calculated, or
			alternative tests be used? options are “wilcoxon” for Wilcoxon
			(also known as Breslow), “tarone-ware” for Tarone-Ware,
			“peto” for Peto test. Defaults to logrank test.

	Returns:
		pd.DataFrame: a Dataframe with each group a row, size of the group,
			number of events, median survival time and 95% confidence interval.
	"""
	from lifelines import KaplanMeierFitter
	from gelpack.survival import add_at_risk_counts
	# from lifelines.plotting import add_at_risk_counts point towards fork.
	from lifelines.utils import median_survival_times  # do we use this?
	from matplotlib import pyplot as plt
	
	kmf = KaplanMeierFitter()
	# ngroup = range(0,len(pd.unique(data['group'])))
	out_d = []
	if plotting is not None:
		# ax = plt.subplot(111, box_aspect=0.3) allow adjusting of aspect ratio
		fig, ax = plt.subplots()
		plt.title(plt_title)
	df = data.copy()
	# this renaming section forces users to rename the strata before the function is run.
	# perhaps approach differently?
	df.replace({'group':map_dict}, inplace=True)  

	stats = logrank(
		data=data,
		mapping=map_dict,
		mult_test='multivariate',
		weightings=weightings
		)
	# print(stats)
	aggfit = {}
	for g in strata:
		aggfit[g]=KaplanMeierFitter()
		s = (df['group'] == g)
		if interval =='days':
			aggfit[g].fit(
				df['survival'].dt.days[s],
				df['status'][s],
				label=g
				)
		elif interval == 'months':
			aggfit[g].fit(
				df['survival'][s]/np.timedelta64(1, 'M'),
				df['status'][s],
				label=g
				)
		out_d.append(
			{
				'group' : g,
				'n' : aggfit[g].event_observed.size,
				'events' : sum(df['status'][s]),
				'median' : aggfit[g].median_survival_time_,
				'upper_ci' : median_survival_times(
					aggfit[g].confidence_interval_
					).iloc[:,1].item(),
				'lower_ci' : median_survival_times(
					aggfit[g].confidence_interval_
					).iloc[:,0].item()
			}
		)
		# aggfit[g]=fit
		if plotting is not None: 
			if col is not None: 
				aggfit[g].plot_survival_function(color=col[g]).plot(ax=ax)
			else:
				aggfit[g].plot_survival_function().plot(ax=ax)
	# add the stats to the plot:
	if plotting is not None:
		if at_risk_counts:
			add_at_risk_counts(
				*list(aggfit.values()),
				ax=ax, 
				mask=mask,
				rows_to_show=['At risk', 'Censored'])
		if weightings == None:
			test='Logrank'
		else: 
			test = weightings
		
		plt.text(
			x=0.8, 
			y=0.8, 
			# transform=ax.transAxes,
			s=(f'''{test} p-value: {stats.p_value:.4}'''),
			fontsize=10,
			ha='right', va='top', transform=ax.transAxes)
		ax.set_xlabel(f'time ({interval})')

	if plotting == 'show':
		plt.tight_layout()
		plt.show()
	elif plotting == 'save':
		if output is None:
			raise InterruptedError('No output directory given.')
		# TODO allow differently named save (or use title)
		plt.savefig(output+'_kmplot.png', bbox_inches='tight', dpi=300)
		plt.close()
		plt.clf()
	else:
		return stats,aggfit
	
	outdf = pd.DataFrame(out_d)
	if table:
		if output is None:
			raise InterruptedError('No output directory given.')
		outdf.to_csv(output+'_surv_table.csv')
	# note lifelines uses the greenwood formulation for confidence intervals
	# these can be be larger than the bounds of survival and may differ
	# from the confidence intervals calculated by the R survival package.
	print((
		f'''strata {'':<8} | n {'':<6} | events {'':<6} | median {'':<6} |'''
		f''' 0.95UCL{'':<6}| 0.95LCL{'':<6}\n'''
		'''--------------------------------------------------------'''
		'''-----------------------\n'''),
		end=''
	)
	outdf.apply(
		lambda x:
			print(
				(f'''{x['group']:<10} \t| {x['n']:<8} | {x['events']:<13} |'''
				f'''{x['median']:<14} | {x['upper_ci']:<13}| '''
				f'''{x['lower_ci']:<2}\n'''),
				end=''
			),
			axis=1
	)
	# if not plotting:
	# 	return kmf.plot_survival_function(), outdf
	return outdf


def stepwise_coxph_regression(survival_df, time_col, event_col, to_test_cov):
	# check if covariates in survival_df are of type in, str or bool.
	trans=[]
	for cov in to_test_cov:
		if (survival_df[cov].dtype.kind) not in (['f','i', 'b']):
			trans.append(cov)
	if len(trans) > 0:
		survival_df = pd.get_dummies(
			data=survival_df,
			columns=trans,
			drop_first=False	
			)
	
	# make sure the survival column is int
	if survival_df[time_col].dtype.kind != 'i':
		try: 
			survival_df[time_col] = survival_df[time_col].dt.days
		except:
			raise AttributeError('Time column is neither an integer nor datetime.')
	import math
	# run univariate analysis for each covariate:
	incl_features = []
	p_values = []
	hazard_ratios = []
	hr_upper = []
	hr_lower = []

	# the one hot encoding above leads to us being unable to just use the
	# covariates argument  ..
	all_covariates = survival_df.drop([event_col,time_col],axis=1).columns.tolist()
	# perform the regression analysis with the first covariate.
	# then add covaria
	use_cov = []
	tmp_cov = all_covariates[0]
	# run the model with the first covariate, 
	# if its p<0.05: add it to use_cov and run the next use_cov+temp_cov.
	# how do we measure model improvement?


	# if the model improves again, add it to use_cov.
	for cov in all_covariates[1:]:
		cph.fit(
			survival_df[[event_col,time_col,cov]],
			duration_col=time_col,
			event_col=event_col)
		p_value = cph._compute_p_values()[-1]
		p_values.append(p_value)
		# c_index = cph.score(
		# 	survival_df[[event_col,time_col,cov]],
		# 	scoring_method="concordance_index"
		# 	)
		hazard_ratio = cph.hazard_ratios_.tolist()[-1]
		hazard_ratios.append(hazard_ratio)
		hazard_ratios_ci = cph.confidence_intervals_
		print(f'testing covariante: {cov}')
		print(f'p_value = {p_value}')
		print(f'hazard_ratios {hazard_ratio}')
		
		if p_value < 0.05:
			incl_features.append(cov)
		# if p_value[-1] < 0.1:
		try:
			lower = math.exp(hazard_ratios_ci.iloc[-1, 0])
		except OverflowError:
			lower = float('inf')
		hr_lower.append(lower)
		try:
			upper = math.exp(hazard_ratios_ci.iloc[-1, 1])
		except OverflowError:
			upper = float('inf')
		hr_upper.append(upper)

	
	
	univariate_df = pd.DataFrame(
		{
			'covariate':all_covariates,
			'p_value':p_values,
			'hazard_ratio':hazard_ratios,
			'hazard_ratio_upper_95':hr_upper,
			'hazard_ratio_lower_95':hr_lower
		}
		)

def coxph_regression(survival_df, time_col, event_col, covariates, formula=None):
	"""Performs Cox Proportional hazards regression. The script attempts to
	clean the data by forcing it to integer, one hot encoding strings and 
	removing highly correlative covariates

	Args:
		survival_df (pd.DataFrame): A dataframe of survival data, includes the
		time_col, event_col and covariates.
		time_col (str): Name of the column corresponding to the time till event 
		in the survival_df. Can be int or Datetime.
		event_col (str): Name of the column corresponding to the event-indication
		can only  an integer (0,1).
		covariates (list): list of strings with names corresponding to columns in
		survival_df

	Returns:
		univariate_df: a pd.DataFrame with outputs of the univariate coxph analysis.
		CoxPHFitter: a lifelines CoxPHFit associated with the data.
	"""
	
	import numpy as np
	from lifelines import CoxPHFitter
	# survival_df=survival_data.copy()
	# time_col='survival'
	# event_col='status'
	# covariates=['sex','group','age_at_consent','predicted_ancestry']
	survival_df=survival_df[[time_col, event_col] + covariates]
	
	# check if covariates in survival_df are of type in, str or bool.
	trans=[]
	for cov in covariates:
		if (survival_df[cov].dtype.kind) not in (['f','i', 'b']):
			trans.append(cov)
	if len(trans) > 0:
		survival_df = pd.get_dummies(
			data=survival_df,
			columns=trans,
			drop_first=False	
			)

	# make sure the survival column is int
	if survival_df[time_col].dtype.kind != 'i':
		try: 
			survival_df[time_col] = survival_df[time_col].dt.days
		except:
			raise AttributeError('Time column is neither an integer nor datetime.')


	# import matplotlib.pyplot as plt
	# import seaborn as sns	
	# plt.figure()
	# # TODO output plotslike KM.
	# dataplot = sns.heatmap(survival_df.corr(), cmap="YlGnBu", annot=True)
	# plt.show()

	# check correlation, remove variables that are highly corr.
	corr_matrix = survival_df.corr().abs()
	# Select upper triangle of correlation matrix
	upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
	# Find features with correlation greater than 0.95
	to_drop = [column for column in upper.columns if any(upper[column] > 0.95)]
	print(f'dropping correlative features: {*to_drop,}')
	# Drop features 
	survival_df.drop(to_drop, axis=1, inplace=True)

	cph = CoxPHFitter()
	import math
	# run univariate analysis for each covariate:
	incl_features = []
	p_values = []
	hazard_ratios = []
	hr_upper = []
	hr_lower = []

	# the one hot encoding above leads to us being unable to just use the
	# covariates argument..
	all_covariates = survival_df.drop([event_col,time_col],axis=1).columns.tolist()
	for cov in all_covariates:
		cph.fit(
			survival_df[[event_col,time_col,cov]],
			duration_col=time_col,
			event_col=event_col)
		p_value = cph._compute_p_values()[-1]
		p_values.append(p_value)
		# c_index = cph.score(
		# 	survival_df[[event_col,time_col,cov]],
		# 	scoring_method="concordance_index"
		# 	)
		hazard_ratio = cph.hazard_ratios_.tolist()[-1]
		hazard_ratios.append(hazard_ratio)
		hazard_ratios_ci = cph.confidence_intervals_
		print(f'testing covariante: {cov}')
		print(f'p_value = {p_value}')
		print(f'hazard_ratios {hazard_ratio}')
		
		if p_value < 0.05:
			incl_features.append(cov)
		# if p_value[-1] < 0.1:
		try:
			lower = math.exp(hazard_ratios_ci.iloc[-1, 0])
		except OverflowError:
			lower = float('inf')
		hr_lower.append(lower)
		try:
			upper = math.exp(hazard_ratios_ci.iloc[-1, 1])
		except OverflowError:
			upper = float('inf')
		hr_upper.append(upper)

	
	
	univariate_df = pd.DataFrame(
		{
			'covariate':all_covariates,
			'p_value':p_values,
			'hazard_ratio':hazard_ratios,
			'hazard_ratio_upper_95':hr_upper,
			'hazard_ratio_lower_95':hr_lower
		}
		)
		
	# now perform multivariate analysis with covariates that are interesting.
	mult_survival_df= survival_df[[time_col,event_col]+incl_features]

	try:
		cph.fit(mult_survival_df, duration_col=time_col, event_col=event_col)
	except:
		print('###################################################')
		print('Convergence error, trying again with penalizer=0.01')
		print('###################################################')
		cph = CoxPHFitter(penalizer=0.01)
		cph.fit(mult_survival_df, duration_col=time_col, event_col=event_col)

	# print summary
	cph.print_summary()  # access the individual results using cph.summary

	return univariate_df, cph


def add_at_risk_counts(
    *fitters,
    labels= None,
    rows_to_show=None,
    ypos=-0.6,
    xticks=None,
    ax=None,
    at_risk_count_from_start_of_period=False,
    mask=False,
    **kwargs,
):
    """
    Add counts showing how many individuals were at risk, censored, and observed, at each time point in
    survival/hazard plots.

    Tip: you probably want to call ``plt.tight_layout()`` afterwards.

    Parameters
    ----------
    fitters:
      One or several fitters, for example KaplanMeierFitter, WeibullFitter,
      NelsonAalenFitter, etc...
    labels:
        provide labels for the fitters, default is to use the provided fitter label. Set to
        False for no labels.
    rows_to_show: list
        a sub-list of ['At risk', 'Censored', 'Events']. Default to show all.
    ypos:
        make more positive to move the table up.
    xticks: list
        specify the time periods (as a list) you want to evaluate the counts at.
    at_risk_count_from_start_of_period: bool, default False.
        By default, we use the at-risk count from the end of the period. This is what other packages, and KMunicate suggests, but
        the same issue keeps coming up with users. #1383, #1316 and discussion #1229. This makes the adjustment.
    ax:
        a matplotlib axes
    mask:
        mask counts lower than 5 with '<5'

    Returns
    --------
      ax:
        The axes which was used.

    Examples
    --------
    .. code:: python

        # First train some fitters and plot them
        fig = plt.figure()
        ax = plt.subplot(111)

        f1 = KaplanMeierFitter()
        f1.fit(data)
        f1.plot(ax=ax)

        f2 = KaplanMeierFitter()
        f2.fit(data)
        f2.plot(ax=ax)

        # These calls below are equivalent
        add_at_risk_counts(f1, f2)
        add_at_risk_counts(f1, f2, ax=ax, fig=fig)
        plt.tight_layout()

        # This overrides the labels
        add_at_risk_counts(f1, f2, labels=['fitter one', 'fitter two'])
        plt.tight_layout()

        # This hides the labels
        add_at_risk_counts(f1, f2, labels=False)
        plt.tight_layout()

        # Only show at-risk:
        add_at_risk_counts(f1, f2, rows_to_show=['At risk'])
        plt.tight_layout()

    References
    -----------
     Morris TP, Jarvis CI, Cragg W, et al. Proposals on Kaplan–Meier plots in medical research and a survey of stakeholder views: KMunicate. BMJ Open 2019;9:e030215. doi:10.1136/bmjopen-2019-030215

    """
    
    from matplotlib import pyplot as plt
    from lifelines.plotting import move_spines, remove_spines, remove_ticks, is_latex_enabled
    if ax is None:
        ax = plt.gca()
    fig = kwargs.pop("fig", None)
    if fig is None:
        fig = plt.gcf()
    if labels is None:
        labels = [f._label for f in fitters]
    elif labels is False:
        labels = [None] * len(fitters)
    if rows_to_show is None:
        rows_to_show = ["At risk", "Censored", "Events"]
    else:
        assert all(
            row in ["At risk", "Censored", "Events"] for row in rows_to_show
        ), 'must be one of ["At risk", "Censored", "Events"]'
    n_rows = len(rows_to_show)

    # Create another axes where we can put size ticks
    ax2 = plt.twiny(ax=ax)
    # Move the ticks below existing axes
    # Appropriate length scaled for 6 inches. Adjust for figure size.
    ax_height = (
        ax.get_position().y1 - ax.get_position().y0
    ) * fig.get_figheight()  # axis height
    ax2_ypos = ypos / ax_height

    move_spines(ax2, ["bottom"], [ax2_ypos])
    # Hide all fluff
    remove_spines(ax2, ["top", "right", "bottom", "left"])
    # Set ticks and labels on bottom
    ax2.xaxis.tick_bottom()
    # Set limit
    min_time, max_time = ax.get_xlim()
    ax2.set_xlim(min_time, max_time)
    # Set ticks to kwarg or visible ticks
    if xticks is None:
        xticks = [xtick for xtick in ax.get_xticks() if min_time <= xtick <= max_time]
    ax2.set_xticks(xticks)
    # Remove ticks, need to do this AFTER moving the ticks
    remove_ticks(ax2, x=True, y=True)
    
    ticklabels = []

    for tick in ax2.get_xticks():
        lbl = ""

        # Get counts at tick
        counts = []
        
        for f in fitters:
            # this is a messy:
            # a) to align with R (and intuition), we do a subtraction off the at_risk column
            # b) we group by the tick intervals
            # c) we want to start at 0, so we give it it's own interval
            if at_risk_count_from_start_of_period:
                event_table_slice = f.event_table.assign(at_risk=lambda x: x.at_risk)
            else:
                event_table_slice = f.event_table.assign(
                    at_risk=lambda x: x.at_risk - x.removed
                )
            if not event_table_slice.loc[:tick].empty:
                event_table_slice = (
                    event_table_slice.loc[:tick, ["at_risk", "censored", "observed"]]
                    .agg(
                        {
                            "at_risk": lambda x: x.tail(1).values,
                            "censored": "sum",
                            "observed": "sum",
                        }
                    )  # see #1385
                    .rename(
                        {
                            "at_risk": "At risk",
                            "censored": "Censored",
                            "observed": "Events",
                        }
                    )
                    .fillna(0)
                )
                
                counts.extend([int(c) for c in event_table_slice.loc[rows_to_show]])
            else:
                counts.extend([0 for _ in range(n_rows)])
        if n_rows > 1:
            if tick == ax2.get_xticks()[0]:
                max_length = len(str(max(counts)))
                for i, c in enumerate(counts):
                    if i % n_rows == 0:
                        if is_latex_enabled():
                            lbl += (
                                ("\n" if i > 0 else "")
                                + r"\textbf{%s}" % labels[int(i / n_rows)]
                                + "\n"
                            )
                        else:
                            lbl += (
                                ("\n" if i > 0 else "")
                                + r"%s" % labels[int(i / n_rows)]
                                + "\n"
                            )
                    l = rows_to_show[i % n_rows]
                    s = (
                        "{}".format(l.rjust(10, " "))
                        + (" " * (max_length - len(str(c)) + 3))
                        + "{:>{}}\n".format(
                            str(c) if not mask or c==0 else (str(c) if int(c) >= 5 else '<5'), max_length)
                    )
                    lbl += s.format(c)
            else:
                # Create tick label
                lbl += ""
                for i, c in enumerate(counts):
                    if i % n_rows == 0 and i > 0:
                        lbl += "\n\n"
                    s = "\n{}"
                    if not mask:
                        lbl += s.format(c)
                    else:
                        lbl +=s.format(c if c >= 5 or c==0 else '<5')
        else:
            # if only one row to show, show in "condensed" version
            if tick == ax2.get_xticks()[0]:
                max_length = len(str(max(counts)))

                lbl += rows_to_show[0] + "\n"

                for i, c in enumerate(counts):
                    s = (
                        "{}".format(labels[i].rjust(10, " "))
                        + (" " * (max_length - len(str(c)) + 3))
                        + "{:>{}}\n".format(
                            str(c) if not mask or c==0 else (str(c) if int(c) >= 5 else '<5'), max_length)
                    )
                    lbl += s

                    # if not mask:
                    #     lbl += s.format(c)
                    # else:
                    #     lbl +=s.format(c if c >= 5 else '<5')
            else:
                # Create tick label
                lbl += ""
                for i, c in enumerate(counts):
                    s = "\n{}"
                    if not mask:
                        lbl += s.format(c)
                    else:
                        lbl +=s.format(c if c >= 5 or c==0 else '<5')
        ticklabels.append(lbl)
    # Align labels to the right so numbers can be compared easily
    ax2.set_xticklabels(ticklabels, ha="right", **kwargs)

    return ax

######################
# main
######################
if __name__ == '__main__':

	options = argparser()

	if options.out[-1] != '/':
		options.out=options.out+'/'
	Path(options.out).mkdir(parents=True, exists_ok=True)
	
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
		options.imputate_flag)
	c.quer_ons()  # get survival data : c.ons
	c.quer_hes()  # query HES data for date of last follow up : c.hes
	c.quer_dod()  # get date of diagnosis :c.dod
	c.merge_dod()  # match date of diagnosis with cohort: c.pid_diag, c.no_diag
	# dod_impute only does something if options.imputate_flag is true (--impute)
	c.dod_impute()  # impute date of diagnosis from average per disease type c.full_diag
	c.surv_time()  # use ons, hes, dod and pid_diag for survival data.
	# c.surv_dat will be the ultimate survival table.
	
	####################
	# include SNV data
	####################
	# we create strata for multiple genes with assign_groups; full:
	# BRACA mut, PIK3CA mut
	# BRACA mut, PIK3CA wt
	# BRACA wt, PIK3CA mut
	# BRACA wt, PIK3CA wt
	surv_dat = query_ctd(  # loading in snvdb can be slow.
		df = c.surv_dat,
		version=options.version,
		genes=options.genes
		)


	if not options.strata == 'full':
		grouped_dat = assign_groups(
			dataframe=surv_dat, 
			vars=options.genes, 
			type=options.strata
			)
		map_dict = {'group_1':'Mut','group_2':'WT'}
	else:
		grouped_dat, mapping = assign_groups(
			dataframe=surv_dat, 
			vars=options.genes, 
			type=options.strata
			)
		map_dict = create_name_map(
			mapping=mapping,
			genes=options.genes
			)

	###
	# calculate KM survival
	###
	dat = km_survival(
		data=grouped_dat,
		strata=map_dict.values(),
		map_dict=map_dict,
		output=options.out,
		plt_title=options.plttitle,
		plotting='save'
		)



