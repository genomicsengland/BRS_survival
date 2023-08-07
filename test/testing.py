import pandas

cohort1 = {
	'terms':[
		"Early onset and familial Parkinson's Disease",
		"Complex Parkinsonism (includes pallido-pyramidal syndromes)"
		],
	'hpo':[
		'HP:0001300',
		'HP:0002548',
		],
	'icd10':[
		'Q780',
		'M80',
		'M81',
		'M82'
		],
	'cancer_terms':[
		'COLORECTAL',
		]
		}
	
version = "main-programme/main-programme_v17_2023-03-30"

cohort = Cohort(featdict=cohort1, version = version)


cohort.get_term_pids()
cohort.pids['dterm']
cohort.dterms
cohort.dterm_table

cohort.get_hpo_pids()
cohort.hpo_table
cohort.pids

cohort.get_cterm_pids()
cohort.cterm_table
cohort.pids['cterm']
cohort.pids

cohort.get_icd10_pids(limit=['all'])
cohort.icd10s
cohort.pids['icd10']
cohort.icd10_table

# get additional cohort info:
cohort.age()
cohort.age_table


cohort.ancestry()
cohort.ancestry_table

cohort.sex()
cohort.sex_table

cohort.quer_ons()
cohort.mortality_table


cohort.summary_stats()
cohort.summary

cohort.ontology_stats()
cohort.ont_summary
cohort.ont_vcount
cohort.ont_vcount['icd10_full']
cohort.ont_vcount['icd10_simple']
cohort.icd10_overlap_matrix


# te inflated counts are probably coming from the overlap. Where multiple ICD-10 codes are differentiationg patients multiple times.
# icd_coh1 = icd10
cohort1 = c1
icd_coh1 = icd10.copy()

icd_coh1.code.value_counts()
# this can be one of the experts. (as an excel sheet)
icd_coh1.groupby(['participant_id', 'code']).size()


icd_coh1.loc[icd_coh1['code']=='Q780']

icd_coh1.loc[(
    (icd_coh1['code']=='Q780')
    & (~icd_coh1['participant_id'].isin(icd_coh2['participant_id']))
    )]

cohort1 = c2

icd_coh2 = icd10.copy()

icd_coh2.code.value_counts()

icd_coh1.groupby(['code','participant_id']).size().unstack()
icd10['code'].value_counts()

### inlcude in cohort.py.

icd10.loc[icd10['participant_id']==111001038]

# how come I'm finding more of the Q780 participants in the multi search
# compared to the single search?
len(icd_coh1.loc[icd_coh1['code']=='Q780', 'participant_id'].unique())
len(icd_coh2['participant_id'].unique())

# welke patients zitten in coh1 die niet in coh2 zitten?

missing_pats = icd_coh1.loc[(
    (icd_coh1['code']=='Q780')
    & (~icd_coh1['participant_id'].isin(icd_coh2['participant_id']))
    )]

# lets look at the diagcodes for these participants if we load in c1

(icd10['diag_all']
    .reset_index(drop=True)
    .str.extractall(p)
    )

icd10.loc[icd10['participant_id'].isin(
    [115017283,111000073, 111000762,]),'diag_all'].str.extractall(p)].unstack()

icd10['diag_all'].str.extractall(p)

# it was the reset_index that was missing.