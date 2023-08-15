import pandas as pd
import gelpack
from gelpack import cohort

gel.Cohort()

from gelpack.gel_utils import lab_to_df
from gelpack.cohort import Cohort

## create a cohort 
cohort1 = {
	'terms':[
		"Osteogenesis imperfecta",
		],
	'icd10':[
		'Q780',
		],
		}
	
version = "main-programme/main-programme_v17_2023-03-30"

cohort = Cohort(
	featdict=cohort1, 
	version = version, 
	name='Bone')


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

cohort.concat_all()
cohort.all_age


############# testing cancer cohort creation ################
version = "main-programme/main-programme_v17_2023-03-30"

# input a dict of hpo-terms and/or icd-10 codes
# we can select samples based on ICD-10 codes.
# cancer_terms
# study abbreviations
# 
# and modify these with histology codes. (incl, excl)
# 
cohort2_dict = {
	'icd10':[
		'C349',
		],
	'cancer_terms':[
		'LUNG',
        ],
    'cancer_abbr':[
        'LUAD',
        'LUSC',
        'LUG'
        ],
# this can get a bit murky - as histology-codes are tied to particular disease types.
# if these are not included we can assume we want all the histology types for a study abbreviation.
    'abbr_histology_incl':[
        ('LUG','8046/3')
        ],
    'abbr_histology_excl':[]
}

cohort2 = Cohort(
	featdict=cohort2_dict, 
	version=version,
	name='NSCLC')


cohort2.get_icd10_pids(limit=['all'])
cohort2.icd10_table
cohort2.get_cterm_pids()
cohort2.cterm_table
cohort2.age()
cohort2.ancestry()
cohort2.sex()
cohort2.quer_ons()
cohort2.concat_all()

# visualise the single cohort.

import gel_vis
fig2 = vis_cohorts(cohort, show=False)
fig2.show()


# get all participant ids associated with 
# icd-10 -> we don't know if we have a sample and/or if this sample is for this specific icd-10 code
# cancer_term -> cancer_participant_disease: not that accurate, and we don't kow if we a participant has a sample.
# cancer_abbr -> only availabe in cancer_analysis. They are more accurate than cancer_terms + we have a sample.
#              -> we may want to specify a specfic histology to include / exclude for a cohort. 
#                   - are these only associated with ICD-10 and/ or study_abbreviations? 
#                   - can we do this based on cancer analysis or do we need to grab our secondary data resources?