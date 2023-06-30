
#### Christian J. Bouwens
#### Miruna Carmen Barbu
#### BRSC team
#### generate a GEL cohort from HPO terms and ICD-10 codes
#### can be applied for cancer as well as rare disease cases.
#### last update: 2023.06.30

# 1. adopt the code for querying a Rare disease cohort
# 2. put it in cohort - class/object that can be passed to other functions.


import pandas as pd
import numpy as np
import labkey
from gelsurvival import gel_utils  # import the lab_to_df function.

# data release version
version = "main-programme/main-programme_v17_2023-03-30"

# input a dict of hpo-terms and/or icd-10 codes
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
        'G20',
        'G21',
        'G22'
        ]

}
    
icd10_lookup =  pd.read_csv(
    '/Users/christianbouwens/Documents'
    '/Internal/survival/surv_git/src/gelpack/coding19.tsv',
    sep='\t',
    usecols=['coding', 'meaning']
)

# mask function:
def maskfunc(x):
    np.where(1 <= x <= 4, -999, x )


###### load in the required tables #####
## participant id and normalised specific disease
rd_sql = (f'''
    SELECT
        DISTINCT participant_id,
        normalised_specific_disease
    FROM
        rare_diseases_participant_disease
    WHERE normalised_specific_disease IN {*cohort1['terms'],}
    ''')
rare_disease = lab_to_df(
    sql_query=rd_sql,
    dr=version
    )
    
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
        normalised_hpo_id IN {*cohort1['hpo'],} 
    ''')
hpo_terms = lab_to_df(
    sql_squery=hpo_sql,
    dr=version
    )

## participant_id and cancer_participant_disease
cancer_sql = (f'''
    SELECT 
        DISTINCT participant_id,
        cancer_disease_type
    FROM 
        cancer_participant_disease
    WHERE 
        cancer_disease_type IN {cohort1['cancer_terms']}''')
cancer_disease = lab_to_df(
    sql_query=cancer_sql,
    dr=version
)

  


## icd-10 diagnosis
hes_tables=[
	('apc'),
	('ae'),
	('cc')
	]

query = (f'''
	SELECT 
    	DISTINCT x.participant_id, 
        diag_all 
	FROM 
		hes_x as x 
	WHERE
        diag_all LIKE {
            " OR diag_all LIKE ".join(f"'%{i}%'" for i in cohort1['icd10'])}''')  # nested f string awh yeah.



q_apc, q_ae, q_cc = [
	query.replace('x',i) for i in hes_tables
	]


apc_icd10 = lab_to_df(
	sql_query=q_apc,
	dr=version
	)

ae_icd10= lab_to_df(
	sql_query=q_ae,
	dr=version
	)


cc_icd10= lab_to_df(
	sql_query=q_cc,
	dr=version
	)

dfs = [apc_icd10, ae_icd10, cc_icd10]
# (do we need to merge here? or just keep the independent tables?)

# extract codes from mortality


# mental health


