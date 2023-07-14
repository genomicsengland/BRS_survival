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
        'G22',
        'J690',
        'C50'
        ],
    'cancer_terms':[
        'LUNG',
        'BREAST'
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
# SQL will interpret double apostophe as single
term_clean = [i.replace("'", "''") for i in cohort1['terms']]
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
        normalised_hpo_id 
            IN ({', '.join(f"'{i}'" for i in cohort1['hpo'])}) 
    ''')
hpo_terms = lab_to_df(
    sql_query=hpo_sql,
    dr=version
    )

## participant_id and cancer_participant_disease
# have to put in a join clause as we cannot have a trailing comma for 
# a single disease_type.
format_ca_types =  ', '.join(f"'{i}'" for i in cohort1['cancer_terms'])
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
    dr=version
    )


## icd-10 diagnosis

# gathering all the sources of icd-10 codes, concatening and cleaning them up
# in the end we only want the relevant icd-10 codes (and one per participant).
# keeping in mind this top down approach is not perse relevant for cancer.
# which prefers a bottom-up approach- where we start from the cancer samples
#   we have, as opposed to starting from diagnosis codes.
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
            " OR diag_all LIKE ".join(f"'%{i}%'" for i in cohort1['icd10'])}''')  # nested f string awh yeah.

q_apc, q_ae, q_op = [
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
op_icd10= lab_to_df(
	sql_query=q_op,
	dr=version
	)


# extract codes from mortality
# only way to rename the long column to diag_all is by preceding the WHERE
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
            " OR diag_all LIKE ".join(f"'%{i}%'" for i in cohort1['icd10'])}''')
mortality = lab_to_df(
    sql_query=mortsql,
    dr=version
    )

dfs = [apc_icd10, ae_icd10, op_icd10, mortality]
icd10 = pd.concat(dfs)

# seperate the diag fields and keep unique+matching icd-10 codes per participant
# have to perform some regex gymnastics as extractall can't expand=True
p=r""+'|'.join(f"({i})"for i in cohort1['icd10'])+r""
icd10[[*cohort1['icd10']]] =(icd10['diag_all']
    .str.extractall(p)
    .groupby(level=0)
    .first()
    )
icd10 = (pd.melt(
    icd10,
    id_vars=['participant_id'],
    value_vars=[*cohort1['icd10']]
    )
    .dropna(axis=0)
    .drop_duplicates(['participant_id','value'])
    .drop(['variable'], axis=1)
	.rename({'value':'code'}, axis=1)
    )

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
                .join(f"'%{i}%'" for i in cohort1['icd10'])
            }
        OR REGEXP_REPLACE(diag_mhmdsec, 'b.', '') LIKE {
             " OR REGEXP_REPLACE(diag_mhmdsec, 'b.', '') LIKE "
             .join(f"'%{i}%'" for i in cohort1['icd10'])
        }
        '''.replace('b', '\\')

mhmd_sql, mhmdds_sql = [
    mh_sql.replace('x',i) for i in ['mhmd_v4_event', 'mhldds_event']
    ]

mhml = lab_to_df(
    sql_query=mhmd_sql,
    dr=version
    )

mhldds = lab_to_df(
    sql_query=mhmdds_sql,
    dr=version
    )

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
             .join(f"'%{i}%'" for i in cohort1['icd10'])
        })'''.replace(';', '\\')


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
    	dr=version
    	) for i in mhsds_sql_repl]

mhsds_diag = pd.concat(mhsds_diag_list).reset_index(drop=True)

mhsds_diag['code'] = (mhsds_diag['diag']
	.str.replace('.','', regex=False)
	.str.extract(r'([A-Z][0-9]+)')
)
mhsds_diag = (mhsds_diag
    .drop_duplicates(['participant_id','code'])
    .drop(['diag'],axis=1)
    )


#### merge mental health and ICD10 terms
icd_10 = pd.concat([icd10, mhsds_diag]).drop_duplicates(['participant_id','code'])



####### icd-10 cancer-specific
# we should display some kind of warning this is just to gather all known
# participants which may have had or get cancer. But not a proper approach
# to identify patients for which we have tumour samples.
# cancer invest sample pathology

sqlstr =  f'''
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
             .join(f"'%{i}%'" for i in cohort1['icd10'])
        })'''.replace(';', '\\')


## the R script is filtering out participants who may have multiple cancers here.
path_sql = f'''
    WITH table AS
    (
        SELECT
            participant_id,
            x as code
            
        FROM
            z
    )
    SELECT 
        DISTINCT participant_id,
            code
    FROM
        table
    WHERE 
        (REGEXP_REPLACE(code, '\\.', '') LIKE {
            " OR REGEXP_REPLACE(code, ';.', '') LIKE "
             .join(f"'%{i}%'" for i in cohort1['icd10'])
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
    )
    ]


cancer_sql_repl  = [path_sql
    .replace('x', i[1])
    .replace('z', i[0]) for i in c_table_str]


cancer_pathology, part_tumour, c_registry, rts, sact = [
	lab_to_df(
    	sql_query=i,
    	dr=version
    	) for i in cancer_sql_repl]

cancer_code_list= [
	lab_to_df(
    	sql_query=i,
    	dr=version
    	) for i in cancer_sql_repl]



cancer_diag = pd.concat(cancer_code_list).reset_index(drop=True)

cancer_diag['code'] = (cancer_diag['diag']
	.str.replace('.','', regex=False)
	.str.extract(r'([A-Z][0-9]+)')
)
cancer_diag = (cancer_diag
    .drop_duplicates(['participant_id','code'])
    .drop(['diag'],axis=1)
    )

# probably put in a different class?
# at this point you are no longer "gathering" participants
# but filling in details.
# Age

# Sex

# ancestry

# if we set the table strings in a dictionary per version - which is loaded 
# upon setting the version -> easy to do version control.