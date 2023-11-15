import pandas as pd
import numpy as np

from gelpack.gel_utils import lab_to_df
from gelpack.cohort import Cohort
version = "main-programme/main-programme_v17_2023-03-30"
## create a cohort 
cohort1 = {
	'terms':[
		"Osteogenesis imperfecta",
		],
	'icd10':[
		'Q780',
		'M81',
		'M80',
		'M82'
		],
		}


cohort1 = {
	'icd10':[
		'M81',
		'M80',
		'M82'
		],
		}

cohort = Cohort(
	featdict=cohort1, 
	version = version, 
	name='osteoporosis')

# cohort.get_term_pids()
# cohort.pids['dterm']
# cohort.dterms
# cohort.dterm_table

# cohort.get_hpo_pids()
# cohort.hpo_table
# cohort.pids

# cohort.get_cterm_pids()
# cohort.cterm_table
# cohort.pids['cterm']
# cohort.pids

cohort.get_icd10_pids(limit=['all'])
cohort.icd10s
cohort.pids['icd10']
cohort.icd10_table
cohort.pids
# get additional cohort info:
cohort.age()
cohort.age_table

cohort.ancestry()
cohort.ancestry_table

cohort.sex()
cohort.sex_table

cohort.mortality()
cohort.mortality_table

# cohort.omics_sample_metadata()
# cohort.omics_sample_data
# cohort.omics_sample_counts
# cohort.omics_sample_location
cohort.rd_sample_data()
cohort.rd_samples
cohort.platekeys

cohort.concat_all()
cohort.all_mortality
cohort.all_data
# cohort.all_omics_metadata.drop_duplicates(keep='first')
cohort.all_sample_counts
cohort.all_data  # 4770 partcipants


# test filtering out a few participants:
# cohort.custom_pids(pids_lst_file=[111001622,111001638,111002005, 111002048], action='exclude')
# now we can further filter the cohort.
# lets say for ICD-10 codes M81 and M80 and Q82 we only want those age 30 and below? 
# (for completeness sake we could do this by event dates - what if they got that diagnosis when they were 30, but are 40 at enrolment?)
filter = '''age_at_consent > 30'''
cohort.feature_filter(query_string=filter, action='exclude')

cohort.all_data # 1263 data point
len(cohort.all_data.participant_id.unique()) 
cohort.pids # 214 participants.

cohort.feature_tables

cohort.summary_stats() 

cohort.pids  # is this the only 
cohort.concat_all()
cohort.pids
cohort.all_pids
cohort.summary  # are the summary stats not updated after I apply a feature filter?
cohort.all_age.age_at_consent.value_counts()
cohort.all_age.loc[cohort.all_age['age_at_consent']>30]

cohort.all_pids.to_csv('/Users/christianbouwens/Documents/Projects/Relation_therapeutics/parts.tsv')

from gelpack.gel_vis import vis_cohorts, simple_count_plt, simple_hist_plt
import seaborn as sns
import matplotlib.pyplot as plt
plt.clear()
plt.clf()
fig = plt.figure()
vis_cohorts(cohort, figure = fig, mask = True, show=False)
fig.tight_layout()
fig.savefig('/Users/christianbouwens/Documents/Projects/Relation_therapeutics/young_osteo_summary.png', dpi=600)

cohort2 = {
	'terms':[
		"Osteogenesis imperfecta",
		],
	'icd10':[
		'Q780'
		]
		}
# imperfecta:
cohort2 = Cohort(
	featdict=cohort2,
	version=version,
	name='imperfecta'
	)
cohort2.get_term_pids()
cohort2.get_icd10_pids()
cohort2.age()
cohort2.mortality()
cohort2.ancestry()
cohort2.sex()
cohort2.concat_all()
cohort2.summary_stats()
cohort2.summary  
plt.clear()
plt.clf()
fig = plt.figure()
vis_cohorts(cohort2, figure = fig, mask = True, show=False)
fig.tight_layout()
fig.savefig('/Users/christianbouwens/Documents/Projects/Relation_therapeutics/imperfecta.png', dpi=600)




plt.clear()
plt.clf()
fig = plt.figure()
vis_cohorts([cohort,cohort2], figure = fig, mask = True, show=True)
fig.tight_layout()
fig.savefig('/Users/christianbouwens/Documents/Projects/Relation_therapeutics/combined_osteo_imperfecta.png', dpi=600)


comb = {
	'terms':[
		"Osteogenesis imperfecta",
		],
	'icd10':[
		'Q780',
		'M81',
		'M80',
		'M82'
		],
		}
# combined
comb_cohort = Cohort(
	featdict=comb,
	version=version,
	name='Osteoporsosis and Osteogensis Imperfecta'
	)
comb_cohort.get_term_pids()
comb_cohort.get_icd10_pids()
comb_cohort.age()
comb_cohort.mortality()
comb_cohort.ancestry()
comb_cohort.sex()
comb_cohort.concat_all()
comb_cohort.all_data
filter = ('diag.str.contains("M8") and age_at_consent > 30') # and age_at_consent > 30')
# comb_cohort.all_data.query(filter, engine='python')
comb_cohort.feature_filter(filter,action='exclude')
comb_cohort.summary_stats()
comb_cohort.summary

comb_cohort.all_age


cohort.summary
cohort2.summary
comb_cohort.summary




plt.clear()
plt.clf()
fig = plt.figure()
vis_cohorts(comb_cohort, figure = fig, mask = True, show=True)
fig.tight_layout()
fig.savefig('/Users/christianbouwens/Documents/Projects/Relation_therapeutics/added_osteo_imperfecta.png', dpi=600)




##### the gmc trust and or availability of omics data should be baked in the script now.

sqlstr = (f'''
	SELECT
		participant_id, 
		registered_at_gmc_trust,
		registered_at_gmc_ods_code
	FROM
		participant
	WHERE
		participant_id IN {*cohort.all_pids,} '''
)
pids_ods = lab_to_df(sqlstr, dr=version)
pids_ods = cohort.all_gmc_registration
pids_ods
full_reverse = pd.merge(cohort.all_omics_metadata, pids_ods, on='participant_id', how='left', indicator=True)

ods_names = (pids_ods[[
		'registered_at_gmc_ods_code','registered_at_gmc_trust'
		]]
	.drop_duplicates()
	.sort_values('registered_at_gmc_ods_code')
	.dropna()
	)

rna = full_reverse.loc[
	full_reverse['sample_type']=='RNA Blood'
	]
prot = full_reverse.loc[
	full_reverse['sample_type']=='Serum'
	]
edta_prot  = full_reverse.loc[
	full_reverse['sample_type']=='EDTA Plasma'
	]
pids_ods['rna_avail'] = np.where(pids_ods['participant_id'].isin(rna['participant_id']),'available', 'unavailable')
pids_ods['protein_avail'] = np.where(pids_ods['participant_id'].isin(prot['participant_id']),'available', 'unavailable')
pids_ods['edta_protein_avail'] = np.where(pids_ods['participant_id'].isin(edta_prot['participant_id']),'available', 'unavailable')

rnaperc = pids_ods['rna_avail'].value_counts(dropna=False, normalize=True)
protperc = pids_ods['protein_avail'].value_counts(dropna=False, normalize=True)
edtaprotperc = pids_ods['edta_protein_avail'].value_counts(dropna=False, normalize=True)

percentages = pd.concat([rnaperc, protperc, edtaprotperc],axis=1)
percentages.T



## plotting


plt.close()
plt.clf()
fig = plt.figure(figsize=(10,6))
gs = fig.add_gridspec(ncols=2, nrows=2)

axs1 = fig.add_subplot(gs[0,:])
axs2 = fig.add_subplot(gs[1,0])
axs3 = fig.add_subplot(gs[1,1])
simple_count_plt(
	table=full_reverse.drop_duplicates(
		['participant_id','registered_at_gmc_ods_code'],
		keep='first'),
	x='registered_at_gmc_ods_code',
	ax=axs1,
	mask=True)
axs1.set_ylabel('Participants')
axs1.set_xlabel(None)
axs1.set_title('GMC ODS registeration')


percentages.T.plot(
	kind='bar', 
	stacked=True,
	ax=axs2)
axs2.set_xticklabels(['RNA blood','Serum','Plasma'])
axs2.set_ylabel('fraction of cohort')
axs2.set_xlabel(None)
axs2.set_title('Availability',y=1.15)
sns.move_legend(axs2,
	"lower center", bbox_to_anchor=(.45, 1), ncol=3, title=None, frameon=False,)

tab = axs3.table(
	cellText=ods_names.values, 
	colLabels=['ODS code', 'GMC name'], 
	loc='center',
	cellLoc='center')
tab.auto_set_font_size(False)
tab.auto_set_column_width(col=[0,1])
tab.scale(1.5, 1)
tab.set_fontsize(7.5)
for i in range(0,len(ods_names.values)):
    tab[i, 0]._loc = 'center'
    tab[i, 0]._text.set_horizontalalignment('center') 

for i in range(0,len(ods_names.values)):
    tab[i, 1]._loc = 'left'
    tab[i, 1]._text.set_horizontalalignment('left') 
axs3.axis('off')
# axs3.set_xlabel=None
# axs3.set_ylabel=None
# axs3.set_axis_off=True
# axs3.set_frame_on=False

plt.tight_layout()
plt.savefig('/Users/christianbouwens/Documents/Projects/Relation_therapeutics/omics.png', dpi=300)

rna.participant_id.value_counts()


plt.close()
plt.clf()
fig = plt.figure()
gs = fig.add_gridspec(ncols=2,nrows=1)
axs1 = fig.add_subplot(gs[0,0])
axs2 = fig.add_subplot(gs[0,1])

bins = max(rna.aliquots)-min(rna.aliquots)
simple_hist_plt(
	table=rna.drop_duplicates(['participant_id'],keep='first'),
	x='aliquots',
	binwidth=1,
	ax=axs1,
	mask=True)
axs1.xaxis.set_label_position('bottom')
axs1.set_xlabel('Number of Aliquots')
axs1.set_title('RNA blood samples')

simple_hist_plt(
	table=prot.drop_duplicates(['participant_id'],keep='first'),
	x='aliquots',
	binwidth=1,
	ax=axs2,
	mask=True)
axs2.xaxis.set_label_position('bottom')
axs2.set_xlabel('Number of Aliquots')
axs2.set_title('Serum')
plt.tight_layout()
plt.savefig('/Users/christianbouwens/Documents/Projects/Relation_therapeutics/aliquots.png', dpi=400, width=8, height = 12)

# testc = Cohort(
# 	featdict=None, 
# 	version=version,
# 	name='test',
# 	participants = [])
# testc.pids.keys()
# testc.age()
# testc.sex()
# testc.mortality()
# testc.ancestry()

##### the issie is in concat_cohort()!
testc.concat_all()  
testc.all_data 
testc.summary_stats()
testc.summary
testc.pids.keys()
cohort.all_data
# filters:
# (subset,feature_filter,action(include/exclude))
filters = {
	'filter1':{
		'diag_subset':'(M81|M80|M82)',
		'feature_filter':'age_at_consent > 30',
		'action':'exclude'
	}
}

filter2 = '''
	diag.str.contains('Q78') and sex == 'Female'
	'''
# include / exclude
# we could assemble the query string based on this dictionary,
# or we just ask for the user to use pd.query format.


# we should have more here as osteoeneiss imperfecta is not filtered.
cohort.all_tables = cohort.all_data.loc[~(
	(cohort.all_data['age_at_consent'] > 30) & (cohort.all_data.diag.str.contains('(M81|M80|M82)'))
	)]

un_participants = filtered_master_table.participant_id.drop_duplicates()


un_participants_global30 = fifth_merge.loc[~(
	(fifth_merge['age_at_consent'] > 30)
	),'participant_id'].drop_duplicates()
# now recreate the cohort with these Pids. 

# easiest is filter current cohort. or create a new one.

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

from gelpack.gel_vis import vis_cohorts
fig2 = vis_cohorts(cohort, mask = False, show=True)
fig2.show()



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
        ]
	}
# this can get a bit murky - as histology-codes are tied to particular disease types.
# if these are not included we can assume we want all the histology types for a study abbreviation.
#     'abbr_histology_incl':[
#         ('LUG','8046/3')
#         ],
#     'abbr_histology_excl':[]
# }

# cohort_allca = Cohort(featdict=)

cohort2 = Cohort(
	featdict=cohort2_dict, 
	version=version,
	name='NSCLC')

cohort2.featdict

cohort2.get_icd10_pids(limit=['all'])
cohort2.icd10_table
cohort2.get_cterm_pids()
cohort2.cterm_table

cohort2.pids.keys()

cohort2.age()
cohort2.ancestry()
cohort2.sex()
cohort2.mortality()
cohort2.concat_all()

cohort2.ca_sample_data()
cohort2.platekeys  # these may currently still have platekeys where participants were not recruited for the specific disease type.
cohort2.limit_samples_to_cancer_type
# visualise the single cohort.
cohort2.sample_tables
cohort2.cancer_samples['icd10'].columns
cohort2.concat_all()
cohort2.all_data

import gel_vis
fig2 = vis_cohorts(cohort, show=False)
fig2.show()

vars(cohort)
# get all participant ids associated with 
# icd-10 -> we don't know if we have a sample and/or if this sample is for this specific icd-10 code
# cancer_term -> cancer_participant_disease: not that accurate, and we don't kow if we a participant has a sample.
# cancer_abbr -> only availabe in cancer_analysis. They are more accurate than cancer_terms + we have a sample.
#              -> we may want to specify a specfic histology to include / exclude for a cohort. 
#                   - are these only associated with ICD-10 and/ or study_abbreviations? 
#                   - can we do this based on cancer analysis or do we need to grab our secondary data resources?


# TODO what if I want to start a cohort from all cancer_analysis samples?

sqlstr = ('''
	SELECT 
		tumour_sample_platekey,
	FROM
		cancer_analysis
	WHERE
		study_abbreviation = 'BRCA'
	'''
)

ca = lab_to_df(sqlstr, version)

len(ca.participant_id.unique())

ca_cohort = Cohort(participants=ca.participant_id, version = version)
ca_cohort = Cohort(participants='all_cancer', version = version)
ca_cohort.ca_sample_data()
ca_cohort.cancer_samples




bca_cohort = Cohort(
	platekeys=list(ca['tumour_sample_platekey']), 
	version = version,
	name='breast_cancer')

bca_cohort.pids
bca_cohort.ca_sample_data()
bca_cohort.concat_all()
sqlstr = (
    f'''
    SELECT
        ca.participant_id,
        ca.tumour_sample_platekey,
        ca.disease_type,
        ca.tumour_clinical_sample_time,
        ca.tumour_type,
        av.diagnosisdatebest,
        av.er_status,
        av.pr_status,
        av.her2_status
    FROM lists.cancer_analysis ca
    INNER JOIN lists.av_tumour av
    ON ca.participant_id = av.participant_id
    AND av.diagnosisdatebest =
                (
                SELECT 
                    MAX(av_e.diagnosisdatebest)
                FROM 
                    lists.av_tumour av_e 
                WHERE 
                    av_e.participant_id = av.participant_id 
                    AND 
                    ca.tumour_clinical_sample_time >= av_e.diagnosisdatebest
                )
    WHERE ca.participant_id IN {*bca_cohort.all_pids,}
    ''')

len(bca_cohort.platekeys)

av_stat = lab_to_df(sqlstr, dr=version)

av_stat.her2_status.value_counts(dropna=False)
av_stat.er_status.value_counts(dropna=False)
av_stat.pr_status.value_counts(dropna=False)

av_stat.loc[
	(
		(av_stat['her2_status'].isin(['N','X']))
	& (av_stat['er_status'].isin(['P','Pm','B']) | av_stat['pr_status'].isin(['P','Pm','B']))
	)
]


from gelpack.cohort import Cohort
from gelpack.survival import cohort_surv, km_survival, logrank

version = "main-programme/main-programme_v17_2023-03-30"

lung_cohort = Cohort(
    featdict={'cancer_terms':['LUNG']},
    version=version,
    name='lung_cancer')
lung_cohort.get_cterm_pids()
lung_cohort.ca_sample_data()
lung_cohort.concat_all()
lung_cohort.all_data
lung_cohort.feature_filter(
    query_string=('nsv4_somatic_small_variants_annotation_vcf.isna()')
    )

colo_cohort = Cohort(
    featdict={'cancer_terms':['COLORECTAL']},
    version=version,
    name='colorectal cancer')
colo_cohort.get_cterm_pids()
colo_cohort.ca_sample_data()
colo_cohort.concat_all()
colo_cohort.all_data
colo_cohort.feature_filter(
    query_string=('nsv4_somatic_small_variants_annotation_vcf.isna()')
    )

print(
    f'''size of Lung cohort: {len(lung_cohort.all_pids)} and size of colorectal cohort: {len(colo_cohort.all_pids)}''')


lung_cohort.feature_filter(query_string='participant_id.isin([211002184, 211002690, 213000337, 217000638, 221000303])')
colo_cohort.feature_filter(query_string='participant_id.isin([211002184, 211002690, 213000337, 217000638, 221000303])')

lung_cohort.age()
lung_cohort.sex()
lung_cohort.ancestry()

colo_cohort.age()
colo_cohort.sex()
colo_cohort.ancestry()


cohort, na_counts, survival_data, mapping = cohort_surv(cohorts=[lung_cohort, colo_cohort])
cohort.ca_sample_data()
cohort.cancer_samples
cohort.concat_all()
cohort.all_cancer_samples
cohort.select_single_ca_sample()
cohort.sample_tables

cohort.all_cancer_samples
cohort.all_cancer_samples.columns
cohort.all_data.study_abbreviation.value_counts()
sig_cols = ['signature_'+str(i) for i in range(1,30)]

total_data = pd.merge(
	survival_data,
	cohort.all_cancer_samples[['participant_id','tmb']+sig_cols],
	how='left',
	on='participant_id'
)
## coxPH can't have NAs. how do we account for that?
use_dat = total_data[[
	'group',
	'participant_id',
	'survival',
	'status',
	'age_at_consent',
	'sex',
	'predicted_ancestry',
	'tmb'] 
	+ sig_cols]

use_dat[sig_cols] = use_dat[sig_cols].fillna(0)
use_dat.loc[use_dat['tmb'].isna()]
use_dat= use_dat.dropna()


univariate_df,cph = coxph_regression(
	survival_df=use_dat,
	time_col='survival',
	event_col='status',
	covariates=['age_at_consent','sex','predicted_ancestry','tmb','group']+sig_cols)


import matplotlib.pyplot as plt

fig = plt.figure()
cph.plot()
plt.show()

import matplotlib.pyplot as plt

fig = plt.figure()
cph.plot()
plt.show()

kmsurvival(
	data=survival_data,
	strata=mapping.values(),
	map_dict=mapping,
	output=None,
	plt_title='Comparing lung and colorectal cancer.',
	plotting='show',
	table=False
)

# link covariates to survival_data. we should be able to do it in cohort_surv. (when the cohorts have features.)


# IGSF8,() NBR1, SYVN1, TMEM129, SND1 _ -> autophaggy of MHC-I