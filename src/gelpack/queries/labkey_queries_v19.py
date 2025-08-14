# dictionary of required SQL queries.
# this file needs to be edited when tables change during data releases.

# Centralised SQL templates used by Cohort + Survdat.

# Placeholders:
#   {participants} / {platekeys} / {disease_types} / {abbrs} – must be IN (...) strings
#   {like_diag} – an OR-joined list of LIKE clauses against the canonical alias "diag",
#                 e.g. "(diag LIKE '%C34%' OR diag LIKE '%D05%')"
#
# Notes:
# - Every template that needs a “LIKE” now exposes a canonical column "diag".
#	this is to improve the experience when updating column names / tables.
# - any stripping of dots (for instance for ICD10 codes), ishandled inside the CTE when creating "diag".
# - If an upstream table has multiple diagnosis columns, they’re UNIONed into one "diag".
#
# the loader should render:
#   sql = QUERIES[key].format(participants=in_clause, like_diag=or_like("diag", terms), ...)
# and can batch the {participants} placeholder if needed.

QUERIES = {
    # --------------------- Generic / utility ---------------------
    "cancer_parts": """
        SELECT participant_id
        FROM cancer_analysis
    """,

    "pids_by_platekeys": """
        SELECT participant_id, plate_key
        FROM plated_sample
        WHERE plate_key IN {platekeys}
    """,

    # --------------------- Rare disease / HPO / Cancer terms ---------------------
    "rd_participant_disease_by_terms": """
        SELECT DISTINCT participant_id, normalised_specific_disease
        FROM rare_diseases_participant_disease
        WHERE normalised_specific_disease IN {terms}
    """,

    "hpo_pids_by_ids": """
        SELECT DISTINCT participant_id, normalised_hpo_id, normalised_hpo_term
        FROM rare_diseases_participant_phenotype
        WHERE hpo_present = 'Yes'
          AND normalised_hpo_id IN {hpo_ids}
    """,

    "cancer_participant_disease_by_types": """
        SELECT DISTINCT participant_id, cancer_disease_type AS disease_type
        FROM cancer_participant_disease
        WHERE cancer_disease_type IN {cancer_types}
    """,

    "cancer_analysis_by_disease_type": """
        SELECT DISTINCT participant_id, disease_type
        FROM cancer_analysis
        WHERE disease_type IN {cancer_types}
    """,

    "cancer_analysis_by_abbr": """
        SELECT DISTINCT participant_id, study_abbreviation
        FROM cancer_analysis
        WHERE study_abbreviation IN {abbrs}
    """,

    # --------------------- ICD-10 (HES / Mort / Mental Health / Cancer sources) ---------------------
    # HES: apc / ae / op – expose canonical diag + event_date (dot-stripped)
    "hes_apc_diag_like": """
        WITH base AS (
            SELECT DISTINCT
                participant_id,
                REGEXP_REPLACE(diag_all, '\\.', '') AS diag,
                disdate AS event_date
            FROM hes_apc
            WHERE participant_id IN {participants}
        )
        SELECT participant_id, diag, event_date
        FROM base
        WHERE {like_diag}
    """,

    "hes_ae_diag_like": """
        WITH base AS (
            SELECT DISTINCT
                participant_id,
                REGEXP_REPLACE(diag_all, '\\.', '') AS diag,
                arrivaldate AS event_date
            FROM hes_ae
            WHERE participant_id IN {participants}
        )
        SELECT participant_id, diag, event_date
        FROM base
        WHERE {like_diag}
    """,

    "hes_op_diag_like": """
        WITH base AS (
            SELECT DISTINCT
                participant_id,
                REGEXP_REPLACE(diag_all, '\\.', '') AS diag,
                apptdate AS event_date
            FROM hes_op
            WHERE participant_id IN {participants}
        )
        SELECT participant_id, diag, event_date
        FROM base
        WHERE {like_diag}
    """,

    # Mortality (ICD10 multiple cause) – canonical diag + event_date
    "mortality_icd10_like": """
        WITH base AS (
            SELECT DISTINCT
                participant_id,
                REGEXP_REPLACE(icd10_multiple_cause_all, '\\.', '') AS diag,
                date_of_death AS event_date
            FROM mortality
            WHERE participant_id IN {participants}
        )
        SELECT participant_id, diag, event_date
        FROM base
        WHERE {like_diag}
    """,

    # Mental health – collapse primary + secondary into one canonical "diag"
    "mhmd_event_icd10_like": """
        WITH base AS (
            SELECT participant_id, REGEXP_REPLACE(ic_eve_primarydiagnosis,   '\\.', '') AS diag
            FROM mhmd_v4_event
            WHERE participant_id IN {participants}
            UNION ALL
            SELECT participant_id, REGEXP_REPLACE(ic_eve_secondarydiagnosis, '\\.', '') AS diag
            FROM mhmd_v4_event
            WHERE participant_id IN {participants}
        )
        SELECT DISTINCT participant_id, diag
        FROM base
        WHERE {like_diag}
    """,

    "mhldds_event_icd10_like": """
        WITH base AS (
            SELECT participant_id, REGEXP_REPLACE(primarydiagnosis,   '\\.', '') AS diag
            FROM mhldds_event
            WHERE participant_id IN {participants}
            UNION ALL
            SELECT participant_id, REGEXP_REPLACE(secondarydiagnosis, '\\.', '') AS diag
            FROM mhldds_event
            WHERE participant_id IN {participants}
        )
        SELECT DISTINCT participant_id, diag
        FROM base
        WHERE {like_diag}
    """,

    # MHSDS family – keep diagscheme filters; expose canonical "diag"
    "mhsds_medical_history_previous_diagnosis_icd10": """
        WITH base AS (
            SELECT DISTINCT
                participant_id,
                REGEXP_REPLACE(prevdiag, '\\.', '') AS diag,
                diagschemeinuse AS diagscheme
            FROM mhsds_medical_history_previous_diagnosis
            WHERE participant_id IN {participants}
        )
        SELECT participant_id, diag
        FROM base
        WHERE diagscheme IN ('2','02')
          AND {like_diag}
    """,

    "mhsds_provisional_diagnosis_icd10": """
        WITH base AS (
            SELECT DISTINCT
                participant_id,
                REGEXP_REPLACE(provdiag, '\\.', '') AS diag,
                diagschemeinuse AS diagscheme
            FROM mhsds_provisional_diagnosis
            WHERE participant_id IN {participants}
        )
        SELECT participant_id, diag
        FROM base
        WHERE diagscheme IN ('2','02')
          AND {like_diag}
    """,

    "mhsds_primary_diagnosis_icd10": """
        WITH base AS (
            SELECT DISTINCT
                participant_id,
                REGEXP_REPLACE(primdiag, '\\.', '') AS diag,
                diagschemeinuse AS diagscheme
            FROM mhsds_primary_diagnosis
            WHERE participant_id IN {participants}
        )
        SELECT participant_id, diag
        FROM base
        WHERE diagscheme IN ('2','02')
          AND {like_diag}
    """,

    "mhsds_secondary_diagnosis_icd10": """
        WITH base AS (
            SELECT DISTINCT
                participant_id,
                REGEXP_REPLACE(secdiag, '\\.', '') AS diag,
                diagschemeinuse AS diagscheme
            FROM mhsds_secondary_diagnosis
            WHERE participant_id IN {participants}
        )
        SELECT participant_id, diag
        FROM base
        WHERE diagscheme IN ('2','02')
          AND {like_diag}
    """,

    "mhsds_care_activity_icd10": """
        WITH base AS (
            SELECT DISTINCT
                participant_id,
                REGEXP_REPLACE(codefind, '\\.', '') AS diag,
                findschemeinuse AS diagscheme
            FROM mhsds_care_activity
            WHERE participant_id IN {participants}
        )
        SELECT participant_id, diag
        FROM base
        WHERE diagscheme IN ('1','01')
          AND {like_diag}
    """,

    "mhsds_indirect_activity_icd10": """
        WITH base AS (
            SELECT DISTINCT
                participant_id,
                REGEXP_REPLACE(codefind, '\\.', '') AS diag,
                findschemeinuse AS diagscheme
            FROM mhsds_indirect_activity
            WHERE participant_id IN {participants}
        )
        SELECT participant_id, diag
        FROM base
        WHERE diagscheme IN ('1','01')
          AND {like_diag}
    """,

    # Cancer-source ICD pattern – make {code_col} canonical as "diag"
    "cancer_icd10_from_table": """
        WITH base AS (
            SELECT DISTINCT
                participant_id,
                REGEXP_REPLACE({code_col}, '\\.', '') AS diag
            FROM {table}
            WHERE participant_id IN {participants}
        )
        SELECT participant_id, diag
        FROM base
        WHERE {like_diag}
    """,

    # --------------------- Age / Sex / Ancestry ---------------------
    "participant_age": """
        SELECT participant_id, year_of_birth, date_of_consent
        FROM participant
        WHERE participant_id IN {participants}
    """,

    "aggregate_gvcf_ancestry": """
        SELECT
            participant_id,
            pred_african_ancestries,
            pred_south_asian_ancestries,
            pred_east_asian_ancestries,
            pred_european_ancestries,
            pred_american_ancestries
        FROM aggregate_gvcf_sample_stats
        WHERE participant_id IN {participants}
    """,

    "participant_sex": """
        SELECT participant_id, participant_phenotypic_sex AS sex
        FROM participant
        WHERE participant_id IN {participants}
    """,

    # --------------------- Mortality (status / DoD) ---------------------
    "death_details_by_pids": """
        SELECT DISTINCT participant_id, death_date
        FROM death_details
        WHERE participant_id IN {participants}
    """,

    "mortality_by_pids": """
        SELECT participant_id, date_of_death
        FROM mortality
        WHERE participant_id IN {participants}
    """,

    "rare_diseases_pedigree_member_dod": """
        SELECT participant_id, date_of_death
        FROM rare_diseases_pedigree_member
        WHERE participant_id IN {participants}
    """,

    # --------------------- Rare Disease samples ---------------------
    "rd_samples_by_platekeys": """
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
        FROM rare_disease_interpreted
        WHERE plate_key IN {platekeys}
    """,

    "rd_samples_by_participants": """
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
        FROM rare_disease_interpreted
        WHERE participant_id IN {participants}
    """,

    # --------------------- Cancer Analysis samples ---------------------
    "cancer_samples_by_platekeys": """
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
            ca.tumour_clinical_sample_time,
            ca.preparation_method,
            ca.tumour_purity,
            ca.coverage_homogeneity,
            ca.somatic_coding_variants_per_mb AS tmb,
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
        FROM cancer_analysis ca
        LEFT JOIN cancer_100K_genomes_realigned_on_pipeline_2 da
            ON ca.tumour_sample_platekey = da.tumour_sample_platekey
        LEFT JOIN cancer_staging_consolidated csc
            ON ca.tumour_sample_platekey = csc.tumour_sample_platekey
        WHERE ca.tumour_sample_platekey IN {platekeys}
    """,

    "cancer_samples_by_participants": """
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
            ca.tumour_clinical_sample_time,
            ca.preparation_method,
            ca.tumour_purity,
            ca.coverage_homogeneity,
            ca.somatic_coding_variants_per_mb AS tmb,
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
        FROM cancer_analysis ca
        LEFT JOIN cancer_100K_genomes_realigned_on_pipeline_2 da
            ON ca.tumour_sample_platekey = da.tumour_sample_platekey
        LEFT JOIN cancer_staging_consolidated csc
            ON ca.tumour_sample_platekey = csc.tumour_sample_platekey
        WHERE ca.participant_id IN {participants}
    """,

    # --------------------- Omics sample metadata + GMC registration ---------------------
    "omics_metadata_by_participants": """
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
        FROM laboratory_sample samp
        LEFT JOIN laboratory_sample_omics_availability omics
            ON samp.laboratory_sample_id = omics.laboratory_sample_id
        WHERE samp.participant_id IN {participants}
    """,

    "participant_gmc_registration": """
        SELECT participant_id, registered_at_gmc_trust, registered_at_gmc_ods_code
        FROM participant
        WHERE participant_id IN {participants}
    """,

    # --------------------- Survdat: HES last seen + AV treatment ---------------------
    "hes_apc_last_seen": """
        SELECT participant_id, MAX(epistart) AS last_seen
        FROM hes_apc
        WHERE participant_id IN {participants}
        GROUP BY participant_id
    """,
    "hes_ae_last_seen": """
        SELECT participant_id, MAX(arrivaldate) AS last_seen
        FROM hes_ae
        WHERE participant_id IN {participants}
        GROUP BY participant_id
    """,
    "hes_op_last_seen": """
        SELECT participant_id, MAX(apptdate) AS last_seen
        FROM hes_op
        WHERE participant_id IN {participants}
        GROUP BY participant_id
    """,
    "hes_cc_last_seen": """
        SELECT participant_id, MAX(ccdisdate) AS last_seen
        FROM hes_cc
        WHERE participant_id IN {participants}
        GROUP BY participant_id
    """,
    "av_treatment_last_seen": """
        SELECT participant_id, MAX(eventdate) AS last_seen
        FROM av_treatment
        WHERE participant_id IN {participants}
        GROUP BY participant_id
    """,

    # --------------------- Survdat: diagnosis (DoD sources) ---------------------
    "cpt_diagnosis_by_pids": """
        SELECT DISTINCT
            participant_id,
            diagnosis_date,
            diagnosis_icd_code,
            morphology_icd_code
        FROM cancer_participant_tumour
        WHERE participant_id IN {participants}
    """,
    "av_tumour_diagnosis_by_pids": """
        SELECT DISTINCT
            participant_id,
            diagnosisdatebest,
            site_icd10_O2_3char,
            site_coded_3char,
            histology_coded
        FROM av_tumour
        WHERE participant_id IN {participants}
    """,
    "cancer_registry_by_pids": """
        SELECT DISTINCT
            participant_id,
            event_date,
            cancer_site,
            cancer_type,
            cancer_behaviour
        FROM cancer_registry
        WHERE participant_id IN {participants}
    """,

    # --------------------- Survdat: general CA pulls ---------------------
    "cancer_analysis_minimal": """
        SELECT
            participant_id,
            tumour_sample_platekey,
            study_abbreviation,
            disease_type,
            tumour_clinical_sample_time
        FROM cancer_analysis
    """,

    # CLI option for primary tumours by disease_type(s)
    "ca_primary_participants_by_disease_types": """
        SELECT participant_id, disease_type
        FROM cancer_analysis
        WHERE tumour_type = 'PRIMARY'
          AND disease_type IN {disease_types}
    """,
    "ca_primary_participants_all": """
        SELECT participant_id, disease_type
        FROM cancer_analysis
        WHERE tumour_type = 'PRIMARY'
    """,
}
