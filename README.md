# Survival Python

## GEL survival analysis
This python package contains functions to extract diagnosis dates, date of death and last follow ups with the clinic to calculate and plot a Kaplan-Meier survival analysis (stratifying on Domain variants).

## Data sources and approach:
### This package requires a set up .netrc to link with LabKey.

Cancer participants:
* PRIMARY cancers in cancer_analysis
* keeping unique participant_ids + disease_types.

Date of death:
* Earliest date between death_details and mortality tables

Date of diagnosis:
* Sources:
- cancer_participant_tumour
- av_tumour
- cancer_registry
* Approach:
- Translate ICD-10 code to a gel_disease_type
- Sort tables by date and keep only the first diagnosis date for each disease type.
- merge sources and keep unique values per participant.
- merge dod with cancer_analysis matching both the participant_id and disease_types

* inferring missing date of diagnosis from averages only occurs if the --infer flag is True.

Date of last follow up:
- the last known interaction within HES data.


## License
MIT license

