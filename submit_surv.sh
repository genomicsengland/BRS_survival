#!/bin/bash
#BSUB -q short
#BSUB -P bio
#BSUB -o ./survival.%J.out
#BSUB -e ./survival.%J.err
#BSUB -J survival
#BSUB -R "rusage[mem=10000] span[hosts=1]"
#BSUB -M 8000
#BSUB -n 1
#BSUB -cwd .module purge

source /resources/conda/miniconda3/envs/pd_lk_survival

python3 survival.py \
--genes KRAS TP53 \
-s and \
--imputate \
-o BRS_test/




# 3 use cases: 
# 1. run everything from the script e.g.: choose disease types (or pan cancer) and domain variants to stratify on.
#       output: survival table (median, events, 9% confidence interval) + KM figure
# 2. for selected participant ids and stratifications - grab survival time and censoring/death
#       output - surv_dat[['participant_id', survival time, censoring]] - to merge with own events.
# 3. just run the KM or CoxPH if you have a table with surv/event/strat
#       output same as 1 without the pre-processing.



