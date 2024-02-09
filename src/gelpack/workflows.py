#### Christian J. Bouwens
#### BRSC team
#### ingest and clean outputs from GEL workflows
#### currently supporting: complex variant workflow
#### last update: 2023.09.20

import matplotlib.pyplot as plt
import pyranges as pr
import pandas as pd
import numpy as np

############################################
def ingest_cv_workflow(
    input_path, 
    gene_dict, 
    type='sv',
    maxsize=20e6, 
    svtype=['INS','DEL','DUP','BND','GAIN','LOSS'],
    cohort=None):
    """reads the output of the complex variant workflow, sets the column names
    and applies some basic filtering to the detected variants.

    Args:
        input_path (str): path to the SV or CNV mergedResult files.
        gene_dict (dictionary): a dictionary of relevant elements in the gene of
            interest
        type (str, optional): is the input a structural variant (sv) or
            copy number variant (cnv)?. Defaults to 'sv'.
        maxsize (int, optional): Set the maximum size of variants,
            variants above this size will be filtered out. Defaults to 20e6.
        svtype (list, optional): limits the variants to these structural
            variant types. Defaults to ['INS','DEL','DUP']. (under development)
        cohort (Class, optional): Add a gelpack Cohort to limit the found variants to the 
            samples and attach the output to the cohort. Defaults to None.

    Returns:
        _type_: _description_
    """
    # this function requires to run on the cluster if you want to normalise copy numbers.
    import pandas as pd
    import numpy as np
    import pyranges as pr
    if type=='sv':
        columns = [
            'platekey',
            'chr',
            'start',
            'variant_id',
            'REF',
            'ALT',
            'quality',
            'filter',
            'svtype',
            'size',
            'end',
            'PR',
            'SR']
    elif type == 'cnv':
        columns = [
            'platekey',
            'chr',
            'start',
            'variant_id',
            'REF',
            'ALT',
            'filt',
            'end',
            'svtype',
            'CN',
            'MCC',
            'GT']
    
    data = pd.read_csv(
        input_path, 
        sep='\t', 
        header=None, 
        names = columns)
    
    if type == 'cnv':
        data['size'] = data['end'] - data['start']
    else:
        data.loc[data['size']=='.','size'] = 1
        data.loc[data['end']=='.','end'] = pd.NA
        data['end'].fillna(data['start']+1,inplace=True)
        data['size'] = data['size'].astype(np.int64)
        data['size'] = abs(data['size'])

    # get a gene / element range to compare SV overlap.
    starts = gene_dict['start']
    ends = gene_dict['end']
    chr = gene_dict['chr']
    
    element_range = pr.PyRanges(
        chromosomes=chr,
        starts=starts,
        ends=ends
        )

    vars = data[['chr','start','end','platekey']].reset_index(drop=True)
    vars.columns = "Chromosome Start End Name".split()
    var_pr= pr.PyRanges(vars)
    var_overlap = var_pr.overlap(element_range)

    # filte variants based on overlap with the element of interest
    # And subsequently size
    # we may run into issues here when working with fusion breakpoints.
    
    if var_overlap.df.shape[0] > 0:
        data_filt = data.loc[
            (data['platekey'].isin(var_overlap.df['Name']))
            & (data['size'] <= maxsize)]
        if not cohort:
            return data_filt
        else:
            if type == 'sv':
                cohort.sv = data_filt
            else:
                cohort.cnv = data_filt
    else: 
        print('no overlapping data found')


def sv_lof_filt(sv_data, gene_start, gene_end):
    """Further applies SV filtering based on breakpoints. Genes that get 
    inverted/duplicated with two breakpoinst outside the 5'-3' region of a gene
    will be unlikely to be a loss of function variant. Here we apply a filter for
    BNDs, and DUP svtpes. (note manta no longer returns INV in dragen V3.2).

    Args:
        sv_data (pd.Dataframe): a pandas dataframe filtered for variants that 
            overlap at least one exon of the gene of interest. 
        gene_start (int): 5' start site (genomic coordinates) of the gene of interest.
        gene_end (int): 3' end site (genomic coordinates) of the gene of interest

    Returns:
        pd.dataframe: a pandas dataframe with additional filtering applied for 
            variants that are likely affecting gene function.
    """
    sv_data['start'] = sv_data['start'].astype(np.int64)
    sv_data['end'] = sv_data['end'].astype(np.int64)
    sv_filt_affect_gene = sv_data.loc[~(
            (
                sv_data['svtype'].isin(['DUP'])
                & ((sv_data['start'] <= gene_start) & (sv_data['end'] >= gene_end))
            )
            |
            (
                sv_data['svtype'].isin(['BND'])
                        & (
                            (sv_data['start'] < gene_start)
                            | (sv_data['start'] > gene_end)
                            )
            )
        )]
    return sv_filt_affect_gene