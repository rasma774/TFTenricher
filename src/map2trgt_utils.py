import pandas as pd
import numpy as np

np.random.seed(0)

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2020 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'

def trrust_genes(TFs, weighted=False, silent=False, top_n_genes=None):
    """


    Parameters
    ----------
    TFs : list
        List of transcription factors to map to target genes.

    Returns
    -------
    List of TRRUST target genes.

    """
    # Load TRRUST
    pw = __file__.split('/src')[0]

    if not silent:
        print('loading TRRUST')
    TRRUST = pd.read_csv(
        pw + '/data/TRRUST/trrust_rawdata.human.tsv',
        sep='\t',
        header=None,
        )
    if not silent:
        print('Done')

    # As of now, we dont use the direction or publications
    TRRUST = TRRUST.iloc[:, :2]


    in_TFs = np.in1d(TRRUST[0].unique(), TFs)
    in_TRRUST = np.in1d(TFs, TRRUST[0].unique())


    if not silent:
        print(str(100*np.sum(~in_TRRUST)/len(in_TRRUST)) + '% of TFs are not in TRRUST')
        print(str(100*np.sum(in_TFs)/len(in_TFs))[:5] + '% of TRRUST TFs were in the TF list')

    target_genes = TRRUST[TRRUST[0].isin(TFs)][1].values
    unique_targets, counts = np.unique(target_genes, return_counts=True)

    if ~weighted:
        return unique_targets

    targets = pd.DataFrame((unique_targets, counts)).transpose()
    targets = targets.set_index(targets.columns[0]).iloc[:, 0]
    return targets


def correlation_genes(TFs, thresh=0.95, silent=False, top_n_genes=None):
    """


    Parameters
    ----------
    TFs : list
        List of transcription factors to map to target genes..

    Returns
    -------
    Panda series of correlating target genes summed over TFs.

    """

    pw = __file__.split('/src')[0]
    if not silent:
        print('loading corr')
    corr = pd.read_pickle(pw + '/data/pickles/correlations.p')
    if not silent:
        print('Done')


    in_TFs = corr.index.isin(TFs)
    in_corr = np.in1d(TFs, corr.index)

    if in_corr.sum() == 0:
        raise Exception('No input TFs are in correlation matrix')


    if not silent:
        print(str(100*np.sum(~in_corr)/len(in_corr)) + '% of TFs are not found')
        print(str(100*np.sum(in_TFs)/len(in_TFs))[:5] + '% of correlation table TFs were in the TF list')


    # Since the self-correlation is one, we need to remove the input TFs from
    # the set
    corr_out = corr.iloc[:, ~corr.columns.isin(TFs)]
    target_genes = corr_out[in_TFs].abs().sum()

    if top_n_genes is not None:
        return target_genes.sort_values()[::-1][:top_n_genes].index

    if thresh == -1:
        return target_genes.index.values

    cval_dist = []
    for _ in range(40):
        randtfs = np.random.choice(corr.index, size=(in_corr.sum()), replace=False)
        ctmp = corr.iloc[:, ~corr.columns.isin(randtfs)][corr.index.isin(randtfs)].abs().sum()
        cval_dist.append(np.sort(ctmp)[int(len(ctmp)*thresh)])

    target_genes_adj = target_genes[target_genes >= np.max(cval_dist)]
    return target_genes_adj.index.values



def STRING_ppi(TFs, FDR=0.95, Npermut=100, silent=False, top_n_genes=None):
    """


    Parameters
    ----------
    TFs : TYPE
        DESCRIPTION.
    multiple_testing_correct : TYPE, optional
        DESCRIPTION. The default is True.
    thresh : TYPE, optional
        DESCRIPTION. The default is 0.95.

    Returns
    -------
    None.

    """
    pw = __file__.split('/src')[0]
    if not silent:
        print('loading STRING PPI...')
    ppi = pd.read_pickle(pw + '/data/string_links.p')
    if not silent:
        print('Done')

    in_TFs = np.in1d(ppi.index.unique(), TFs)
    in_STRINGdb = np.in1d(TFs, ppi.index.unique())

    if not silent:
        print(str(100*np.sum(~in_STRINGdb)/len(in_STRINGdb)) + '% of TFs are not in STRINGdb')
        print(str(100*np.sum(in_TFs)/len(in_TFs))[:5] + '% of STRINGdb TFs were in the TF list')

    summed_score = ppi[ppi.index.isin(TFs)].groupby('target_SYMBOL')['combined_score'].sum()
    summed_score = pd.DataFrame(summed_score)

    if top_n_genes is not None:
        return summed_score.sort_values('combined_score')[::-1][:top_n_genes].index

    # TODO: The test here is not stringent at all. The more TFs, the more power,
    # and the more targets we get. Tested random 400 TFs, got 6500 significant genes.
    # But may be reasonable... ~25% of all possible TFs gave ~25% of all genes.
    # Should add option to just select top N genes

    if FDR == -1:
        return summed_score
    raise Exception('STRINGdb without no top_n_genes not supported')


