import numpy as np
import scipy.stats as sts

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2020 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'

def benjaminihochberg_correction(p, FDR=0.05):
    """


    Parameters
    ----------
    p : np.array
        array of p-values of independent tests.
    FDR : float, optional
        False discovey rate. 0 < FDR < 1. The default is 0.05.

    Returns
    -------
    passes_FDR
        A vector of same length as input parameter p, with bool value True
        where the test passed a BH FDR correction.

    """
    sorted_p = np.sort(p)
    rank = np.arange(1, len(p)+1)

    BH_crit = (rank/rank[-1])*FDR

    if np.any(sorted_p < BH_crit):
        thresh = sorted_p[sorted_p < BH_crit][-1]
        return p <= thresh
    else:
        return p < 0

def bonferroni_correction(p, FDR=0.05):
    """
    Bonferroni correction for multiple testing. Takes in a vector of p-values
    and returns true where p_i < alpha_corrected, where alpha_corrected = 0.05/N
    and N = number of independent tests.

    Parameters
    ----------
    p : list or numpy array
        Probabilities of independent statistical tests stored in a vector.

    Returns
    -------
    bonferroni_correction : numpy array
        A boolean vector of same shape as p where True indicates that the test
        passed a bonferroni correction.

    See Also
    --------
    stat_utils.benjaminihochberg_correction : often less stringent test correction

    """
    return np.array(p) < (FDR/len(p))


def _stringdb_bootstrap(summed_score, ppi, nTFs, FDR=0.05, N=100):
    unique_string_tfs = ppi.index.unique()
    for _ in range(N ):
        TFrand = np.random.choice(unique_string_tfs, size=nTFs, replace=False)
        random = ppi[ppi.index.isin(TFrand)].groupby('target_SYMBOL')['combined_score'].sum()
        summed_score[_] = random

    # Nan here just means not found, so should be 0
    summed_score[np.isnan(summed_score)] = 0

    means = summed_score.iloc[:, 1:].mean(1)
    stds = summed_score.iloc[:, 1:].std(1)
    Z = (means - summed_score.combined_score)/stds
    p = sts.norm.cdf(Z)
    is_sign = benjaminihochberg_correction(p, FDR=FDR)
    return is_sign

def _fisher_approx(a,b,c,d):
    """
    Estimates very small p-values in a Fisher exact test using logarithms.

    Parameters
    ----------
    a, b, c, d : int
        According to the table [[a, b], [c, d]]

    Returns
    -------
    p : float
        -log10 estimation of the p value.

    """
    def get_fact(num1, num2):
        return [i for i in range(2,num1+num2+1)]

    ab_fac = get_fact(a,b)
    cd_fac = get_fact(c,d)
    ac_fac = get_fact(a,c)
    bd_fac = get_fact(b,d)

    all_upper = []
    for tmp in [ab_fac, cd_fac, ac_fac, bd_fac]:
        for tmp_element in tmp:
            all_upper.append(tmp_element)

    all_lower = []
    for num in [a,b,c,d,a+b+c+d]:
        for i in range(2,num+1):
            all_lower.append(i)

    return (np.sum(np.log10(all_upper)) - np.sum(np.log10(all_lower)))
