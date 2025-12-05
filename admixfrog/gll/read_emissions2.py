import numpy as np
from scipy.optimize import minimize
from numba import njit
from math import lgamma, exp


@njit
def binom_pmf(O, N, p):
    res = np.power(p, O) * np.power(1.0 - p, N - O)
    for i, (o, n) in enumerate(zip(O, N)):
        if o > 0 and n > 1:
            res[i] *= exp(lgamma(n + 1) - lgamma(o + 1) - lgamma(n - o + 1))
    return res


@njit   #for debug
def p_reads_given_gt_gllmode(O, N, Pcont, c, error, n_obs):    
    """calculates P(O | G); probabilty of anc/derived reads given genotype
    per read group

    """
    n_gt = 3
    read_emissions = np.ones((n_obs, n_gt))
    for g in range(3):
        p = c * Pcont + (1 - c) * g / 2
        p = p * (1 - error) + (1 - p) * error
        read_emissions[:, g] = binom_pmf(O, N, p)

    return read_emissions

import numpy as np


"""
for N observations, assuming we have ref_err and alt_err, each with shape (N,2), examples are:
ref_err = [[err1.1, err1.2], [err2.1...], [Nan], ...]
alt_err = [[Nan], [Nan], [err3.1], ...]

we want to have L: [P(O|G=1), P(O|G=0), P(O|G=2)] for each observation 


example:
for ref_alt record, we have alleles [0,0, 1] if we have 2 refs and 1 alt
with errors [err_ref1, err_ref2, err_alt1]
or we have alleles [0, 0] with errors [err_ref1, err_ref2]


for g=0
P(read_1 | g=0) * P(read_2 | g=0) * P(read_3 | g=0) ... 
for g=1
P(read_1 | g=1) * P(read_2 | g=1) * ... 
for g=2
P(read_1 | g=2) * P(read_2 | g=2) * ... 

read_i = p if allele is alt
read_i = 1-p if allele is ref

"""


def gl_per_read_from_err(ref_err, alt_err, Pcont, c):
    ref_err = np.asarray(ref_err, dtype=float)
    alt_err = np.asarray(alt_err, dtype=float)
    if ref_err.ndim == 0:
        ref_err = np.array([]) if np.isnan(ref_err) else np.array([ref_err])
    if alt_err.ndim == 0:
        alt_err = np.array([]) if np.isnan(alt_err) else np.array([alt_err])
    ref_err = ref_err[~np.isnan(ref_err)]
    alt_err = alt_err[~np.isnan(alt_err)]

    n_ref = len(ref_err)
    n_alt = len(alt_err)

    if n_ref + n_alt == 0:
        raise ValueError("Both ref_err and alt_err are empty (only NaNs).")
    alleles = np.concatenate([
        np.zeros(n_ref, dtype=int),
        np.ones(n_alt, dtype=int),
    ])
    errors = np.concatenate([ref_err, alt_err])

    L = np.zeros(3, dtype=float)

    for g in range(3):
        p_true = c * Pcont + (1.0 - c) * g / 2.0  # same as original code
        p_read = p_true * (1.0 - errors) + (1.0 - p_true) * errors  #same as original code
        per_read = np.where(alleles == 1, p_read, 1.0 - p_read)
        L[g] = per_read.prod()
    return L


def gl_per_read_multi(ref_err_list, alt_err_list, Pcont, c):
    if len(ref_err_list) != len(alt_err_list):
        raise ValueError("ref_err_list and alt_err_list must have the same length")
    n_snp = len(ref_err_list)
    L_all = np.zeros((n_snp, 3), dtype=float)
    for i in range(n_snp):
        L_all[i, :] = gl_per_read_from_err(
            ref_err_list[i],
            alt_err_list[i],
            Pcont[i],
            c[i]
        )
    return L_all

def p_reads_given_gt(*args, gt_mode=False, **kwargs):
    if gt_mode:
        return p_reads_given_gt_gtmode(*args, **kwargs)
    else:
        return p_reads_given_gt_gllmode(*args, **kwargs)


@njit
def p_reads_given_gt_gtmode(O, N, Pcont, c, error, n_obs):
    """calculates P(O | G); probabilty of anc/derived genotype given input genotype"""
    n_gt = 3
    read_emissions = np.ones((n_obs, n_gt))
    for g in range(n_gt):
        read_emissions[O / N == g / 2, g] = 1 - 2 * error[O / N == g / 2]
        read_emissions[O / N != g / 2, g] = error[O / N != g / 2]

    return read_emissions


@njit
def read2snp_emissions(read_emissions, n_snps, ix):
    n_gt = 3
    snp_emissions = np.ones((n_snps, n_gt))
    for i, row in enumerate(ix):
        snp_emissions[row] *= read_emissions[i]
    return snp_emissions


def p_snps_given_gt(P, c, error, IX, gt_mode=False, position_based_error=False):
    """calculates probabilty of observed read data given genotype"""
    
    if position_based_error:
        read_emissions = gl_per_read_multi(P.ref_err, P.alt_err, P.P_cont, c)
    else:
        read_emissions = p_reads_given_gt(
        P.O, P.N, P.P_cont, c, error, IX.n_obs, gt_mode=gt_mode
    )
    return read2snp_emissions(read_emissions, IX.n_snps, IX.OBS2SNP)
