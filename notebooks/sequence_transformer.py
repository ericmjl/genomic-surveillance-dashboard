"""
Functions that transform a sequence into its numerical representation.
"""
from Bio import SeqIO
from molecular_weight import molecular_weights
from isoelectric_point import isoelectric_points
from scipy.interpolate import interp1d

import pandas as pd
import numpy as np

reflengths = dict()
reflengths['protease'] = 99
reflengths['rt'] = 560


def to_numeric_rep(sequence, rep='mw'):
    """
    Converts a linear sequence to its molecular weight representation.

    Parameters:
    ===========
    - sequence: (str) the amino acid string.
    - rep: (str) one of ['mw', 'pKa']
    """
    assert isinstance(sequence, str), 'sequence must be a string.'
    assert 'X' not in sequence, 'X not allowed.'
    allowed_rep = ['mw', 'pka']
    assert rep in allowed_rep, 'rep must be one of {0}'.format(allowed_rep)

    if rep == 'mw':
        numeric_dict = molecular_weights
    if rep == 'pKa':
        numeric_dict = isoelectric_points

    df = pd.DataFrame(np.array([i for i in sequence]))\
        .T\
        .replace(numeric_dict.keys(), numeric_dict.values())

    return df.ix[0].values


def standardize_sequence(rep, protein='protease'):
    """
    Standardizes the sequence to a particular protein's length.

    Parameters:
    ===========
    - rep: (np.array) the n-variable length array.
    - protein: (str) one of ['protease', 'rt']
    """
    # assert isinstance(rep, np.array), 'rep must be a numpy array'
    assert isinstance(protein, str), 'protein must be a string'
    assert protein in reflengths.keys(), 'protein must be one of {0}'.format(
        reflengths.keys())

    interpolator = interp1d(np.arange(rep.size), rep, fill_value="extrapolate")
    ref_size = reflengths[protein]
    interp_arr = interpolator(np.linspace(0, ref_size, ref_size))
    return interp_arr
