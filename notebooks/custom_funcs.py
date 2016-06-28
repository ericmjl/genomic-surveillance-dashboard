import numpy as np
import pandas as pd
from Bio import SeqIO
from molecular_weight import molecular_weights
from isoelectric_point import isoelectric_points
from sklearn.cross_validation import train_test_split

allowed_drugnames = ['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV',
                     '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'EFV', 'NVP',
                     'ETR', 'RPV',
                     ]


def read_data(protein):
    """
    Reads in the data for the protein.
    """
    drug_col_vals = {'protease': 8,
                     'nnrt': 4,
                     'nrt': 6}

    assert protein in drug_col_vals.keys()

    data = pd.read_csv('../data/hiv-{0}-data.csv'.format(protein),
                       index_col='SeqID')

    drug_cols = data.columns[0:drug_col_vals[protein]]
    feat_cols = data.columns[drug_col_vals[protein]:]

    return data, drug_cols, feat_cols


def replace_ambiguous_letters_with_nan(df, feat_cols):
    """
    Within the dataframe `df`, replaces all cells amongst the feature columns
    that have ambiguous letters.

    Assumes that `df` has the structure:
    - Index = some sequence identifier (must be unique)
    - Columns 1-->N: regressable numbers.
    - Columns N-->end: sequence feature matrix.

    Inputs:
    =======
    - df: The `pandas` DataFrame that stores the HIV data.
    - feat_cols: The columns that correspond to the amino acid positions.

    Returns:
    ========
    - new_df: The `pandas` DataFrame that has all double or triple letter
              cells replaced with `np.nan`.
    """
    assert isinstance(df, pd.DataFrame)

    new_df = df.copy()
    for col in feat_cols:
        new_df[col] = df[col].apply(lambda x: np.nan if len(str(x)) > 1 else x)
    return new_df


def replace_dashes_with_canonical_letters(df, consensus_map, feat_cols):
    """
    Replaces all the dashes with the canonical letter for that string.

    Cleans the data by:
    - imputing consensus amino acids into the sequence.
    - removing sequences with early stop codons.
    - removing sequences with ambiguous sequences (tough to deal with).

    Parameters:
    ===========
    - data: (pandas DataFrame) the sequence feature matrix with drug
            resistance metadata.
    - feat_cols: (pandas columns) the drug resistance columns.
    - consensus_map: (dict) a dictionary going from position to consensus
                     sequence.

    Returns:
    ========
    - data: (pandas DataFrame) the cleaned sequence feature matrix with drug
            resistance metadata.
    """

    # Defensive programming checks
    assert isinstance(df, pd.DataFrame)
    assert isinstance(consensus_map, dict)

    # Impute consensus sequence
    for i, col in enumerate(feat_cols):
        # Replace '-' with the consensus letter.
        df[col] = df[col].replace({'-': consensus_map[i]})

    df = df.replace({'X': np.nan})
    df = df.replace({'.': np.nan})

    return df


def read_consensus(protein_name):
    """
    Reads in the consensus sequence, makes a map of position to letter.
    """

    # Defensive programming checks.
    allowed_names = ['protease', 'nnrt', 'nrt']
    assert protein_name in allowed_names,\
        "protein_name must be in {0}".format(allowed_names)

    if protein_name == 'protease':
        handle = '../data/hiv-protease-consensus.fasta'
    elif protein_name in ['nnrt', 'nrt']:
        handle = '../data/hiv-rt-consensus.fasta'
    consensus = SeqIO.read(handle, 'fasta')
    consensus_map = {i: letter for i, letter in enumerate(str(consensus.seq))}

    return consensus_map


def drop_na_from_data(df, drug_name, feat_cols):

    # Defensive programming checks
    assert isinstance(df, pd.DataFrame)
    assert isinstance(drug_name, str)
    assert drug_name in allowed_drugnames
    assert drug_name in df.columns, "{0} not in data.".format(drug_name)

    columns = list(feat_cols)
    columns.append(drug_name)

    return df[columns].dropna()


def get_cleaned_data(protein_name, drug_name):
    """
    A composition of the above functions. Expresses the "business logic"
    behind the data preprocessing steps.
    """
    data, drug_cols, feat_cols = read_data(protein_name)
    consensus_map = read_consensus(protein_name)
    data = replace_ambiguous_letters_with_nan(data, feat_cols)
    data = replace_dashes_with_canonical_letters(data,
                                                 consensus_map,
                                                 feat_cols)

    data = drop_na_from_data(data, drug_name, feat_cols)
    data[drug_name] = data[drug_name].apply(lambda x: np.log10(x))

    return data, feat_cols


def to_numeric_rep(df, feat_cols, rep='mw'):
    if rep == 'mw':
        numeric_dict = molecular_weights
    if rep == 'pKa':
        numeric_dict = isoelectric_points
    for col in feat_cols:
        df[col] = df[col].replace(numeric_dict.keys(),
                                  numeric_dict.values())

    return df


def test_data_integrity(data):
    """
    Data integrity tests! Make sure to use this to check that the data are
    suitable for input into scikit-learn.
    """
    assert '-' not in data
    assert '.' not in data
    assert 'X' not in data
    assert not data.dropna().isnull().values.any()
