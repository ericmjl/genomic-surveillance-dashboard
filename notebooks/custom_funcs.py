import numpy as np
import pandas as pd
from Bio import SeqIO
from gsdash.molecular_weight import molecular_weights
from gsdash.isoelectric_point import isoelectric_points
from sklearn.cross_validation import train_test_split

allowed_drugnames = ['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV',
                     '3TC', 'ABC', 'AZT', 'D4T', 'DDI', 'TDF', 'EFV', 'NVP',
                     'ETR', 'RPV',
                     ]
drug_col_vals = {'protease': 8,
                 'nnrt': 4,
                 'nrt': 6}


def read_data(protein, sparse=True):
    """
    Reads in the data for the protein.

    Has two options:
    - sparse=True: loads the version that has dashes replacing consensus
                   sequence.
    - sparse=False: loads version has actual letters in each position from
                    consensus sequence.

    Returns:
    ========
    - data: (pd.DataFrame) protein sequence matched with drug resistance data
            (for all drugs)
    - drug_cols: (iterable) list of drug columns names.
    - feat_cols: (iterable) list of feature column names (corresponding to
                 sequence position)
    """
    assert protein in drug_col_vals.keys()

    if sparse:
        path = 'data/hiv-{0}-data-sparse.csv'.format(protein)
        sep = '\t'
    else:
        path = 'data/hiv-{0}-data.csv'.format(protein)
        sep = ','

    data = pd.read_csv(path, index_col='SeqID', sep=sep)

    drug_cols = data.columns[0:drug_col_vals[protein]]
    feat_cols = data.columns[drug_col_vals[protein]:]
    for col in feat_cols:
        data[col] = data[col].str.upper()

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

    def is_allowed(x):
        if x not in ['#', '~']:
            return True
        else:
            return False

    new_df = df.copy()
    for col in feat_cols:
        new_df[col] = df[col].apply(lambda x: np.nan if len(str(x)) > 1 else x)
        new_df[col] = new_df[col].apply(lambda x: np.nan
                                        if not is_allowed(x) else x)
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


def read_consensus(drug_class):
    """
    Reads in the consensus sequence, makes a map of position to letter.
    """

    # Defensive programming checks.
    assert drug_class in drug_col_vals.keys(),\
        "drug_class must be in {0}".format(drug_col_vals.keys())

    if drug_class == 'protease':
        handle = 'data/hiv-protease-consensus.fasta'
    elif drug_class in ['nnrt', 'nrt']:
        handle = 'data/hiv-rt-consensus.fasta'
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


def get_protein_drug_data(drug_class, sparse=False):
    """
    A composition of the above functions.

    Returns a dataframe that may contain NaN cells. Contains all of the drug
    resistance data (floating point columns) matched with the sequences. Drug
    resistance columns are not log10 transformed. Also returns the drug
    columns and feature columns.
    """
    data, drug_cols, feat_cols = read_data(drug_class, sparse=sparse)
    consensus_map = read_consensus(drug_class)
    data = replace_ambiguous_letters_with_nan(data, feat_cols)
    data = replace_dashes_with_canonical_letters(data,
                                                 consensus_map,
                                                 feat_cols)

    return data, drug_cols, feat_cols


def sparse_to_dense_data(drug_class):
    """
    A function that converts the sparsely-represented data to a dense form.
    """
    data, drug_cols, feat_cols = get_protein_drug_data(drug_class, sparse=True)
    data.to_csv('data/hiv-{0}-data.csv'.format(drug_class))


def get_cleaned_data(drug_class, drug_name):
    """
    A composition of the above functions. Expresses the "business logic"
    behind the data preprocessing steps.

    Returns a clean dataframe with no NaN cells. Positions 1 to (n-1) are the
    amino acids, while the last position is a single column of floating point
    numbers.
    """
    data, drug_cols, feat_cols = get_protein_drug_data(drug_class)
    data = drop_na_from_data(data, drug_name, feat_cols)
    data[drug_name] = data[drug_name].apply(lambda x: np.log10(x))

    return data, feat_cols


def to_numeric_rep(df, feat_cols, rep='mw'):
    df_new = df.copy()
    allowed_rep = ['mw', 'pKa']
    assert rep in allowed_rep, 'kwarg rep must be in {0}'.format(allowed_rep)

    if rep == 'mw':
        numeric_dict = molecular_weights
    if rep == 'pKa':
        numeric_dict = isoelectric_points
    for col in feat_cols:
        df_new[col] = df_new[col].replace(numeric_dict.keys(),
                                          numeric_dict.values())

    return df_new


def test_data_integrity(data):
    """
    Data integrity tests! Make sure to use this to check that the data are
    suitable for input into scikit-learn.
    """
    assert '-' not in data
    assert '.' not in data
    assert 'X' not in data
    assert not data.dropna().isnull().values.any()


def to_train_test_split(data, feat_cols, drug_name, test_size=0.1):
    X = data[feat_cols]
    Y = data[drug_name]

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y,
                                                        test_size=test_size)

    return X, Y, X_train, X_test, Y_train, Y_test
