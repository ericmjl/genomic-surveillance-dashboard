import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


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
    """
    assert isinstance(df, pd.DataFrame)
    assert isinstance(consensus_map, dict)

    """
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
    # Impute consensus sequence
    for i, col in enumerate(feat_cols):
        # Replace '-' with the consensus letter.
        df[col] = df[col].replace({'-': consensus_map[i+1]})

    df = df.replace({'X': np.nan})
    df = df.replace({'.': np.nan})

    return df


def read_consensus(handle):
    """
    Reads in the consensus sequence, makes a map of position to letter.
    """
    consensus = SeqIO.read(handle, 'fasta')
    consensus_map = {i: letter for i, letter in enumerate(str(consensus.seq))}

    return consensus_map


def drop_na_from_data(df, drug_name, feat_cols):

    # Defensive programming checks
    assert isinstance(df, pd.DataFrame)
    assert isinstance(drug_name, str)
    assert drug_name in df.columns, "{0} not in data.".format(drug_name)

    columns = list(feat_cols)
    columns.append(drug_name)

    return df[columns].dropna()
