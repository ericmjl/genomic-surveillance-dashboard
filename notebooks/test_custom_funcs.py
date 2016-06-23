import custom_funcs as cf
import pandas as pd

data_protease = pd.read_csv('../data/hiv-protease-data.csv', index_col='SeqID')
protease_drug_cols = data_protease.columns[0:8]
protease_feat_cols = data_protease.columns[8:]

data_nnrt = pd.read_csv('../data/hiv-nnrt-data.csv', index_col='SeqID')
nnrt_drug_cols = data_nnrt.columns[0:4]
nnrt_feat_cols = data_nnrt.columns[4:]

data_nrt = pd.read_csv('../data/hiv-nrt-data.csv', index_col='SeqID')
nrt_drug_cols = data_nrt.columns[0:6]
nrt_feat_cols = data_nrt.columns[6:]


def test_data_integrity():
    assert data_protease.shape == (1808, 107)
    assert data_nnrt.shape == (1765, 244)
    assert data_nrt.shape == (1498, 246)


def test_replace_ambiguous_letters_with_nan():
    df = cf.replace_ambiguous_letters_with_nan(data_protease,
                                               protease_feat_cols)
    assert df.dropna(subset=protease_feat_cols).shape == (871, 107)


def test_replace_dashes_with_canonical_letters():
    """
    TODO.
    """
    pass


def test_read_consensus():
    protease_consensus = cf.read_consensus(
        '../data/hiv-protease-consensus.fasta')
    assert len(protease_consensus) == 99


def test_drop_na_from_data():
    """
    TODO
    """
    df = cf.drop_na_from_data(data_protease, 'FPV', protease_feat_cols)

    # There should be no NaN values.
    assert not df.dropna().isnull().values.any()
