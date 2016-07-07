import custom_funcs as cf
import pytest

data_prot, prot_drug_cols, prot_feat_cols = cf.read_data('protease')
data_prot_full, _, _ = cf.read_data('protease', sparse=False)
data_nnrt, nnrt_drug_cols, nnrt_feat_cols = cf.read_data('nnrt')
data_nnrt_full, _, _ = cf.read_data('nnrt', sparse=False)
data_nrt, nrt_drug_cols, nrt_feat_cols = cf.read_data('nrt')
data_nrt_full, _, _ = cf.read_data('nrt', sparse=False)

prot_consensus = cf.read_consensus('protease')
nnrt_consensus = cf.read_consensus('nnrt')
nrt_consensus = cf.read_consensus('nrt')


def test_read_data():
    assert data_prot.shape == data_prot_full.shape
    assert data_nnrt.shape == data_nnrt_full.shape
    assert data_nrt.shape == data_nrt.shape


def test_replace_ambiguous_letters_with_nan():
    df = cf.replace_ambiguous_letters_with_nan(data_prot,
                                               prot_feat_cols)
    assert df.dropna(subset=prot_feat_cols).shape == (871, 107)


def test_replace_dashes_with_canonical_letters():
    df = cf.replace_dashes_with_canonical_letters(data_prot,
                                                  prot_consensus,
                                                  prot_feat_cols)

    assert '-' not in df
    assert '.' not in df
    assert 'X' not in df


def test_read_consensus():
    assert len(prot_consensus) == 99
    assert len(nnrt_consensus) == len(nrt_consensus)

    with pytest.raises(AssertionError):
        cf.read_consensus('bogus_name')
    with pytest.raises(AssertionError):
        cf.read_consensus('rt')


def test_to_numeric_rep():
    numeric_prot = cf.to_numeric_rep(data_prot, prot_feat_cols)
    assert not numeric_prot[prot_feat_cols].isnull().values.any()
    numeric_nnrt = cf.to_numeric_rep(data_nnrt, nnrt_feat_cols, 'pKa')
    assert not numeric_nnrt[nnrt_feat_cols].isnull().values.any()


def test_drop_na_from_data():
    df = cf.drop_na_from_data(data_prot, 'FPV', prot_feat_cols)
    # There should be no NaN values.
    assert not df.dropna().isnull().values.any()


def test_to_train_test_split():
    splitted_data = cf.to_train_test_split(data_prot_full, prot_feat_cols,
                                           'ATV')

    X, Y, X_train, X_test, Y_train, Y_test = splitted_data

    assert len(X) == len(Y)
    assert len(X_train) == len(Y_train)
    assert len(X_test) == len(Y_test)
    assert len(X_train) + len(X_test) == len(X)


def test_get_cleaned_data():
    df, feat_cols = cf.get_cleaned_data('protease', 'FPV')
    assert df.shape == (726, 100)
    assert len(feat_cols) == 99
