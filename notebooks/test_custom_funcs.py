import custom_funcs as cf
import pytest

data_prot, prot_drug_cols, prot_feat_cols = cf.read_data('protease')
data_nnrt, nnrt_drug_cols, nnrt_feat_cols = cf.read_data('nnrt')
data_nrt, nrt_drug_cols, nrt_feat_cols = cf.read_data('nrt')

prot_consensus = cf.read_consensus('protease')
nnrt_consensus = cf.read_consensus('nnrt')
nrt_consensus = cf.read_consensus('nrt')


def test_read_data():
    assert data_prot.shape == (1808, 107)
    assert data_nnrt.shape == (1765, 244)
    assert data_nrt.shape == (1498, 246)


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


def test_drop_na_from_data():
    df = cf.drop_na_from_data(data_prot, 'FPV', prot_feat_cols)
    # There should be no NaN values.
    assert not df.dropna().isnull().values.any()


def test_get_cleaned_data():
    df = cf.get_cleaned_data('protease', 'FPV')
    assert df.shape == (726, 100)
