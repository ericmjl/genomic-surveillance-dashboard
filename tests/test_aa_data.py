from gsdash.isoelectric_point import isoelectric_points as ip
from gsdash.molecular_weight import molecular_weights as mw

def test_mixed_amino_acids():
    """
    Tests that the mixed amino acids B and Z are the average
    of D/N and E/Q respectively.
    """
    assert ip['B'] == (ip['D'] + ip['N']) / 2
    assert mw['B'] == (mw['D'] + mw['N']) / 2
    assert ip['Z'] == (ip['E'] + ip['Q']) / 2
    assert mw['Z'] == (mw['E'] + mw['Q']) / 2
