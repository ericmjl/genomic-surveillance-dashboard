from sequence_transformer import to_numeric_rep, standardize_sequence

new_sequence1 = 'PQITLNQRTLVTIPIGGDLKEALLKTGADDTVLEEMNLPGRMKPKMIGGIGGFIKVRQYD'\
    'QILIEICGRKAIGTVLVGPTPVNIIGRNLLTQVGCTLNFP'

new_sequence2 = 'PQITLWQRPLVTIKIGGQLKEALLKTGADDTVLEEMNLPGRWKPKMIRGIGGFIKMRQYD'\
    'KILIEICGHQAIGTVLVGPTPRNIIGDNLLTIIGCTLN'

arr1 = to_numeric_rep(new_sequence1)
arr2 = to_numeric_rep(new_sequence2)

arr1_pka = to_numeric_rep(new_sequence1, 'pKa')
arr2_pKa = to_numeric_rep(new_sequence2, 'pKa')

std1 = standardize_sequence(arr1, 'protease')
std2 = standardize_sequence(arr2, 'protease')


def test_to_numeric_rep():
    assert arr1.size == len(arr1)
    assert arr2.size == len(arr2)


def test_standardize_sequence():
    assert std1.shape == std2.shape
