from galaxy_predictor import CodonCounter


def test_initialize():
    cc = CodonCounter()
    assert cc.get_table() == {
        "AAA": 0,
        "AAC": 0,
        "AAG": 0,
        "AAT": 0,
        "ACA": 0,
        "ACC": 0,
        "ACG": 0,
        "ACT": 0,
        "AGA": 0,
        "AGC": 0,
        "AGG": 0,
        "AGT": 0,
        "ATA": 0,
        "ATC": 0,
        "ATG": 0,
        "ATT": 0,
        "CAA": 0,
        "CAC": 0,
        "CAG": 0,
        "CAT": 0,
        "CCA": 0,
        "CCC": 0,
        "CCG": 0,
        "CCT": 0,
        "CGA": 0,
        "CGC": 0,
        "CGG": 0,
        "CGT": 0,
        "CTA": 0,
        "CTC": 0,
        "CTG": 0,
        "CTT": 0,
        "GAA": 0,
        "GAC": 0,
        "GAG": 0,
        "GAT": 0,
        "GCA": 0,
        "GCC": 0,
        "GCG": 0,
        "GCT": 0,
        "GGA": 0,
        "GGC": 0,
        "GGG": 0,
        "GGT": 0,
        "GTA": 0,
        "GTC": 0,
        "GTG": 0,
        "GTT": 0,
        "TAA": 0,
        "TAC": 0,
        "TAG": 0,
        "TAT": 0,
        "TCA": 0,
        "TCC": 0,
        "TCG": 0,
        "TCT": 0,
        "TGA": 0,
        "TGC": 0,
        "TGG": 0,
        "TGT": 0,
        "TTA": 0,
        "TTC": 0,
        "TTG": 0,
        "TTT": 0,
    }


def test_count_sequence():
    cc = CodonCounter()
    cc.count_sequence("CGTCACA")
    table = cc.get_table()
    assert (
        table["CGT"] == 1
        and table["GTC"] == 1
        and table["TCA"] == 1
        and table["CAC"] == 1
        and table["ACA"] == 1
    )


def test_count_file():
    cc = CodonCounter()
    cc.count_file("tests/Nanobdella_aerobiophila.fasta")
    table = cc.get_table()
    assert (
        table["AAA"] == 35418
        and table["AAC"] == 6871
        and table["AAG"] == 12172
        and table["AAT"] == 34964
    )


def test_reset():
    cc = CodonCounter()
    cc.count_sequence("ACGTACCTG")
    cc.reset()
    assert sum(cc.get_table().values()) == 0


def test_empty_sequence():
    cc = CodonCounter()
    cc.count_sequence("")
    assert sum(cc.get_table().values()) == 0


def test_invalid_file():
    cc = CodonCounter()
    not_found = False
    try:
        cc.count_file("tests/invalid.fasta")
    except FileNotFoundError as e:
        print("File not found")
        not_found = True
    assert not_found
