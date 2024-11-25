from galaxy_predictor import CodonCounter
from galaxy_predictor import GalaxyPredictor


def test_simple_prediction():
    cc = CodonCounter()
    cc.count_file("tests/Nanobdella_aerobiophila.fasta")
    table = cc.get_table()
    predictor = GalaxyPredictor("models/archaea.pkl")
    predicted_class = predictor.predict(table)[0]
    assert predicted_class == 1
    assert predictor.get_label(predicted_class) == "CENTER"


def test_invalid_model():
    cc = CodonCounter()
    table = cc.get_table()

    not_found = False
    try:
        predictor = GalaxyPredictor("models/invalid.pkl")
        predictor.predict(table)
    except FileNotFoundError as e:
        not_found = True
    assert not_found
