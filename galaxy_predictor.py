#!/usr/bin/env python3
import Bio.SeqIO
import joblib
import os


class CodonCounter:
    """
    Object to count codons in FASTA files
    """

    def __init__(self):
        self.table = {key: 0 for key in self._generate_all_codons()}

    def _generate_all_codons(self):
        """
        Generate all possible codons

        Returns:
            list: a list of all possible codons
        """
        codons = []
        for i in "ACGT":
            for j in "ACGT":
                for k in "ACGT":
                    codons.append(i + j + k)
        return codons

    def reset(self):
        """
        Reset the count table
        """
        self.table = {key: 0 for key in self.table.keys()}

    def count_file(self, path):
        """
        Count codons in a FASTA file

        Args:
            path (str): the path to the FASTA file

        Raises:
            FileNotFoundError: if the file is not found
        """
        try:
            for record in Bio.SeqIO.parse(path, "fasta"):
                self.count_sequence(str(record.seq))
        except FileNotFoundError:
            print(f"Error reading file {path}")

    def count_sequence(self, sequence):
        """
        Count codons in a sequence

        Args:
            sequence (str): the sequence to count
        """
        if len(sequence) == 0:
            return
        for i in range(0, len(sequence) - 2, 1):
            self.table[sequence[i : i + 3]] += 1

    def get_table(self):
        return self.table


class GalaxyPredictor:
    labels = {
        0: "LOWER ARM",
        1: "CENTER",
        2: "UPPER ARM",
    }

    def __init__(self, model_path=None):
        self.model_path = model_path
        self.model = None

    def _load_model(self):
        if os.path.exists(self.model_path):
            self.model = joblib.load(self.model_path)
        else:
            raise FileNotFoundError("Model not found at path: " + self.model_path)

    def predict(self, toto):
        """
        Predicts the galaxy region for a given input.

        Returns:
            int: the predicted class

        Raises:
            FileNotFoundError: if the model is not found
        """
        if not self.model:
            # load the model
            self._load_model()
        values = [toto[key] for key in sorted(toto.keys())]
        return self.model.predict([values])

    def get_label(self, value):
        """
        Get the label for a given class

        Args:
            value (int): the class value

        Returns:
            str: the class label
        """
        return self.labels[value]


if __name__ == "__main__":
    cc = CodonCounter()
    cc.count_file("tests/Nanobdella_aerobiophila.fasta")
    table = cc.get_table()
    predictor = GalaxyPredictor("models/archaea.pkl")
    predicted_class = predictor.predict(table)[0]
    print(predicted_class, predictor.get_label(predicted_class))
