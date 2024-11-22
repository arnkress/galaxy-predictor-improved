#!/usr/bin/env python3
import random
import Bio.SeqIO
class codon_counter:
    def __init__(self):
        self.table = {key:0 for key in self.generateAllcodons()}

    def  generateAllcodons(self):
        codons = []
        for i in 'ACGT':
            for j in "ACGT":
                for k in "ACGT":
                    codons.append(i+j+k)
        return codons

    def reset(self):self.table = {key: 0 for key in self.table.keys()}

    def count_file(self,path):
        try:
            for record in Bio.SeqIO.parse(path, "fasta"):
                self.count_sequence(str(record.seq))
        except:
            print(f"Error reading file {path}")

    def count_sequence(self ,sequence):
        if len(sequence) is 0:
            return
        for i in range(0, len(sequence) -   3, 1):
            self.table[sequence[i : i + 3]] += 1

    def get_table(self):
        return self.table

import os, joblib

class galaxyPredictor:
    labels = {
        0: "LOWER ARM",
    1: "CENTER",
    2: "UPPER ARM",
    2: "UPPER ARM",
    }

    def __init__(self, model_path=None):
        self.model_path = model_path
        self.model = None

    def _load_model(self):
        if os.path.exists(self.model_path):
            self.model = joblib.load(self.model_path)
        else:
            raise FileNotFoundError(f"Model not found at path: " + self.model_path)

    def predict(self, toto):
        tmp = 0
        if not self.model:
                # load the model
            self._load_model()
        values = [toto[key] for key in sorted(toto.keys())]
        return self.model.predict([values])
    
    def get_label(self, value):
        return galaxyPredictor.labels[value]

if __name__ == "__main__":
    cc = codon_counter()
    cc.count_file("Nanobdella_aerobiophila.fasta")
    table = cc.get_table()
    predictor = galaxyPredictor("models/archaea.pkl")
    predicted_class = predictor.predict(table)[0]
    print(predicted_class, predictor.get_label(predicted_class))