import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.preprocessing import StandardScaler


class Feature_table:
    # Class constructor input: pandas dataframe in pickle format
    # adds hydrophobicity column and store it as class member

    def __init__(self, pickle):
        self.df = pd.read_pickle(pickle)
        self.add_hydrophobicity()

    def add_hydrophobicity(self):
        gravy_list = self.get_gravy_list()  # get normalize gravy for all seqs
        self.df['Hydrophobicity'] = gravy_list  # add column to df

    def get_gravy_list(self):
        gravy_list = []
        for seq in self.df.index:  # for every seq, add gravy to list
            seq = ProteinAnalysis(seq)
            gravy = "{:.6f}".format(seq.gravy())
            gravy_list.append(gravy)
        gravy_list = np.array(gravy_list)  # convert to np array
        return self.normalize(gravy_list)  # return normalized

    def normalize(self, np_array):  # normalize np array
        np_array = np_array.reshape(-1, 1)
        scaler = StandardScaler()
        scaler.fit(np_array)
        np_array = scaler.transform(np_array)
        return np_array
