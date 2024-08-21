#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score

def resistancegeneLogisticRegressor(fasta_file, \
                                       expression_file, \
                                             prediction_motif, \
                                             prediction_size, \
                                                      output_file):
    """
    applying a logistic regressor specific for the training on the plant
    resistance genes. This applies the logistic regession based
    on the sequence characteristics of the plant disease resistance
    genes and the corresponding expression profile. A normalized 
    log transformed expression methods can be used.
    @fasta_file: file containing the plant disease resistance genes
    @expression_file: file containing the expression profile for
    those genes. it also takes a prediction motif profile which
    defines the presence and the absence of the resistance genes
    and makes a probability index. The function returns a accuracy score
    and writed the model to an output file. You can put a meshgrid to 
    enable the visualization of the model using the np.array
    You can change the save to a pickle file for the model
    """
    sequence_file_train_read = list(filter(None,[x.strip() for \
                                                 x in open(fasta_file).readlines()]))
    sequence_train_dict = {}
    for i in sequence_file_train_read:
        if i.startswith(">"):
            path = i.strip()
            if i not in sequence_train_dict:
                sequence_train_dict[i] = ""
                continue
        sequence_train_dict[path] += i.strip()
    ids = list(sequence_train_dict.keys())
    sequences = list(sequence_train_dict.values())
    sequence_dataframe = pd.DataFrame([(i,j)for i,j in zip(ids, sequences)]). \
                                          rename(columns = {0: "ids", 1: "sequence"})
    sequence_dataframe["genomic_content"] = sequence_dataframe["sequence"]. \
                                  apply(lambda n: round((n.count("G")+n.count("C")),2)/len(n))
    sequence_expression = []
    with open(expression_file, "r") as expression:
        for line in expression.readlines():
            sequence_expression.append(line.strip().split())
    sequence_dataframe["expression"] = sequence_expression
    sequence_motif = prediction_motif
    sequence_dataframe["classify"] = sequence_dataframe["sequence"]. \
                                     apply(lambda n: "0" if sequence_motif in n else "1")
    sequence_dataframe_2d_matrix = [(i,j) for i,j in zip(sequence_dataframe["genomic_content"] \
                                                .to_list(), sequence_dataframe["expression"].to_list())]
    sequence_scoring_matrix = [i for i in sequence_dataframe["classify"].to_list()]
    test_size = prediction_size
    X_train, X_test, y_train, y_test = train_test_split(sequence_dataframe_2d_matrix, 
                                                        sequence_scoring_matrix,
                                                        test_size = 0.25, 
                                                        random_state = 0)
    scaler_X = StandardScaler()
    X_train = scaler_X.fit_transform(X_train)
    X_test = scaler_X.transform(X_test)
    model = LogisticRegression(solver = 'liblinear', \
                                    random_state = 0).fit(X_train, y_train)
    y_pred = model.predict(X_test)
    model.predict_proba(sequence_dataframe_2d_matrix)[:,1]
    accuracy = accuracy_score(y_test, y_pred)
    with open("outfile", "w") as model_file:
        model_file.write(print(model))
    return print(f"{accuracy*100}")
