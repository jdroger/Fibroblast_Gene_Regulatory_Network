"""
Scripts for generating cross validation folds and datasets
"""

import random

def makeFolds(data, k):
    """
    Given a dataframe of cpm values, randomly 
    creates k folds and returns nested list of 
    column names.
    """
    # randomize columns
    order = data.columns.tolist()
    random.shuffle(order)
    # split into folds (specified by k)
    folds = []
    fold=0
    while fold < k:
        start = fold*k
        end = (fold*k)+4
        folds.append(order[start:end])
        fold = fold+1
    return folds

def assignFolds(folds):
    """
    Given a nested list of strings from makeFolds, 
    assigns testing sets from each list(fold) and 
    assigns rest to training set for each fold. K
    is determined from length of nested list as the
    number of folds.
    """
    training = []
    testing = []
    k = len(folds)
    fold=0
    while fold < k:
        testing.append(folds.pop(fold))
        training.append([y for x in folds for y in x])
        folds.insert(fold, testing[fold])
        fold = fold+1
    return training, testing

