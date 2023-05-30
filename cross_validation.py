#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue, May 30th 2023
Written during WIN FMRIB/OHBA Analysis groups' PyPhyPred hackathon
#add names   
A separate copy of the licence is also included within the repository.

This file includes helper functions for cross validation, cross-validated feature preselection, and function for chance-level predictions.
Names of functions and input/outputs are just illustrative, please change as you see fit  
"""

###############################################################################

#import python packages that you need here; e.g. numpy scipy etc

###############################################################################


def cross_validate(x, cv_method='kfold', splits=5, other_parameters=None, **kwargs):
    """
    this should take data and yield train and test indices.
    cv_method could be kfold, stratified kfold, loo, shuffle_split, timeseries_split
    
    
    Parameters
    ----------
    change/add/remove parameters as required
    x : ndarray, Nsubj x k (k does not matter)
         
    
    Returns
    -------
    train_inds, test_inds: ndarray, vectors of indices
    """
    return #add output here


def feature_preselect(x,y,train_inds, nfeat_select=50, target_type='continuous', criteria='corr', other_parameters=None, **kwargs):
    """
    find nfeat_select of best features in feature matrix X, based on similarity of X columns to y within the training set
    
    Parameters
    ----------
    change/add/remove parameters as required
    x - feature matrix: ndarray, Nsubj x Nfeat
    y - target vector: ndarray, Nsubj x 1
    train_inds: ndarray vector
         
    Returns
    -------
    preselected_inds: ndarray, vector of indices for top nfeat_select features in X
    """
    return #add output here


def chance_level_prediction(x,y,train_inds, other_parameters=None, **kwargs):
    """
    create chance-level prediction accuracy for data. This can be a single number or a null distribution
    see here: https://scikit-learn.org/stable/auto_examples/model_selection/plot_permutation_tests_for_classification.html#sphx-glr-auto-examples-model-selection-plot-permutation-tests-for-classification-py

    Parameters
    ----------
    change/add/remove parameters as required
    x - feature matrix: ndarray, Nsubj x Nfeat
    y - target vector: ndarray, Nsubj x 1
    train_inds: ndarray vector
         
    Returns
    -------
    chance_pred: float or ndarray, vector of null distribution 
    """
    return #add output here