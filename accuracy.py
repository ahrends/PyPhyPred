#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue, May 30th 2023
Written during WIN FMRIB/OHBA Analysis groups' PyPhyPred hackathon
#add names 
A separate copy of the licence is also included within the repository.

This file includes helper functions for computing prediction accuracies for regression and classification:
https://scikit-learn.org/stable/modules/model_evaluation.html

Names of functions and input/outputs are just illustrative, please change as you see fit  
"""

###############################################################################

#import python packages that you need here; e.g. numpy scipy etc

###############################################################################


def accuracy_regression(y_actual , y_pred, accu_type='corr', other_parameters=None, **kwargs):
    """
    this should include options for computing some of the common accuracy metrics that are suitable 
    for regression-based predictions, e.g. correlation, R2 score, F1 score, RMSE
    see: https://scikit-learn.org/stable/modules/model_evaluation.html#regression-metrics 
    
    
    Parameters
    ----------
    change/add/remove parameters as required
    y_actual, y_pred : ndarray, vectors
    actual and predicted values of target variable
    - if accu_type not specified by the user, use the type of 'normalised accuracy' that's most appropriate for
    data type, here are some data types to consider: continuous, binary, multinomial, ordinal, count, cox
    - for binary, multinomial, ordinal variables, account for when different classes are not balanced 

    Returns
    -------
    returns accuracy value (float), and accuracy type (str). 
    """
    return #add output here

def accuracy_classification(y_actual , y_pred, accu_type='corr', other_parameters=None, **kwargs):
    """
    same as accuracy_regression but for classification: 
    https://scikit-learn.org/stable/modules/model_evaluation.html#classification-metrics 
    
    
    Parameters
    ----------
    change/add/remove parameters as required
    y_actual, y_pred : ndarray, vectors
    actual and predicted values of target variable
    - if accu_type not specified by the user, use the type of 'normalised accuracy' that's most appropriate for
    data type, here are some data types to consider: continuous, binary, multinomial, ordinal, count, cox 

    Returns
    -------
    returns accuracy value (float), and accuracy type (str). 
    """
    return #add output here

def accuracy_plots(fig, ax, y_actual , y_pred, fitted_model, figname, plot_type='roc', other_parameters=None, **kwargs):
    """
    include options for a few common plot types that can be used to evaluate performance, 
    e.g. error distribution for regression, roc for classification etc. 
    
    Parameters
    ----------
    change/add/remove parameters as required
    y_actual, y_pred : ndarray, vectors
    fitted_model : sklearn instance of a fitted model that we'd want to evaluate
    
    Returns
    -------
    returns figure and saves it at figname on disk 
    """
    return #add output here
