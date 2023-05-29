#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue, May 30th 2023
Written during WIN FMRIB/OHBA Analysis groups' PyPhyPred hackathon
#add names   
A separate copy of the licence is also included within the repository.

This file includes helper functions for nonlinear regression methods for various data distributions.
Names of functions and input/outputs are just illustrative, please change as you see fit  

TBC
"""

###############################################################################

#import python packages that you need here; e.g. numpy scipy etc

###############################################################################


def nonlinear_regression_simple(x,y,train_inds,test_inds, target_type='continuous', regularisation_type='elasticnet', options={}, **kwargs):
    """
    this should handle linear regression "without" nested cross-validation;
    include options for various regularisations:
    unregularised regression, Ridge, Lasso and ElasticNet
    
    Parameters
    ----------
    change/add/remove parameters as required
    x - feature matrix: ndarray, Nsubj x Nfeat
    y - target vector: ndarray, Nsubj x 1
    train_inds and train_inds: ndarray vectors 
         
    
    Returns
    -------
    vector of y_test_predicted, and fitted regression model. because, why not.
    """
    return #add output here

def linear_regression_nested_cv(x,y,train_inds,test_inds, target_type='continuous', regularisation_type='elasticnet', options={}, **kwargs):
    """
    same as linear_regression_simple but with nested cross validation for hyperparameter optimisation
    
    Parameters
    ----------
    change/add/remove parameters as required
    x - feature matrix: ndarray, Nsubj x Nfeat
    y - target vector: ndarray, Nsubj x 1
    train_inds and train_inds: ndarray vectors 
         
    Returns
    -------
    vector of y_test_predicted, and fitted regression model. because, why not.
    """
    return #add output here

