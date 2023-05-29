#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue, May 30th 2023
Written during WIN FMRIB/OHBA Analysis groups' PyPhyPred hackathon 
#add names 
A separate copy of the licence is also included within the repository.

This file includes helper functions for confound modelling.
Names of functions and input/outputs are just illustrative, please change as you see fit  
note: see preprocessing.py for functions that you may need here for preprocessing confounds (e.g. demeaning etc.)
Team 1 will be in charge of that. 
"""

###############################################################################

#import python packages that you need here; e.g. numpy scipy etc

###############################################################################


def deconfound(x, y, confmat, train_inds, test_inds, deconfound_x=True, deconfound_y=True, other_parameters=None, **kwargs):
    """
    this function does [linear] deconfounding of feature and target matrices
    
    Parameters
    ----------
    change/add/remove parameters as required
    x : ndarray, feature matrix: Nsubj x Nfeature
    y : ndarray, target vector: Nsubj x 1
    confmat: ndarray, confounds matrix: Nsubj x Nconf
    - train_inds and test_inds should be passed on from the layer above
    - allow user to choose if they want to do one-sided or two-sided deconfounding via deconfound_x and deconfound_y
    - note: we need deconfounding done within cross-validation loop; i.e. betas estimated based on train set and applied to test set    
    
    Returns
    -------
    returns: deconfounded x and y , means of original x and y, and betas  
    """
    return #add output here

def reconfound(y_actual, y_pred,  other_parameters=None, **kwargs):
    """
    when actual_y is deconfounded prior to prediction, the values change 
    (e.g. age will be demeaned, gaussianised, confounds will be regressed out etc.)
    we need this function to bring actual and predicted y into to the original space  
    
    Parameters
    ----------
    change/add/remove parameters as required
    y_actual, y_pred : ndarray, each is a vector 
    
    Returns
    -------
    returns: unconfounded  y_actual, y_pred
    """
    return #add output here

def confounds_predict(y, confmat, train_inds, test_inds, other_parameters=None, **kwargs):
    """
    this function is to give us a baseline prediction that is obtained based on confmat only
    note: do we need this as a separate function?
    is there anything distinct about code to make predictions based on confounds vs based on actual x matrix? 
    
    Parameters
    ----------
    change/add/remove parameters as required
    y: ndarray, target vector: Nsubj x 1
    confmat: ndarray, confounds matrix: Nsubj x Nconf
    - train_inds and test_inds should be passed on from the layer above
    - allow user to choose if they want to do one-sided or two-sided deconfounding via deconfound_x and deconfound_y
    - note: we need deconfounding done within cross-validation loop; i.e. betas estimated based on train set and applied to test set    
    
    Returns
    -------
    returns: y_pred_confounds  
    """
    return #add output here