#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue, May 30th 2023
Written during WIN FMRIB/OHBA Analysis groups' PyPhyPred hackathon
#add names   
A separate copy of the licence is also included within the repository.

This file includes helper functions for preprocessing X and y.
names of functions and input/outputs are just illustrative, please change as you see fit  
"""

###############################################################################

#import python packages that you need here; e.g. numpy scipy etc

###############################################################################

def get_target_type(y, other_parameters=None, **kwargs):
    """
    Determines what type a target variable is
    
    Parameters
    ----------
    change/add/remove parameters as required
    y (numpy array / vector): Response variable
    
    Returns
    -------
    returns str: 'continuous', 'binary', 'count', 'multinomial', 'ordinal', add more if suitable 

    """
    return #add output here

def get_regression_type(y_type, other_parameters=None, **kwargs):
    """
    Determines what type of linear regression is most suitable for a data type

    Parameters
    ----------
    change/add/remove parameters as required
    y_type (str): output from get_target_type

    Returns
    -------
    returns str: 'linear', 'logistic', etc., add more as appropriate
    """ 
    return #add output here


def handle_outliers(data, other_parameters=None, **kwargs):
    """
    Outlier detection and handling in y  and X

    Parameters
    ----------
    change/add/remove parameters as required
    data (numpy array): vector or matrix

    Returns
    -------
    returns data_no_outlier (numpy array): vector or matrix
    """ 
    return # add output here

def handle_missing(X,y, other_parameters=None, **kwargs):
    """
    Missing data detection and handling in y  and X

    Parameters
    ----------
    change/add/remove parameters as required
    X (numpy array): feature matrix of size Nsubj x Nfeat
    y (numpy array): target vector of size Nsubj x 1

    Returns
    -------
    X_nomissing (numpy array)
    y_nomissing (numpy array)
    note that both outputs should have the same number of rows 
    """
    return #add outputs here

def standardise_data(x, other_parameters=None, **kwargs):
    """
    include options for a few common standardasation techniques; e.g.:
    demeaning, variance normalisation, zscore, tstat, etc.

    Parameters
    ----------
    change/add/remove parameters as required
    x (numpy array): matrix or vector
    
    Returns
    -------
    x_standard (numpy array)

    """
    return #add outputs here

def gaussianise_nonparametric(x, train_inds, other_parameters=None, **kwargs):
    """
    A nonparametric approach to Gaussianise data (quantile transforms could be good)
    Input data could be feature matrix, confound matrix or target vectors
    Learn transforms on train data and fit transform to both train and test  

    Parameters
    ----------
    change/add/remove parameters as required
    x (numpy array): matrix or vector, x.shape[0]=Nsubj
    train_inds (numpy array): vectors of train inds across Nsubj

    Returns
    -------
    x_gaussianised (numpy array) 
    """
    return #add outputs here

def gaussianise_parametric(x, train_inds, other_parameters=None, **kwargs):
    """
    Same as gaussianise_nonparametric but using a parametric approach: I'm not sure if we need this?
    If we do, which scenarios should default to parametric vs non-parametric approach?
    
    Parameters
    ----------
    change/add/remove parameters as required
    x (numpy array): matrix or vector, x.shape[0]=Nsubj
    train_inds (numpy array): vectors of train inds across Nsubj

    Returns
    -------
    x_gaussianised (numpy array) 
    """
    return #add outputs here

def transforms_other(y, train_inds, other_parameters=None, **kwargs):
    """
    Include suitable transforms for non-continuous data types: e.g. count, multinomial etc.
    This would be mainly useful for target variables, but could be applied to feature matrices and confounds too
    
    Parameters
    ----------
    change/add/remove parameters as required
    y (numpy array): matrix or vector, x.shape[0]=Nsubj
    train_inds (numpy array): vectors of train inds across Nsubj

    Returns
    -------
    y_transformed (numpy array) 
    """
    return #add outputs here

def dim_reduction_linear(x, train_inds=None, other_parameters=None, **kwargs):
    """
    Linear dimensionality reduction techniques: SVD, PCA and ICA
    Learn transforms on train data and fit transform to both train and test
    
    Parameters
    ----------
    change/add/remove parameters as required
    x (numpy array): matrix, typically of size Nsubj x Nfeat (could be feature or confounds matrix)
    train_inds (numpy array, optional): vectors of train inds across Nsubj
    
    If train_inds not specified, learn/apply dim reduction across whole Nsubj

    Returns
    -------
    x_reduced (numpy array) 
    """
    return #add outputs here

def dim_reduction_nonlinear(x, train_inds=None, other_parameters=None, **kwargs):
    """
    one or two nonlinear dimensionality reduction techniques: manifold learning, e.g. diffusion embedding
    Learn transforms on train data and fit transform to both train and test
    
    Parameters
    ----------
    change/add/remove parameters as required
    x (numpy array): matrix, typically of size Nsubj x Nfeat (could be feature or confounds matrix)
    train_inds (numpy array, optional): vectors of train inds across Nsubj
    
    If train_inds not specified, learn/apply dim reduction across whole Nsubj

    Returns
    -------
    x_reduced (numpy array) 
    """
    return #add outputs here