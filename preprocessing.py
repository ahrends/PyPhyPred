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

def get_target_type(y):
    """
    Determines what type a target variable is
    Parameters:
    y (numpy array / vector): Response variable
    Returns:
    str: 'continuous', 'binary', 'count', 'multinomial', 'ordinal', add more if suitable
    """
    return #y_type add output string here

def get_regression_type(y_type):
    """
    Determines what type of linear regression is most suitable for a data type
    Parameters:
    y_type (str): output from get_target_type
    Returns:
    str: 'linear', 'logistic', etc., add more as appropriate
    """ 
    return #reg_type add output string here


def handle_outliers(data):
    """
    Outlier detection and handling in y  and X
    Parameters:
    data (numpy array): vector or matrix
    Returns:
    data_no_outlier (numpy array): vector or matrix
    """ 
    return #reg_type add output string here

def handle_missing(X,y):
    """
    Missing data detection and handling in y  and X
    Parameters:
    X (numpy array): feature matrix of size Nsubj x Nfeat
    y (numpy array): target vector of size Nsubj x 1
    Returns:
    X_nomissing (numpy array)
    y_nomissing (numpy array)
    note that both outputs should have the same number of rows 
    """
    return #X_nomissing, y_nomissing add output string here

def gaussianise_nonparametric(X,nquantiles=None,X_trans_qt=None):
    #Xmat should be nsample x nfeature: it will be gaussianised "per column"
    Xmat_tmp=Xmat.copy()
    if nquantiles is None:
        nquantiles=int(Xmat_tmp.shape[0]/5)
    if X_trans_qt is None:
        qt = QuantileTransformer(n_quantiles=nquantiles, output_distribution='normal',random_state=0)    
        X_trans_qt = qt.fit(Xmat_tmp)
    Xmat_normal=X_trans_qt.transform(Xmat_tmp)
    Xmat_normal = np.nan_to_num(Xmat_normal)
    return Xmat_normal,X_trans_q

def gaussianise_nonparametric(X,nquantiles=None, other_parameters=None, **kwargs):
    """
    this should handle various boxplot-like vis e.g. violinplot, raincloud plots, violinplots etc
    
    Parameters
    ----------
    change/add/remove parameters as required
    prediction_accuracy : ndarray, or list of vectors, or dictionary , or panda dataframe
        Every vector includes accuracy obtained from a different method.
        Would be nice to allow comparisons of groups of methods: 
        e.g. a1/a2/a3 from method a, b1/b2/b3 from method b etc., with boxplots of each group stuck together 
    
    Returns
    -------
    returns figure and saves it at figname on disk 
    """
    return #add output here