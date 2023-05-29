#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue, May 30th 2023
Written during WIN FMRIB/OHBA Analysis groups' PyPhyPred hackathon
#add names  
A separate copy of the licence is also included within the repository.

This file includes helper functions for visualisation of various result summaries.
Names of functions and input/outputs are just illustrative, please change as you see fit  
"""

###############################################################################

#import python packages that you need here; e.g. numpy scipy etc

###############################################################################


def boxplot(fig, ax, prediction_accuracy, figname, other_parameters=None, **kwargs):
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

def density_scatter(fig, ax, y_actual , y_pred, other_parameters=None, **kwargs):
    """
    this should handle scatter plots and density scatter plots for comparing y_actual and y_pred
    
    Parameters
    ----------
    change/add/remove parameters as required
    y_actual, y_pred : ndarray, each is a vector 
    
    Returns
    -------
    returns figure and saves it at figname on disk 
    """
    return #add output here

def bland_altman(fig, ax, accuracy_a , accuracy_b, other_parameters=None, **kwargs):
    """
    this is for comparing prediction accuraries obtained from two methods, a and b 
    
    Parameters
    ----------
    change/add/remove parameters as required
    accuracy_a , accuracy_b : ndarray, each is a vector 
    
    Returns
    -------
    returns figure and saves it at figname on disk 
    """
    return #add output here


def plot_image(fig, ax, accuracy_mat, other_parameters=None, **kwargs):
    """
    receives a matrix of accuracies, and visualises it as a matrix
    
    Parameters
    ----------
    change/add/remove parameters as required
    accuracy_mat : ndarray 
    
    Returns
    -------
    returns figure and saves it at figname on disk 
    """
    return #add output here

def plot_manhatan(fig, ax, accuracy_mat, other_parameters=None, **kwargs):
    """
    like big40 but for phenotype predictions. Ideally with sub-plots when you click on the dot, like the genetic ones
    
    Parameters
    ----------
    change/add/remove parameters as required
    accuracy_mat : ndarray, or dict, or dataframe. each vector in the matrix would be prediction accuracies for a single phenotype
    titles can be received as separate entries or embedded within dict/dataframe
    
    Returns
    -------
    returns figure and saves it at figname on disk 
    """
    return #add output here