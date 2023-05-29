#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue, May 30th 2023
Written during WIN FMRIB/OHBA Analysis groups' PyPhyPred hackathon
#add names   
A separate copy of the licence is also included within the repository.

This file includes helper functions for reading from / writing to disk.
 
"""

###############################################################################

#import python packages that you need here; e.g. numpy scipy etc

###############################################################################

def load_file(filename, fileformat='', **kwargs):
    #this function takes a filename and loads it into a numpy array
    #if the fileformat is not user-specified, determine it automatically based on filename extension
    #throw an error if file format is not supported
    #example formats that could be good to support: .npy, .npz, .txt, .csv, .mat 
    return #output_array

def save_file(filename, fileformat='', **kwargs):
    #this function takes a filename and loads it into a numpy array
    #if the fileformat is not user-specified, determine it automatically based on filename extension
    #throw an error if file format is not supported
    #example formats that could be good to support: .npy, .npz, .txt, .csv, .mat 
    return #output_array