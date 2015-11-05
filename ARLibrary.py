#!/usr/bin/env python

import sys,os
sys.path.append("/opt/Mantid/bin")
from numpy import *
from string import *

import mantid

    

          
def _parseBTPlist(value):
    """
    Helper function to transform a string into a list of integers
    For example "1,2-4,8-10" will become [1,2,3,4,8,9,10]
    It will deal with lists as well, so range(1,4) will still be [1,2,3]
    """
    runs = []
    #split the commas
    parts = str(value).strip(']').strip('[').split(',')
    #now deal with the hyphens
    for p in parts:
        if len(p) > 0:
            elem = p.split("-")
        if len(elem) == 1:
            runs.append(int(elem[0]))
        if len(elem) == 2:
            startelem = int(elem[0])
            endelem   = int(elem[1])
            if endelem < startelem:
                raise ValueError("The element after the hyphen needs to be greater or equal than the first element")
            elemlist  = range(startelem,endelem+1)
            runs.extend(elemlist)
    return runs          
           



    
    
def Load_2_Monitors(fname):
    """
    Load all the monitors and extract only monitor 1 and monitor 2
    and put them in the workspace ows
    GEG 10.15.2013
    """
    wm=mantid.simpleapi.LoadNexusMonitors(Filename=fname)
    nsp=wm.getNumberHistograms()
    if nsp < 2:
       raise ValueError("There are less than 2 monitors")

    # determine which spectrum number has which spectra id
    # from A. T. Savici
    for sp in range(nsp):
       if wm.getSpectrum(sp).getDetectorIDs()[0]==-wm.getInstrument().getNumberParameter('ei-mon1-spec')[0]:
                    sp1=sp
       if wm.getSpectrum(sp).getDetectorIDs()[0]==-wm.getInstrument().getNumberParameter('ei-mon2-spec')[0]:
                    sp2=sp
    #extract 2 good monitors into individual workspaces
    tsp1=mantid.simpleapi.ExtractSingleSpectrum(wm,sp1)
    tsp2=mantid.simpleapi.ExtractSingleSpectrum(wm,sp2)
    #Ensure the spectrumnumber for detector ID -1 is 1 and
    #for detector ID -2 is 2
    tsp1.getSpectrum(0).setSpectrumNo(int(wm.getInstrument().getNumberParameter('ei-mon1-spec')[0]))
    tsp2.getSpectrum(0).setSpectrumNo(int(wm.getInstrument().getNumberParameter('ei-mon2-spec')[0]))
    #create a workspace that is monitor 1 first and monitor 2 second
    mantid.simpleapi.ConjoinWorkspaces('tsp1','tsp2')
    # clean up workspaces
    ows = mantid.simpleapi.CloneWorkspace(InputWorkspace='tsp1')
    mantid.simpleapi.DeleteWorkspace('tsp1')
    
    return ows
