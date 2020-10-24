from EML_EEGutils import *
import pandas as pd
import re

# EEG analysis for EML: preprocessing
reprepro = False # re-do preprocessing. If False, checks outdir and skips participants with a prepro'd file
datadir = os.path.normpath('C:/Users/roso8920/Dropbox (Emotive Computing)/EyeMindLink/Data')
outdir = os.path.normpath('../../Data/EEG_processed/')# save directory for processed EEG data
fnroot = 'EML1_'
participants = range(27,63) # recall that range is exclusive in Python 
# TODO automatic generation of file list from particpants & loop

for p in participants:
    pID= fnroot + '{:03d}'.format(p)  
    print("\n\n\n~~~~~pID: {0}~~~~~".format(pID))
    prepro_exists = False
    # check for existing preprocessed data first
    if os.path.isfile(os.path.normpath(outdir + '/' + pID + '_p.fif')):
        print("\n\n\n{0}: Existing preprocessed file found in {1}".format(pID,os.path.normpath(outdir + '/' + pID + '_p.fif')))
        prepro_exists = True 
        if reprepro:
            print("{0}: Will overwrite existing preprocessed file ...".format(pID))
   
    if reprepro or not prepro_exists: 
        print("\n\n\n{0}: Preprocessing pipeline started...".format(pID))
        pIDpath = os.path.normpath(datadir + "/" + pID + "/EEG" )
        print("Searching folder {0}".format(pIDpath))
        vhdrList = glob.glob(pIDpath + "/LA" +  "*.vhdr") # filenames starting 'LA' for SD card recording
        if not(vhdrList):
            print('No EEG header file found!')
        else:
            print('Will process the following EEG file(s): {0}'.format(vhdrList))

        for f in vhdrList:
            # preprocessing
            e = emleeg(data_dir='')
            e.loadFile(f)
            e.extractInfo()
            
            e.raw_copy = e.raw.copy()
            # Make a copy of the data before preprocessing

            e.robustDetrend(order=11) # uses meegkit
            e.checkChannels() # uses pyprep
            e.raw.set_montage('standard_1005')
            e.raw.interpolate_bads(reset_bads=False) # interpolate bad channels based on neighbours but keep them marked as 'bad'
            e.raw.set_eeg_reference() # ref to average
            e.raw.notch_filter(np.arange(60, 241, 60)) # powerline filter

            # write preprocessed data to file
            e.raw.save(fname=os.path.normpath(outdir + '/' + pID + '_p.fif'))
    else:
        print("{0}: Skipping preprocessing pipeline.".format(pID))
