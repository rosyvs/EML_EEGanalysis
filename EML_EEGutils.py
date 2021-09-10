
import glob
import os
import mne
import meegkit #nbara.github.io/python-meegkit/
import numpy as np
import scipy.io as sio # loads/writes matlab files anod other useful io 
import matplotlib.pyplot as plt
import time
import pyprep_local
from pyprep_local.pyprep.find_noisy_channels import NoisyChannels
from pyprep_local.pyprep.utilities import _union, _set_diff

# EEG preprocessing steps for EML
class emleeg(object):
    def __init__(self, data_dir="./EEG/"):
        self.data_dir = data_dir
        self.sampleRate = np.nan
        self.montage = None

    def loadFile(self, filename):
        self.raw = mne.io.read_raw_brainvision(self.data_dir + filename, preload=True)
        ch_names = self.raw.ch_names
        self.prepro = self.raw.copy() # make a copy of the data which will be modified in preprocessing
        self.prepro.info['description'] = self.data_dir + filename # store filename of orig data in mne obj
        print('Data channels: ' + ' '.join(ch_names))

    def extractInfo(self):
        # Extract some info
        eeg_index = mne.pick_types(self.raw.info, eeg=True, eog=False, meg=False)
        ch_names = self.raw.info["ch_names"]
        self.ch_names_eeg = list(np.asarray(ch_names)[eeg_index])
        self.sampleRate = self.raw.info["sfreq"]
        print("Sample Rate : {0}".format(self.sampleRate))
        print("Number of Samples : {0}".format(self.prepro.n_times))

    def robustDetrend(self,order=3,weights=None):
        # uses meegkit from nbara
        # for some reason I can only get this to work on each channel individually
        print("detrending data using meegkit robust detrending. Go make a cuppa :-)")

        # extract data matrix and info from mne object
        X = self.prepro.get_data().T # get eeg channels
        eegchannels = mne.pick_types(self.raw.info,eeg=True)
        Z=X # this will be the detrended array to put back into an mne object
        # loop over channels
        for c in eegchannels:
            x = X[:,c]
            if weights is None:
                weights=np.ones_like(x)
            z, _, _ = meegkit.detrend.detrend(x, order, weights,basis='polynomials', threshold=5, n_iter=4, show=False)
            Z[:,c] = z
            print("...completed detrending for channel {0} of {1}".format(c+1,len(eegchannels)))

        self.prepro._data = Z.T

    def zapline(self, nremove=2):
        print("Zapping 60Hz line noise with DSS")
        X = self.prepro.get_data().T
        eegchannels = mne.pick_types(self.prepro.info,eeg=True)
        Z=X # this will be the line-zapped array to put back into an mne object
        # print(Z.shape)
        # print(Z[:,eegchannels].shape)
        z, artifact = meegkit.dss.dss_line(X[:,eegchannels],fline=60.0, sfreq=self.sampleRate, nremove=nremove)
        Z[:,eegchannels] = z
        self.prepro._data = Z.T
        return self, artifact 

    def checkChannels(self):
        # use noisy channel detection from PREP pipeline
        # We use the local copy of pyprep from APril 2020 but pyprep has since been updated
        # TODO what has changed? I see the dev has removed a lot of the 'noisy' functions which is exactly what we use here - are there issues? 
        pyprepobj = NoisyChannels(self.prepro)
        print("Detecting bad channels by NaN, flat, deviation, HF noise")
        pyprepobj.find_bad_by_nan_flat()
        pyprepobj.find_bad_by_deviation()
        pyprepobj.find_bad_by_hf_noise()
        #pyprepobj.find_bad_by_correlation()
        bads = pyprepobj.get_bads(verbose = True)
        print("Bad channels: {0} of {1}".format(len(bads), len(mne.pick_types(self.raw.info,eeg=True))))
        self.prepro.info['bads']= bads
        return bads

def chn_name_mapping(ch_name):
    """Map channel names to fit standard naming convention."""
    ch_name = ch_name.strip('.')
    ch_name = ch_name.upper()
    if 'Z' in ch_name:
        ch_name = ch_name.replace('Z', 'z')

    if 'FP' in ch_name:
        ch_name = ch_name.replace('FP', 'Fp')

    return ch_name

def print_and_log(str,logfile='EMLEEG_pylog.txt'):
    with open(logfile, 'a') as o:
        o.write(str)
        o.write("\n")
    print(str)

