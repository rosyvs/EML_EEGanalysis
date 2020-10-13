
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
        print("Reading {0}".format(filename))
        self.raw = mne.io.read_raw_brainvision(self.data_dir + filename, preload=True)
        #self.raw = mne.io.read_raw_edf(self.data_dir + filename, preload=True)
        #self.raw.rename_channels(chn_name_mapping)
        ch_names = self.raw.ch_names
        # self.raw.plot_psd()
        # self.raw.plot(duration=5)
        print(ch_names)

    def extractInfo(self):
        # Extract some info
        eeg_index = mne.pick_types(self.raw.info, eeg=True, eog=False, meg=False)
        ch_names = self.raw.info["ch_names"]
        self.ch_names_eeg = list(np.asarray(ch_names)[eeg_index])
        sample_rate = self.raw.info["sfreq"]
        self.sampleRate = self.raw.info["sfreq"]

        # Make a copy of the data
        self.raw_copy = self.raw.copy()

    def robustDetrend(self,order=3,weights=None):
        # uses meegkit from nbara
        # for some reason I can only get this to work on each channel individually
        print("detrending data using meegkit robust detrending. Go make a cuppa :-)")

        # extract data matrix and info from mne object
        X = self.raw.get_data().T # get eeg channels
        eegchannels = mne.pick_types(self.raw.info,eeg=True)
        Z=X # this will be the detrended array to put back into an mne object
        # loop over channels
        for c in eegchannels:
            x = X[:,c]
            if weights is None:
                weights=np.ones_like(x)
            z, _, _ = meegkit.detrend.detrend(x, order, weights,basis='polynomials', threshold=5, n_iter=4, show=False)
            Z[:,c] = z
            print("...completed detrending for channel {0} of {1}".format(c,len(eegchannels)))

        #self.raw = mne.io.RawArray(Z.T,self.raw.info)
        self.raw._data = Z.T

    def checkChannels(self):
        # use noisy channel detection from PREP pipeline
        # We use the local copy of pyprep from APril 2020 but pyprep has since been updated
        # TODO what has changed? I see the dev has removed a lot of the 'noisy' functions which is exactly what we use here - are there issues? 
        pyprepobj = NoisyChannels(self.raw)
        pyprepobj.find_bad_by_nan_flat()
        pyprepobj.find_bad_by_deviation()
        pyprepobj.find_bad_by_hf_noise()
        #pyprepobj.find_bad_by_correlation()
        bads = pyprepobj.get_bads(verbose = True)

        pyprepobj.printAndWrite("Sample Rate : {0}".format(self.sampleRate))
        pyprepobj.printAndWrite("Number of Samples : {0}".format(self.raw.n_times))
        pyprepobj.printAndWrite("Channels Recorded : {0}".format(len(self.raw.info["ch_names"])))
        pyprepobj.printAndWrite("Number of EEG channels Recorded {0}".format(len(self.ch_names_eeg)))
        pyprepobj.printAndWrite("Number of Bad channels: {0}".format(len(bads)))
        self.raw.info['bads']= bads
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
