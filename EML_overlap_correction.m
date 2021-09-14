% EML overlap correction

clear all; close all
init_unfold

% use only file w reliable trigger
hasTriggerList =readtable('triggerSources.csv');
sublist = find(hasTriggerList.sdcard==1);
sublist = 60:70
exclude = [20, 21,22, 26,73, 77]; % Subj to exclude because no eeg or no trigger etc
sublist = sublist(~ismember(sublist,exclude));

dir_raw = '/Volumes/Blue1TB/EyeMindLink/Data';
dir_pre = fullfile('..','..','Data','EEG_processed') ;
mkdir( dir_pre, 'overlap_corrected')