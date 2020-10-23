from EML_EEGutils import *
import pandas as pd
import re

# EEG analysis for EML: find localiser-evoked potentials 
reprepro = False # re-do preprocessing. If False, attempts to load existing preprocessed data from outdir. 
datadir = os.path.normpath('C:/Users/roso8920/Dropbox (Emotive Computing)/EyeMindLink/Data')
outdir = os.path.normpath('../../Data/EEG_processed/')# save directory for processed EEG data
fnroot = 'EML1_'
participants = range(50,51) # recall that range is exclusive in Python 
# TODO automatic generation of file list from particpants & loop
p=0
pID= fnroot + '{:03d}'.format(participants[p])  
print("\n\n\n~~~~~Processing pID~: {0}".format(pID))
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

    # get events
    log_file = os.path.normpath(datadir + "/" + pID + "/" + pID + '_Trials.txt')
    log_cols = ['date','time','mode','msg1','msg2','event','value']
    trials = pd.read_csv(log_file, sep='\t',skiprows=1, header=None,names=['col1','event','value'])
    trials[['date','time','mode','msg1','msg2']]=trials.col1.str.split(expand=True)
    trials=trials.drop(['col1','mode','msg1','msg2'],axis=1)

    # check log matches EEG markers
    n_offline_events = 25
    assert len(trials)-n_offline_events==len(e.raw.annotations) ,'missing {0} hardware triggers'.format(len(trials)-n_offline_events-len(e.raw.annotations))
    # if there are missing markers, we will have to use timings from the log file.
    #TODO separate script to write compatible events files for the various ways the triggers could be f%$*ed 

    # load condition key for language localiser & give a numeric ID to each condition
    langloc_stim = pd.read_table('langloc_key.txt', sep='|',usecols=['Condition','Stim'])
    langloc_recode_dict = dict(zip(np.unique(langloc_stim.Condition),[31,32,33,34]))
    langloc_stim["value"] = langloc_stim.Condition.apply(lambda x: langloc_recode_dict[x])

    ev_df = trials.iloc[0:-n_offline_events].loc[:,['event','value']] # select online events
    
    # reassign event value based on langloc_stim
    def recode(series):
        trial = series.event
        orig_value = series.value
        for sentence, new_value in langloc_stim.loc[:,['Stim','value']].itertuples(index=False):
            if re.search(sentence,trial):
                return new_value
        return orig_value

    ev_df['value'] = ev_df.loc[:,['event','value']].apply(recode, axis=1)


    ev_dict = dict(zip(ev_df.event, ev_df.value)) # make dict of event desc to values

    e.raw.annotations.description = ev_df.loc[:,'event'].to_numpy(dtype=str) # overwrite default uninformative annotations
    # extract events from annotations
    events, event_id =  mne.events_from_annotations(e.raw,ev_dict) # get mne events
    # plot raw
    e.raw_copy.plot(bad_color='r', title="raw",events=events)
    # plot preprocessed
    e.raw.plot(bad_color='r',title="preprocessed",events=events)
    plt.show(block=False)
    
    # epoch the data, selecting only lang localiser
    lang_eventID = {key: value for key, value in ev_dict.items() if value >30} # get dict of localiser events
    e.lang_epochs = mne.Epochs(e.raw, events, lang_eventID, tmin = -.2, tmax= 1)
    # 
    e.lang_erps = {cond: e.lang_epochs[langloc_recode_dict[cond]].average() for cond in langloc_recode_dict.keys()}
    
    e.av_lang = e.lang_epochs.average()
    e.av_lang.plot(exclude=[],spatial_colors=True)

# TODO
# interpolaate bad channels https://mne.tools/stable/generated/mne.io.Raw.html#mne.io.Raw.interpolate_bads 
    # add word onset events or EOG with https://mne.tools/stable/generated/mne.io.Raw.html#mne.io.Raw.add_events
    # access robust detrending and other external fn with https://mne.tools/stable/generated/mne.io.Raw.html#mne.io.Raw.apply_function
    # autoreject instead of ransac http://autoreject.github.io/ 