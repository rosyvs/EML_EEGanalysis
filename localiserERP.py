#%%
from EML_EEGutils import *
import pandas as pd
import re
from collections import defaultdict
"""
EEG analysis for EML: find language localiser ERP
Stimuli after FEdorenko et al. 2010
Sentences of 8 words shown at fixation, RSVP
Sentence types: words, jabberwocky words, sentences, jabberwocky sentences 
100ms per word || 2 sec between sentences || 10 sec between blocks
"""

rawdir = os.path.normpath('C:/Users/roso8920/Dropbox (Emotive Computing)/EyeMindLink/Data')
preprodir = os.path.normpath('../../Data/EEG_processed/')# save directory for processed EEG data
fnroot = 'EML1_'
participants = range(55,72) # recall that range is exclusive in Python 

group_sent_erps = defaultdict(list)
#%%
for p in participants:
    pID= fnroot + '{:03d}'.format(p)  
    print_and_log("\n\n\n~~~~~pID: {0}~~~~~".format(pID))
    prepro_exists = False
    # check for existing preprocessed data first
    if os.path.isfile(os.path.normpath(preprodir + '/' + pID + '_p.fif')):
        print_and_log("\n\n\n{0}: Existing preprocessed file found in {1}".format(pID,os.path.normpath(preprodir + '/' + pID + '_p.fif')))
        prepro_exists = True 
    else:
        print_and_log("\n\n\n{0}: Skipping analysis. No preprocessed file found in {1}".format(pID,os.path.normpath(preprodir + '/' + pID + '_p.fif')))
        continue

    print_and_log("\n\n\nComputing Language Localiser response for {0}".format(pID))

    e = emleeg(data_dir='')
    e.prepro = mne.io.read_raw_fif(fname=os.path.normpath(preprodir + '/' + pID + '_p.fif'),preload=True)

    # get events
    log_file = os.path.normpath(rawdir + "/" + pID + "/" + pID + '_Trials.txt')
    log_cols = ['date','time','mode','msg1','msg2','event','value']
    trials = pd.read_csv(log_file, sep='\t',skiprows=1, header=None,names=['col1','event','value'])
    trials[['date','time','mode','msg1','msg2']]=trials.col1.str.split(expand=True)
    trials=trials.drop(['col1','mode','msg1','msg2'],axis=1)

    # check log matches EEG markers
    n_offline_events = 25
    if not len(trials)-n_offline_events==len(e.prepro.annotations):
        print_and_log('***WARNING***\n Trial log contains {0} online trials, but .vmrk contains {1} markers.'.format(len(trials)-n_offline_events,len(e.prepro.annotations)))
        print_and_log('Missing {0} hardware triggers from .vmrk'.format(len(trials)-n_offline_events-len(e.prepro.annotations)))
        print_and_log('Skipping {0}'.format(pID))
        continue
    # if there are missing markers, we will have to use timings from the log file.
    #TODO separate script to write compatible events files for the various ways the triggers could be f%$*ed 

    e.prepro.set_eeg_reference() # ref to average
    e.prepro.filter(l_freq=.1,h_freq=None)

    # load condition key for language localiser & give a numeric ID to each condition
    langloc_stim = pd.read_table('langloc_key.txt', sep='|',usecols=['Condition','Stim'])
    LLcond_to_value = dict(zip(np.unique(langloc_stim.Condition),[31,32,33,34]))
    LLvalue_to_cond = {v: k for k, v in LLcond_to_value.items()} 
    langloc_stim["value"] = langloc_stim.Condition.apply(lambda x: LLcond_to_value[x])

    ev_df = trials.iloc[0:-n_offline_events].loc[:,['event','value']] # select online events
    
    # reassign event value based on searching for event string in langloc_stim
    def recode(series):
        trial = series.event
        orig_value = series.value
        for sentence, new_value in langloc_stim.loc[:,['Stim','value']].itertuples(index=False):
            if re.search(sentence,trial):
                return new_value
        return orig_value
    ev_df['value'] = ev_df.loc[:,['event','value']].apply(recode, axis=1)  
    # replace long event description (full sentence) with short condition ID
    ev_df=ev_df.assign(event=ev_df.value.map(LLvalue_to_cond).fillna(ev_df.event))

    ev_dict = dict(zip(ev_df.event, ev_df.value)) # make dict of event desc to values

    e.prepro.annotations.description = ev_df.loc[:,'event'].to_numpy(dtype=str) # overwrite default uninformative annotations
    # extract events from annotations
    events, event_id =  mne.events_from_annotations(e.prepro,ev_dict) # get mne events
    
    # epoch the data, selecting only lang localiser epochs
    # lang_eventID = {key: value for key, value in ev_dict.items() if value >30} # get dict of localiser events
    e.sent_epochs = mne.Epochs(e.prepro, events,LLcond_to_value, tmin = -.2, tmax= 4)
    e.sent_erps = {cond: e.sent_epochs[cond].average() for cond in LLcond_to_value.keys()}
    #mne.viz.plot_evoked_topo(list(e.sent_erps.values())) 

    e.sent_av = e.sent_epochs.average()
    #e.sent_av.plot(exclude=[],spatial_colors=True)

    for c in LLcond_to_value.keys():
        group_sent_erps[c].append(e.sent_erps[c])

#%%
ci=0
gav_sent_erps={}
for c in LLcond_to_value.keys():
    gav_sent_erps[c] = (mne.grand_average(group_sent_erps[c]).apply_baseline((1.8,2)))
    ci+=1
mne.viz.plot_compare_evokeds(gav_sent_erps, xlim=(1.8,3))
#TODO
# add word onset events or EOG with https://mne.tools/stable/generated/mne.io.Raw.html#mne.io.Raw.add_events
# %%
