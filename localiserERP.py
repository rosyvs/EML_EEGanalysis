#%%
from EML_EEGutils import *
import pandas as pd
import re
from collections import defaultdict
import fnmatch
import matplotlib.pyplot as plt

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
participants = range(47,75) # recall that range is exclusive in Python 
group_sent_erps = defaultdict(list)
# load condition key for language localiser & give a numeric ID to each condition
langloc_stim = pd.read_table('langloc_key.txt', sep='|',usecols=['Condition','Stim'])
cond_str = np.unique(langloc_stim.Condition)
LLcond_to_value = dict(zip(cond_str,[31,32,33,34]))
LLvalue_to_cond = {v: k for k, v in LLcond_to_value.items()} 
langloc_stim["value"] = langloc_stim.Condition.apply(lambda x: LLcond_to_value[x])

# I want to give each lang loc condition a different value
# reassign event value based on searching for event string in langloc_stim
def recode(series):
    trial = series.EVENT
    orig_value = series.VAL
    for sentence, new_value in langloc_stim.loc[:,['Stim','value']].itertuples(index=False):
        if re.search(sentence,trial):
            return new_value
    return orig_value
#%%
for p in participants:
    pID= fnroot + '{:03d}'.format(p)  
    print_and_log("\n\n\n~~~~~pID: {0}~~~~~".format(pID))
    prepro_exists = False
    # check for existing preprocessed data first
    if os.path.isfile(os.path.normpath(preprodir + '/' + pID + '_p.fif')):
        print_and_log("\n\n\n{0}: Existing preprocessed file found in {1}"\
            .format(pID,os.path.normpath(preprodir + '/' + pID + '_p.fif')))
        prepro_exists = True 
    else:
        print_and_log("\n\n\n{0}: Skipping analysis. No preprocessed file found in {1}"\
            .format(pID,os.path.normpath(preprodir + '/' + pID + '_p.fif')))
        continue

    print_and_log("\n\n\nComputing Language Localiser response for {0}".format(pID))

    e = emleeg(data_dir='')
    e.prepro = mne.io.read_raw_fif(fname=os.path.normpath(preprodir + '/' + pID + '_p.fif'),preload=True)

    # get events
    event_file = os.path.normpath(rawdir + '/' + pID + '/EEG/events.csv')
    trials = pd.read_csv(event_file)

    # detect whether preprocessed EEG came from streamed or SD card and use the right event samples
    if fnmatch.fnmatch(e.prepro.info['description'],'*LA*.vhdr'): # this was preferred option during prepro
        use_col = 'eegSD_sample_est'
        print_and_log('using EEG data and triggers from SD card')
    elif fnmatch.fnmatch(e.prepro.info['description'],'*EML1_???.vhdr'): 
        use_col = 'eeg_sample_est'
        print_and_log('using EEG data and triggers from streamed recording')
    else:
        print_and_log('***WARNING*** \nNo raw EEG filename recorded in preprocessed file!\nCan''t detect which triggers to use. Skipping {0}'.format(pID))
        continue
    
    e.prepro.set_eeg_reference() # ref to average
    e.prepro.filter(l_freq=.1,h_freq=None)


    trials['VAL'] = trials.loc[:,['EVENT','VAL']].apply(recode, axis=1)  
    # replace long event description (full sentence) with short condition ID
    trials=trials.assign(event=trials.VAL.map(LLvalue_to_cond).fillna(trials.MSG))

    events = trials.loc[:,[use_col,'VAL']]
    events = events.dropna(axis=0)# remove rows with no corresponding EEG sample
    ev_dict = dict(zip(trials.loc[~pd.isna(trials[use_col]),'EVENT'], \
        trials.loc[~pd.isna(trials[use_col]),'VAL'])) # make dict of event desc to values

    events.insert(loc=1,column='output',value=0)
# realign language localiser triggers to start at stimulus, not ISI, by adding 2000 ms to trigger
    events.loc[events['VAL'].isin(cond_str),use_col] =\
    2000+events.loc[events['VAL'].isin(cond_str),use_col]
    events=events.to_numpy().astype(int)                             
    
    # params for rejecting epochs based on amplitude
    reject = dict(eeg=150e-6)
                                                                                                                                                                                
    # epoch the data, selecting only lang localiser epochs
    e.sent_epochs = mne.Epochs(e.prepro, events, LLcond_to_value, tmin = -.2, tmax= 3, on_missing='raise',preload=True)
epochs.drop_bad()


    e.sent_erps = {cond: e.sent_epochs[cond].average(picks='eeg') for cond in LLcond_to_value.keys()}
    #mne.viz.plot_evoked_topo(list(e.sent_erps.values())) 

    cond_counts = [len(e.sent_epochs[cond]) for cond in cond_str]
    # TODO detect whether any conditions are missing, if so, skip this subject
    skip_p=False
    for ci in range(len(cond_counts)):
        if cond_counts[ci]<10:
            print_and_log('Only {0} epochs for condition {1}'.format(cond_counts[ci],cond_str[ci]))
            skip_p = True
    if not skip_p:
        for c in LLcond_to_value.keys():
            group_sent_erps[c].append(e.sent_erps[c])
        e.sent_av = e.sent_epochs.average()
    else: 
        print_and_log('Skipping {0}'.format(pID))
    #e.sent_av.plot(exclude=[],spatial_colors=True)
#%%
ci=0
gav_sent_erps={}
for c in cond_str:
    # detect nans in data
    this = group_sent_erps[c]
    si=0
    for s in this:
        data = s.to_data_frame()
        if data.isnull().values.any():
            print('Missing {0} values in condition {1}, subject {2}'.format(data.isnull().sum(),c,participants[si]))
        si=si+1

    gav_sent_erps[c] = (mne.grand_average(group_sent_erps[c], interpolate_bads = False, drop_bads=False).apply_baseline((-0.2,0)))
    ci+=1
fig, ax = plt.subplots()
gavplot = mne.viz.plot_compare_evokeds(gav_sent_erps,axes = ax)
ax=plt.gca()
ax.set_xlim(-0.2, 1.6) 

# %%
