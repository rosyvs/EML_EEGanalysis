from EML_EEGutils import *
import pandas as pd
import re
"""
EEG analysis for EML: preprocessing pipeline
"""

reprepro = True # re-do preprocessing. If False, checks outdir and skips participants with a prepro'd file
plotTF = False
datadir = os.path.normpath('/Volumes/Blue1TB/EyeMindLink/Data')
outdir = os.path.normpath('../../Data/EEG_processed/preprocessed_fif/')# save directory for processed EEG data
#triggerSources = pd.read_csv('triggerSources.csv') # available triggers (as checked by Rosy)
fnroot = 'EML1_'
participants = range(156,182)

for p in participants:
    pID= fnroot + '{:03d}'.format(p)  
    print_and_log("\n\n\n~~~~~pID: {0}~~~~~".format(pID))
    prepro_exists = False # default 
    use_streamed = False # default
    # check for existing preprocessed data first
    if os.path.isfile(os.path.normpath(outdir + '/' + pID + '_p.fif')):
        print_and_log("\n{0}: Existing preprocessed file found in {1}".format(pID,os.path.normpath(outdir + '/' + pID + '_p.fif')))
        prepro_exists = True 
        if reprepro:
            print_and_log("{0}: Will overwrite existing preprocessed file ...".format(pID))
   
    # (Re)-run prepro
    if reprepro or not prepro_exists: 
        print_and_log("\n{0}: Preprocessing pipeline started...".format(pID))
        pIDpath = os.path.normpath(datadir + "/" + pID + "/EEG" )
        print_and_log("Searching folder {0}".format(pIDpath))
        # Choose EEG file to process: prefer SD card if triggers available, otherwise streamed

        print_and_log('Will first look for an SD card recording (LA*.vhdr)')
        vhdr = glob.glob(os.path.join(pIDpath,"LA*.vhdr")) # filenames starting 'LA' for SD card recording
        if not(vhdr):
            print_and_log('Attempted to find SD card recording (LA*.vhdr). No header file found!')
            use_streamed=True

        elif len(vhdr)>1:
            print_and_log('Multiple SD card header files found (LA*.vhdr), will only process first one. Please check data folder.')
            vhdr = vhdr[0]
        else:
            vhdr = vhdr[0]  # still need to unlist even if it is just a list of length 1
  
        if  use_streamed:
            print_and_log('Will look for a streamed recording (EML1_*.vhdr)')
            vhdr = glob.glob(pIDpath + "/EML1" +  "*.vhdr") # filenames starting 'EML1' for streamed recording
            if not(vhdr):
                print_and_log('No valid header file found for either streamed or SD recording! Skipping {0}'.format(pID))
                continue
            elif len(vhdr)>1:
                print_and_log('Multiple streamed header files found (EML1*.vhdr), will only process first one. Please check data folder.')
                vhdr = vhdr[0]
            else:
                vhdr = vhdr[0]  # still need to unlist even if it is just a list of length 1
        print_and_log('Will process the following EEG file(s): {0}'.format(vhdr))

        # preprocessing
        e = emleeg(data_dir='')
        e.loadFile(vhdr)
        e.extractInfo()
        e.zapline(nremove=4) # use meegkit to remove line noise with minimal distortion
            # improves upon e.raw.notch_filter(np.arange(60, 241, 60)) # powerline filter
        e.robustDetrend(order=11) # uses meegkit to detrend, removing slow drifts

        e.checkChannels() # uses pyprep
        e.prepro.set_montage('standard_1005')
        e.prepro.interpolate_bads(reset_bads=False) # interpolate bad channels based on neighbours but keep them marked as 'bad'
        e.prepro.set_eeg_reference() # ref to average

        # write preprocessed data to file & record which input data file was used
        e.prepro.save(fname=os.path.normpath(outdir + '/' + pID + '_p.fif'),overwrite=True)
        metadata_file = open(os.path.join(outdir, pID + '-info.txt'),'w')
        metadata_file.write(e.prepro.info['description'])

        if plotTF:
            # plt.ion()
            scalings=dict(eeg=50e-6) # y axis scaling
            # plot raw
            f1=e.raw.plot(bad_color='r', title="raw",start=600,scalings=scalings,show=False)
            # plot preprocessed
            f2=e.prepro.plot(bad_color='r',title="preprocessed",start=600,scalings=scalings,show=False)
            # plot PSD raw
            f3=e.raw.plot_psd(fmax=100, picks = e.ch_names_eeg,show=False)
            # plot PSD preprocessed
            f4=e.prepro.plot_psd(fmax=100, picks = e.ch_names_eeg,show=False)

            f1.figsize=[10,6.4]
            f1.savefig(os.path.join(outdir, pID +"-raw.png"))
            f2.figsize=[10,6.4]
            f2.savefig(os.path.join(outdir, pID +"-preprocessed.png"))
            f3.savefig(os.path.join(outdir, pID +"-rawPSD.png"))
            f4.savefig(os.path.join(outdir, pID +"-preprocessedPSD.png"))
            # plt.close('all')
    else:
        print_and_log("{0}: Skipping preprocessing pipeline.".format(pID))

# TODO
    # access robust detrending and other external fn with https://mne.tools/stable/generated/mne.io.Raw.html#mne.io.Raw.apply_function
    # autoreject instead of pyprep? http://autoreject.github.io/ 

# %%
