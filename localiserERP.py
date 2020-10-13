from EML_EEGutils import *

# EEG analysis for EML: find localiser-evoked potentials 

datadir = os.path.normpath('C:/Users/roso8920/Dropbox (Emotive Computing)/EyeMindLink/Data')
fnroot = 'EML1_'
participants = [52]
# TODO automatic generation of file list from particpants
p=0
pID= fnroot + '{:03d}'.format(participants[p])  
print("~~~~~Processing pID~: {0}".format(pID))
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
    print(e.raw_copy)
    e.robustDetrend(order=5)
    e.checkChannels()
    print(e.raw)

        # plot raw
    e.raw_copy.plot(bad_color='r', title="raw")
    # plot preprocessed
    e.raw.plot(bad_color='r',title="preprocessed")
    plt.ion()
    plt.show()
