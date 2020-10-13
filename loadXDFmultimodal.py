# %%

import pyxdf as xd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os



data_fpath = "../../Data/XDF/sub-rosy_ses-S001_task-test_run-001_eeg.xdf"

streams, header = xd.load_xdf(data_fpath)

print("Found {} streams:".format(len(streams)))
for ix, stream in enumerate(streams):
    y = stream['time_series']
    print("Stream {}: {} - type {} - uid {} - shape {} at {} Hz (effective {} Hz)".format(
            ix + 1, stream['info']['name'][0],
            stream['info']['type'][0],
            stream['info']['uid'][0],
            (int(stream['info']['channel_count'][0]), len(stream['time_stamps'])),
            stream['info']['nominal_srate'][0],
            stream['info']['effective_srate'])
        )
        if any(stream['time_stamps']):
            print("\tDuration: {} s".format(stream['time_stamps'][-1] - stream['time_stamps'][0]))
    if isinstance(y, list):
        # list of strings, draw one vertical line for each marker
        for timestamp, marker in zip(stream['time_stamps'], y):
            plt.axvline(x=timestamp)
            print(f'Marker "{marker[0]}" @ {timestamp:.2f}s')
    elif isinstance(y, np.ndarray):
        # numeric data, draw as lines
        plt.plot(stream['time_stamps'], y)
    else:
        raise RuntimeError('Unknown stream format')

    
print("Done.")
plt.show()        


# %%
