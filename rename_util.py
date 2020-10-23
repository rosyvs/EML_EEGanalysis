import os

paths = (os.path.join(root, filename)
        for root, _, filenames in os.walk(r'C:\Users\roso8920\Dropbox (Emotive Computing)\EML Rosy\Data\eml1_058-')
        for filename in filenames)

for path in paths:
    # the '#' in the example below will be replaced by the '-' in the filenames in the directory
    newname = path.replace('eml1_058-', 'EML1_058')
    if newname != path:
        os.renames(path, newname)