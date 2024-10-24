import os

# rename folders and files for when the experiment software was given the wrong pID
# best practice is to copy the misnamed folder from dropbox first to your local directory so if you f&$* it up you haven't rekt the data

old_name = 'EML1_180A'
correct_name =  'EML1_180'

paths = (os.path.join(root, filename)
        for root, _, filenames in os.walk(r'/Users/roso8920/Dropbox (Emotive Computing)/EyeMindLink/Data')
        for filename in filenames)

for path in paths:
#   syntax: path.replace(old name, new name)
    newname = path.replace(old_name, correct_name)
    if newname != path:
        os.renames(path, newname)