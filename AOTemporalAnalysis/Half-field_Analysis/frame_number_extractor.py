# Robert Cooper
# 12-30-2015
#
# This extracts the frame indexes that were good enough to be registered from a *.dmp file,
# And throws them into an [imname]_acceptable_frames.csv file so they can be used in processing.
# 
# It requires:
#    * The .dmp file output from Alfredo Dubra's Demotion software suite. **I realize this makes it VERY specific-
#      I do not promise any amazing things happening as result of using this software!**
#    * The 'mat' file corresponding to the grid calibration- also using Alf Dubra's script.
#    * The dataset you wish to put through the temporal analysis pipeline.
#
#

import os,pickle,csv, sys
import Tkinter as tk
import Tkconstants, tkFileDialog
import numpy as np

root = tk.Tk()

options = {}
options['title'] = 'Select the dmp file folder'
options['parent'] = root

folder_path = tkFileDialog.askdirectory(**options)

root.destroy()

for file in os.listdir(folder_path):
	if file.endswith(".dmp"):
		print file
		#Fix the fact that it was done originally in Windows...
		text = open(os.path.join(folder_path,file), 'rb').read().replace('\r\n', '\n')
		open(os.path.join(folder_path,file), 'wb').write(text)
	
		fhandle = open(os.path.join(folder_path,file))
		pick = pickle.load(fhandle)
		fhandle.close()	

		# Get the frames we used, and save them to a file.
		np.savetxt(os.path.join(folder_path,file[0:-4]+ "_acceptable_frames.csv"),pick['acceptable_frames'],delimiter=',',fmt='%f')

     


