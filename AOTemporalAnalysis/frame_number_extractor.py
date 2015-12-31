# Robert Cooper
# 12-30-2015

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
	text = open(file, 'rb').read().replace('\r\n', '\n')
	open(file, 'wb').write(text)
	
	fhandle = open(file)
	pick = pickle.load(fhandle)
        fhandle.close()	

	# Get the frames we used, and save them to a file.
	np.savetxt(file[0:-4]+ "_acceptable_frames.csv",pick['acceptable_frames'],delimiter=',',fmt='%f')

     


