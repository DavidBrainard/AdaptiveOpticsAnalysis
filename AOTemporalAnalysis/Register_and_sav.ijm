fname= getTitle()
fpath=getInfo("image.directory")

run("StackReg", "transformation=Affine");
run("Z Project...", "projection=[Average Intensity]");
saveAs("Tiff", fpath + substring(fname, 0, lengthOf(fname)-4 )+"_affine_AVG.tif");
selectWindow(fname);
run("AVI... ", "compression=None frame=30 save="+fpath + substring(fname, 0, lengthOf(fname)-4 )+"_affine.avi");
