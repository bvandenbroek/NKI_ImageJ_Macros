// Macro to quantify the size parameters (e.g. circularity) of cells in a movie
//
// Bram van den Broek, the Netherlands Cancer Institute, 2017-2018
// b.vd.broek@nki.nl


#@ Integer (label = "Channel to analyze", style = "spinner", min=0, max=5, value=1) ch
#@ Integer (label = "Lower size limit", style = "spinner", min=0, max=1000, value=30) lower_size_limit
#@ Integer (label = "Upper size limit", style = "spinner", min=0, max=1000, value=300) upper_size_limit
#@ Boolean (label = "Apply median filter", value=false) median_filter
#@ Float (label = "Median filter radius (pixels)", style = "spinner", min=0, max=5, value=0.5) median_radius
#@ String(label = "Threshold method", value="Yen", choices={"Default", "Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"}, style="listbox") threshold_method
#@ File (label = "Output directory", style = "directory") savedir

close("\\Others");
run("Set Measurements...", "area perimeter shape stack redirect=None decimal=3");
run("Clear Results");


setBatchMode(true);

name = getTitle();
dir = getDirectory("image");
run("Duplicate...", "duplicate channels="+ch);
if(median_filter==true) run("Median 3D...", "x="+median_radius+" y="+median_radius+" z="+median_radius);

getDimensions(width, height, channels, slices, frames);
circularity = newArray(frames);
nr_cells = newArray(frames);
if(slices>1) run("Z Project...", "projection=[Max Intensity] all");
rename("channel_"+ch);
run("Grays");
run("Enhance Contrast", "saturated=0.35");
setBatchMode("show");
setTool("freehand");

waitForUser("Select region(s) that will not be included in the analysis (e.g. autofluorescence). Use the SHIFT key to select multiple regions. Press OK to continue.");
if(selectionType!=-1) run("Clear", "stack");
run("Select None");

run("Duplicate...", "duplicate title=mask");
run("Enhance Contrast", "saturated=0.35");

//Apply the selected automatic threshold method and analyze the cells
run("Convert to Mask", "method="+threshold_method+" background=Dark calculate black list");
run("Analyze Particles...", "size="+lower_size_limit+"-"+upper_size_limit+" show=[Bare Outlines] exclude record add display summarize stack");
rename("outlines");

run("Invert", "stack");
close("filtered");
close("mask");
run("Merge Channels...", "c1=channel_"+ch+" c2=outlines create");

Stack.setChannel(2);
run("Green");
Stack.setChannel(1);
run("Grays");
run("Enhance Contrast", "saturated=0.35");
setBatchMode("exit and display");

if (median_filter==true) {
	saveAs("Tiff", savedir+File.separator+name+" - median_"+median_radius+" - outlines");
	saveAs("Results", savedir+File.separator+name+" - median_"+median_radius+" - results_all.txt");
}
else {
	saveAs("Tiff", savedir+name+" - outlines");
	saveAs("Results", savedir+name+" - results_all.txt");
}
selectWindow("Results");
run("Close");

//Retreive data from the summary and save it
selectWindow("Summary of mask");
for(i=0;i<frames;i++) {
	nr_cells[i] = getResult("Count",i);
	circularity[i] = getResult("Circ.",i);
}
saveAs("Results", savedir+File.separator+name+" - median_"+median_radius+" - summary.txt");
run("Close");

//Create plots with nr. of detected cells and average circularity
Plot.create("Cells plot", "time (frames)", "number of detected cells");
Plot.setColor("blue");
Plot.add("line", nr_cells);
Plot.show;
saveAs("Tiff", savedir+File.separator+name+" - nr_cells PLOT");

Plot.create("Circularity plot", "time (frames)", "average circularity");
Plot.setColor("red");
Plot.add("line", circularity);
Plot.show;
saveAs("Tiff", savedir+File.separator+name+" - circularity PLOT");
