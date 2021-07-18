/*
 * Macro to quantify the fluorescent signals in the nuclei for all channels
 * 
 * Brief workflow:
 * - Nuclei are detected using the StarDist convolutional neural network model for fluorescent nuclei.
 * - Nuclear ROIs can be filtered on size and eroded to prevent edge effects.
 * - Background signal in the measurement channel is measured (several options).
 * - N.B. Intensity values shown in the Results table are already background-subtracted.
 * 
 * Input: a folder containing 2D images with at least 2 channels. Multi-series microscopy format files are also supported.
 * 
 * Required Fiji update sites:
 * - StarDist
 * - CSBDeep
 * - CLIJ
 * - CLIJ2
 * - SCF MPI CBG
 * 
 * Author: Bram van den Broek, Netherlands Cancer Institute, December 2020 - March 2021
 * b.vd.broek@nki.nl
 *
 * N.B. This script heavily relies on (mostly) CLIJ2 (https://clij.github.io/) and StarDist (https://imagej.net/StarDist).
 * If you use this script in a publication, please cite them appropriately.
 * 
 * current version: 1.4
 * 
 * 
 * ----------------------------------
 * Changelog, from v0.9 onwards:
 * ----------------------------------
 * v1.0:
 * - Added downsampling possibility before StarDist (better for high resolution images)
 * - Nuclei numbers visible in output images 
 * 
 * v1.1:
 * - Added metrics: nucleus area & total intensity
 * - Added an option to exclude nuclei touching image edges
 * - Perform writing of nuclei numbers in output image with CLIJ2 (different method, faster)
 * - Added some visualization options
 *
 * v1.2, March 2021:
 * - Fully omit the RoiManager by:
 *   * Creating outlines (as overlay) using CLIJ
 *   * Writing cell numbers as overlay at the center of mass of the detected cells
 * - General improvements
 * 
 * v1.3, April 2021:
 * - Support for multi-series files (opened with Bio-Formats)
 * 
 *  * v1.4, July 2021:
 * - Fixed a critical bug! (Not clearing the GPU memory after pushCurrentSlice(). Consecutive channels are messed up.)
 * 
 */

#@ File (label = "Input folder", style = "directory") input
#@ File (label = "Output folder", style = "directory") output
#@ String(value="N.B. The output folder should NOT be a subfolder of the input folder.", visibility="MESSAGE") message
#@ String (label = "Process files with extension", value = ".tif") fileExtension
#@ Integer (label = "Nucleus marker channel", value = 1, min = 1) nucleiChannel
#@ Integer (label = "Pre-Stardist image downscale factor (>1 necessary for high resolution images)", value = 4, min = 1, description="The macro uses StarDist s pre-trained deep learning model for nuclei prediction. The network is not trained on high-resolution images; downscaling helps to correctly identify the nuclei.") downsampleFactor
#@ Float (label = "Pre-StarDist median filter radius (creates more round segmentation) (um)", value = 0.0, min = 0) medianRadius_setting
#@ Float (label = "StarDist nucleus probability threshold (0-1, higher is more strict)", value = 0.5, description="The macro uses StarDist's pre-trained deep learning model for nuclei prediction. Higher thresholds will eliminate less likely nuclei.") probabilityThreshold
#@ Boolean (label = "Exclude nuclei touching the edges of the image", value = false) excludeEdges
#@ Integer (label = "Remove nulei with diameter smaller than (um)", value = 4, min = 0) minNucleusSize_setting
#@ Integer (label = "Remove nulei with diameter larger than (um)", value = 40) maxNucleusSize_setting
#@ Float (label = "Shrink segmented nuclei with (units) - negative means expand", value = 0.5) shrinkSize_setting
#@ String(value="Background handling settings", visibility="MESSAGE") message2
#@ String (label = "Background subtraction method", choices={"Calculate value for every channel and image", "Automatic rolling ball", "Manual fixed value everywhere"}, style="listBox") background_subtraction
#@ Integer (label = "Rolling ball radius (units) (if applicable)", value = 100) rollingBallRadius_setting
#@ Integer (label = "Manual background value (if applicable)", value = 0) background
#@ String(value="Visualization options", visibility="MESSAGE") message3
#@ Boolean (label = "Save output images with overlayed nucleus outlines", value = false) saveImages
#@ Integer (label = "Opacity% of nucleus outlines overlay", value = 100, min=0, max=100) labelOpacity
#@ String (label = "Overlay display option", choices={"Cell numbers and outlines as overlay", "Cell numbers and outlines imprinted as pixels (RGB output)"}, style="listBox") overlayChoice
#@ Boolean (label = "Add cell numbers in output image", value = false) addNumbersOverlay
#@ Integer (label = "Cell numbers font size", value = 12, min = 1) labelFontSize
#@ ColorRGB(label = "Cell number font color", value="yellow") fontColor
#@ Boolean (label = "Hide images during processing", value = false) hideImages


//Settings outside the dialog
thickOutlines = false;			// Dilate the outlines, making them 3 pixels wide instead of 1 pixel)
backgroundPercentile = 0.05;	// For automatic background value per image calculation 
maxTileSize = 2000;				// Maximum StarDist tile size

var nrOfImages = 0;
var current_image_nr = 0;
var processtime = 0;
var nrNuclei = 0;
outputSubfolder = output;		// Initialize this variable

saveSettings();

run("Bio-Formats Macro Extensions");
run("Labels...", "color=white font=" + labelFontSize + " show draw");
run("Set Measurements...", "area mean median min stack redirect=None decimal=3");
run("Input/Output...", "jpeg=85 gif=-1 file=.tsv use_file copy_row save_column save_row");
if(nImages>0) run("Close All");
print("\\Clear");
run("Clear Results");
setBatchMode(true);

resultsTable = "All Results";
if(isOpen("Results_all_files.tsv")) close("Results_all_files.tsv");
if(!isOpen(resultsTable)) Table.create(resultsTable);
else Table.reset(resultsTable);
Table.showRowIndexes(true);

scanFolder(input);
processFolder(input);

selectWindow("Results");
run("Close");
selectWindow(resultsTable);
Table.rename(resultsTable, "Results");
saveAs("Results", output + File.separator + "Results_all_files.tsv");

restoreSettings;



// function to scan folders/subfolders/files to count files with correct fileExtension
function scanFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			scanFolder(input + File.separator + list[i]);
		if(endsWith(list[i], fileExtension))
			nrOfImages++;
	}
	if(nrOfImages == 0) exit("No files with extension '"+fileExtension+"' were found in "+input);
}


// function to scan folders/subfolders/files to find files with correct fileExtension
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i])) {
			outputSubfolder = output + File.separator + list[i];	
			if(!File.exists(outputSubfolder)) File.makeDirectory(outputSubfolder);	//create the output subfolder if it doesn't exist
			processFolder(input + File.separator + list[i]);
		}
		if(endsWith(list[i], fileExtension)) {
			current_image_nr++;
			showProgress(current_image_nr/nrOfImages);
			processFile(input, outputSubfolder, list[i]);
		}
	}
	//	print("\\Clear");
	print("\\Update1:Finished processing "+nrOfImages+" files.");
	print("\\Update2:Average speed: "+d2s(current_image_nr/processtime,1)+" images per minute.");
	print("\\Update3:Total run time: "+d2s(processtime,1)+" minutes.");
	print("\\Update4:-------------------------------------------------------------------------");

}

// The actual processing
function processFile(input, output, file) {
	run("Close All");

	starttime = getTime();
	print("\\Update1:Processing file "+current_image_nr+"/"+nrOfImages+": " + input + File.separator + file);
	print("\\Update2:Average speed: "+d2s((current_image_nr-1)/processtime,1)+" images per minute.");
	time_to_run = (nrOfImages-(current_image_nr-1)) * processtime/(current_image_nr-1);
	if(time_to_run<5) print("\\Update3:Projected run time: "+d2s(time_to_run*60,0)+" seconds ("+d2s(time_to_run,1)+" minutes).");
	else if(time_to_run<60) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. You'd better get some coffee.");
	else if(time_to_run<480) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes ("+d2s(time_to_run/60,1)+" hours). You'd better go and do something useful.");
	else if(time_to_run<1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. ("+d2s(time_to_run/60,1)+" hours). You'd better come back tomorrow.");
	else if(time_to_run>1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. This is never going to work. Give it up!");
	print("\\Update4:-------------------------------------------------------------------------");

	run("Bio-Formats Macro Extensions");	//Somehow necessary for every file(?)
	Ext.setId(input + File.separator + file);
	Ext.getSeriesCount(nr_series);

	if(endsWith(fileExtension, "tif") || endsWith(fileExtension, "jpg")) {	//Use standard opener
		open(input + File.separator + file);
		process_current_series(file);
	}
	else {	//Use Bio-Formats
		for(s = 0; s < nr_series; s++) {
			run("Close All");
			run("Bio-Formats Importer", "open=["+input + File.separator + file+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+s+1);
			seriesName = getTitle();
			seriesName = replace(seriesName,"\\/","-");	//replace slashes by dashes in the seriesName
			print(input + File.separator + seriesName);
	//		outputPath = output + File.separator + substring(seriesNa
			process_current_series(seriesName);
		}
	}
}


function process_current_series(image) {
	image = getTitle();
	Stack.setDisplayMode("grayscale");
	Stack.setChannel(nucleiChannel);
	if (hideImages == false) setBatchMode("show");
	run("Enhance Contrast", "saturated=0.35");

	getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pw, ph);
	minNucleusSize = PI*pow((minNucleusSize_setting / pw / 2),2);	//Calculate the nucleus area as if it were a circle
	maxNucleusSize = PI*pow((maxNucleusSize_setting / pw / 2),2);

	medianRadius = medianRadius_setting / pw;
	rollingBallRadius = rollingBallRadius_setting / pw;
	shrinkSize = shrinkSize_setting / pw;
	
	if(nucleiChannel > channels && current_image_nr==1) {
		showMessage("Warning: The image has less channels ("+channels+") than the selected nuclei channel ("+nucleiChannel+"). Proceeding with the default value (1).");
		nucleiChannel = 1;
	}
	detect_nuclei(image, nucleiChannel);
	labelmap = getLabelMaps_GPU(image, unit);
	roiManager("reset");

	//Calculate backgrounds
	selectWindow(image);
	backgrounds = newArray(channels);
	showStatus("Determining backgrounds...");
	if(background_subtraction == "Automatic rolling ball") {
		for(c=1; c<=channels; c++) {
			Stack.setChannel(c);
			showStatus("Subtracting background channel "+c+"...");
			showProgress(c, channels+1);
			run("Subtract Background...", "rolling="+rollingBallRadius+" stack");
			backgrounds[c-1] = "auto subtracted";
		}
	}
	else if(background_subtraction == "Calculate value for every channel and image") {
		for(c=1; c<=channels; c++) {
			backgrounds[c-1] = get_background(image, c, backgroundPercentile);
		}
	}
	else if(background_subtraction == "Manual fixed value everywhere") {
		for(c=1; c<=channels; c++) {
			backgrounds[c-1] = background;
		}
	}
	selectWindow(image);
	Stack.setChannel(nucleiChannel);
	run("Enhance Contrast", "saturated=0.35");

	measureROIs(image, labelmap, backgrounds);	//parameters: image, channel

	//Save the output image (+ cosmetics)
	filename = substring(file, 0, lastIndexOf(file, "."));
	if (saveImages == true) {
		selectWindow(image);
		setSlice(nucleiChannel);
		setMetadata("Label", "nucleus marker");
		if(overlayChoice == "Cell numbers and outlines imprinted as pixels (RGB output)") run("Flatten", "stack");
		saveAs("Tiff", output + File.separator + filename+"_analyzed");
	}
	endtime = getTime();
	processtime = processtime+(endtime-starttime)/60000;
}


function detect_nuclei(image, nucleiChannel) {
	selectWindow(image);
	run("Duplicate...", "duplicate channels=" + nucleiChannel + " title=nuclei");

	if(downsampleFactor != 0) {
		run("Duplicate...", "duplicate channels=" + nucleiChannel + " title=nuclei_downscaled");		
		run("Bin...", "x="+downsampleFactor+" y="+downsampleFactor+" bin=Average");	
	}
	if(medianRadius > 0) run("Median...", "radius="+medianRadius); 
	getDimensions(width, height, channels, slices, frames);

	starDistTiles = pow(floor(maxOf(width, height)/maxTileSize)+1,2);	//Determine the nr. of tiles
	// Run StarDist and output to the ROI manager (creating a label image works only when not operating in batch mode, and that is slower and more annoying.)
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'nuclei_downscaled', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.60000000000001', 'probThresh':'"+probabilityThreshold+"', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'"+starDistTiles+"', 'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");

	//Scale up again
	if(downsampleFactor != 0) {
		showStatus("Scaling ROIs...");
		for(i=0;i<roiManager("count");i++) {
			if(i%100==0) showProgress(i/roiManager("count"));
			roiManager("Select",i);
			run("Scale... ", "x="+downsampleFactor+" y="+downsampleFactor);
			roiManager("update");
		}
	}
	close("nuclei");
	close("nuclei_downscaled");
}


function getLabelMaps_GPU(image, unit) {
//Create a band around each nucleus using CLIJ2

	run("CLIJ2 Macro Extensions", "cl_device=");
	Ext.CLIJ2_clear();
	// In case another GPU needs to be selected:
	//Ext.CLIJ2_listAvailableGPUs();
	//availableGPUs = Table.getColumn("GPUName");
	//run("CLIJ2 Macro Extensions", "cl_device=" + availableGPUs[1]);

	//Create labelmap
	run("ROI Manager to LabelMap(2D)");
	run("glasbey_on_dark");
	labelmap_nuclei_raw = getTitle();
	Ext.CLIJ2_push(labelmap_nuclei_raw);

	//exclude labels on edges
	if(excludeEdges) Ext.CLIJ2_excludeLabelsOnEdges(labelmap_nuclei_raw, labelmap_nuclei);
	else labelmap_nuclei = labelmap_nuclei_raw;

	//Filter on area
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei, nucleiStarDist);	//count nuclei detected by StarDist
	run("Clear Results");
	Ext.CLIJ2_statisticsOfBackgroundAndLabelledPixels(labelmap_nuclei, labelmap_nuclei); //Somehow if you put (image, labelmap) as arguments the pixel count is wrong
	Ext.CLIJ2_pushResultsTableColumn(area, "PIXEL_COUNT");

	Ext.CLIJ2_excludeLabelsWithValuesOutOfRange(area, labelmap_nuclei, labelmap_nuclei_filtered, minNucleusSize, maxNucleusSize);
	Ext.CLIJ2_release(labelmap_nuclei);

	//Shrink nuclei
	if(shrinkSize > 0) {
		Ext.CLIJ2_shrinkLabels(labelmap_nuclei_filtered, labelmap_final, shrinkSize, false);
		Ext.CLIJ2_release(labelmap_nuclei_filtered);
	}
	if(shrinkSize < 0) {
		Ext.CLIJ2_dilateLabels(labelmap_nuclei_filtered, labelmap_final, shrinkSize);
		Ext.CLIJ2_release(labelmap_nuclei_filtered);
	}
	else if(shrinkSize == 0) labelmap_final = labelmap_nuclei_filtered;
	
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_final, nrNuclei);	//get the number of nuclei after filtering
	run("Clear Results");
	Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_final, labelmap_final_ordered);	//Renumber the cells from top to bottom
	Ext.CLIJ2_release(labelmap_final);
	Ext.CLIJ2_statisticsOfLabelledPixels(labelmap_final_ordered, labelmap_final_ordered); //Somehow if you put (image, labelmap) as arguments the pixel count is wrong
	print(image + " : " +nucleiStarDist+" nuclei detected by StarDist ; "+nucleiStarDist - nrNuclei+" nuclei with diameter outside ["+d2s(minNucleusSize_setting,0)+" - "+d2s(maxNucleusSize_setting,0)+"] range "+unit+" were removed.");

	Ext.CLIJ2_pull(labelmap_final_ordered);
	run("glasbey_on_dark");

	if(saveImages == true) {
		//Detect label outlines, dilate them and add as overlay to the original image
		Ext.CLIJ2_detectLabelEdges(labelmap_final_ordered, labelmap_edges);
		if(thickOutlines) {
			Ext.CLIJ2_dilateBox(labelmap_edges, labelmap_edges_dilated);
		}
		else labelmap_edges_dilated = labelmap_edges;
		if(thickOutlines) {	
			Ext.CLIJ2_dilateLabels(labelmap_final_ordered, labelmap_final_ordered_extended, 1);
		}
		else labelmap_final_ordered_extended = labelmap_final_ordered;
		Ext.CLIJ2_mask(labelmap_final_ordered_extended, labelmap_edges_dilated, labelmap_outlines);
		Ext.CLIJ2_release(labelmap_edges);
		if(thickOutlines) Ext.CLIJ2_release(labelmap_final_ordered_extended);
		if(thickOutlines) Ext.CLIJ2_release(labelmap_edges_dilated);
		Ext.CLIJ2_pull(labelmap_outlines);
		run("glasbey_on_dark");
		selectWindow(image);
		run("Add Image...", "image="+labelmap_outlines+" x=0 y=0 opacity="+labelOpacity+" zero");	//Add labelmap to image as overlay
		if(addNumbersOverlay) {
			setFont("SansSerif", labelFontSize, "antialiased");
			color = color_to_hex(fontColor);
			setColor(color);
			for (i = 0; i < nrNuclei; i++) {
				x = getResult("MASS_CENTER_X", i);
				y = getResult("MASS_CENTER_Y", i);
				Overlay.drawString(i+1, x - labelFontSize/2, y + labelFontSize/2);
			}
		}
	}

	return labelmap_final_ordered;
}

function get_background(image, channel, percentile) {
	selectWindow(image);
	Stack.setChannel(channel);
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	total = 0;
	bin=0;
	while (total < nPixels*percentile) {
		total += histogram[bin];
		bin++;
	} 
	setThreshold(0,bin-1);
	background = getValue("Median limit");
	resetThreshold();
	return background;
}


function measureROIs(image, labelmap, backgrounds) {
	//create data arrays
	file_name_image = newArray(nrNuclei);
	cell_nr_image = newArray(nrNuclei);
	nuc_area_image = newArray(nrNuclei);

	for(i=0; i<nrNuclei; i++) {
		file_name_image[i] = image;
		cell_nr_image[i] = i+1;
	}
	//measure intensity (only background corrected if 'Rolling ball' is selected!)
	selectWindow(image);

	//Measure, get the results and put in 'All Results' table
	if(nrNuclei != 0) {
		selectWindow(resultsTable);
		nRows = Table.size;
		for(i=nRows; i<nRows+nrNuclei; i++) {
			Table.set("file name", i, image, resultsTable);
			Table.set("cell nr", i, cell_nr_image[i-nRows], resultsTable);
		}
		for (c = 1; c <= channels; c++) {
			run("Clear Results");
			selectWindow(image);
			Stack.setChannel(c);
			run("Enhance Contrast", "saturated=0.35");
			run("Clear Results");
			Ext.CLIJ2_pushCurrentSlice(image);	//Or maybe copySlice if image is still in the GPU memory?
			Ext.CLIJ2_statisticsOfLabelledPixels(image, labelmap);
			Ext.CLIJ2_release(image);
			mean_image_channel = Table.getColumn("MEAN_INTENSITY", "Results");
			stdDev_image_channel = Table.getColumn("STANDARD_DEVIATION_INTENSITY", "Results");
			nuc_area_image = Table.getColumn("PIXEL_COUNT", "Results");
			
			for(i=nRows; i<nRows+nrNuclei; i++) {
				if(c==1) Table.set("Nucleus area ("+unit+"^2)", i, nuc_area_image[i-nRows]*pw*pw, resultsTable);	//in units^2
				if(background_subtraction == "Automatic rolling ball") {
					Table.set("Mean ch"+c, i, mean_image_channel[i-nRows], resultsTable);
				}
				else Table.set("Mean ch"+c, i, mean_image_channel[i-nRows] - backgrounds[c-1], resultsTable);
				Table.set("StdDev ch"+c, i, stdDev_image_channel[i-nRows], resultsTable);
				if(background_subtraction == "Automatic rolling ball") {
					Table.set("Total Int ch"+c, i, nuc_area_image[i-nRows]*pw*pw * (mean_image_channel[i-nRows]), resultsTable);
				}
				else Table.set("Total Int ch"+c, i, nuc_area_image[i-nRows]*pw*pw * (mean_image_channel[i-nRows] - backgrounds[c-1]), resultsTable);
				Table.set("Background ch"+c, i, backgrounds[c-1], resultsTable);
			}
		}
	}
	Table.update;
}


function color_to_hex(color) {
	colorArray = split(color,",,");
	hexcolor = "#" + IJ.pad(toHex(colorArray[0]),2) + IJ.pad(toHex(colorArray[1]),2) + IJ.pad(toHex(colorArray[2]),2);
	return hexcolor;
}
