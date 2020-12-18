/*
 * Macro to quantify the fluorescent signals in the nuclei and cytoplasm in 2D images, and compute the ratios.
 * 
 * Brief workflow:
 * - Nuclei are detected using the StarDist convolutional neural network model for fluorescent nuclei.
 * - Cytoplasm is defined by a band with a certain size around the nuclei ROIs.
 * - Overlap between nuclei and/or cytoplasm of neighboring cells is automatically removed
 * - Background signal in the measurement channel is measured or subtracted (several options).
 * - Background-subtracted values and ratios are shown in the Results table.
 * 
 * Author: Bram van den Broek, Netherlands Cancer Institute, October 2020
 * b.vd.broek@nki.nl
 * 
 * Input: a folder containing 2D images with at least 2 channels.
 * 
 * Required Fiji update sites:
 * - StarDist
 * - CSBDeep
 * - CLIJ
 * - CLIJ2
 * - SCF MPI CBG
 * - IJPB-plugins (MorphoLibJ)
 * 
 * N.B. This script heavily relies on (mostly) CLIJ2 (https://clij.github.io/) and StarDist (https://imagej.net/StarDist).
 * If you use this script in a publication, please cite them appropriately.
 * 
 */


#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "Process files with extension", value = ".tif") fileExtension
#@ Integer (label = "Nuclei channel", value = 1) nucleiChannel
#@ Integer (label = "Measurement channel", value = 2) measurementChannel
#@ Float (label = "Pre-StarDist median filter radius (units (e.g. um))", value = 0.99) medianRadius_setting
#@ Float (label = "StarDist nucleus probability threshold", value = 0.5) probabilityThreshold
#@ Integer (label = "Remove nulei with diameter smaller than (units)", value = 4) minNucleusSize_setting
#@ Integer (label = "Remove nulei with diameter larger than (units)", value = 40) maxNucleusSize_setting
#@ Float (label = "Gap size between nucleus and cytoplasm (units)", value = 1) gapSize_setting
#@ Float (label = "Band size around the nucleus (units)", value = 1.5) bandSize_setting
#@ String (choices={"Automatic rolling ball","Automatic value per image","Manual fixed value for all images"}, style="listBox") background_subtraction
#@ Integer (label = "Rolling ball radius (units) (if applicable)", value = 50) rollingBallRadius_setting
#@ Integer (label = "Manual background value (if applicable)", value = 0) background
#@ Boolean (label = "Hide images during processing", value = false) hideImages
#@ Boolean (label = "Invert ratio in output?", value = false) invertRatio
#@ Boolean (label = "Save analyzed image?", value = false) saveImages

backgroundPercentile = 0.05;	// For Automatic value per image calculation 
maxTileSize = 2000;				// Maximum StarDist tile size


var nrOfImages = 0;
var current_image_nr = 0;
var processtime = 0;
var nrNuclei = 0;
outputSubfolder = output;	//initialize this variable

saveSettings();

run("Set Measurements...", "area mean median min stack redirect=None decimal=3");
run("Input/Output...", "jpeg=85 gif=-1 file=.tsv use_file copy_row save_column save_row");
if(nImages>0) run("Close All");
print("\\Clear");
run("Clear Results");
roiManager("reset");
roiManager("Show all without labels");
setBatchMode(true);

resultsTable = "All Results";
if(isOpen("Results_all_files.tsv")) close("Results_all_files.tsv");
if(!isOpen(resultsTable)) Table.create(resultsTable);
else Table.reset(resultsTable);
Table.showRowIndexes(true);
//startTime = getTime();

scanFolder(input);
processFolder(input);

selectWindow("DAPI-intensity-measurements");
run("Close");
selectWindow("Measurement-intensity-measurements");
run("Close");
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


function processFile(input, output, file) {
	run("Close All");

	starttime = getTime();
	print("\\Update1:Processing file "+current_image_nr+"/"+nrOfImages+": " + input + file);
	print("\\Update2:Average speed: "+d2s((current_image_nr-1)/processtime,1)+" images per minute.");
	time_to_run = (nrOfImages-(current_image_nr-1)) * processtime/(current_image_nr-1);
	if(time_to_run<5) print("\\Update3:Projected run time: "+d2s(time_to_run*60,0)+" seconds ("+d2s(time_to_run,1)+" minutes).");
	else if(time_to_run<60) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. You'd better get some coffee.");
	else if(time_to_run<480) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes ("+d2s(time_to_run/60,1)+" hours). You'd better go and do something useful.");
	else if(time_to_run<1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. ("+d2s(time_to_run/60,1)+" hours). You'd better come back tomorrow.");
	else if(time_to_run>1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. This is never going to work. Give it up!");
	print("\\Update4:-------------------------------------------------------------------------");

	open(input + File.separator + file);
	Stack.setDisplayMode("grayscale");
	if (hideImages == false) setBatchMode("show");
	run("Enhance Contrast", "saturated=0.35");

	getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pw, ph);
	minNucleusSize = PI*pow((minNucleusSize_setting / pw / 2),2);	//Calculate the nucleus area as if it were a circle
	maxNucleusSize = PI*pow((maxNucleusSize_setting / pw / 2),2);
	medianRadius = medianRadius_setting / pw;
	gapSize = gapSize_setting / pw;
	bandSize = bandSize_setting / pw;
	rollingBallRadius = rollingBallRadius_setting / pw;

	image = getTitle();
	detect_nuclei(image, nucleiChannel);
	labelmapNames = getLabelMaps_GPU(image, unit);

	if(background_subtraction == "Automatic rolling ball") {
		selectWindow(image);
		run("Subtract Background...", "rolling="+rollingBallRadius+" stack");
		background = 0;	//overwrite manual background variable
	}
	else if(background_subtraction == "Automatic value per image") {
		background = get_background(image, measurementChannel, backgroundPercentile);
	}
	selectWindow(image);
	Stack.setChannel(nucleiChannel);
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(measurementChannel);
	run("Enhance Contrast", "saturated=0.35");

	measureROIs(image, nucleiChannel, measurementChannel, labelmapNames, background);	//parameters: image, channel

	//Save the image
	filename = substring(file, 0, lastIndexOf(file, "."));
	if (saveImages == true) {
		selectWindow(image);
		run("Flatten", "stack");
		saveAs("Tiff", output + File.separator + filename+"_analyzed");
	}

	endtime = getTime();
	processtime = processtime+(endtime-starttime)/60000;
}


function detect_nuclei(image, nucleiChannel) {
	selectWindow(image);
	run("Duplicate...", "duplicate channels=" + nucleiChannel + " title=nuclei");
	if(medianRadius > 0) run("Median...", "radius="+medianRadius+" stack");

	getDimensions(width, height, channels, slices, frames);
	starDistTiles = pow(floor(maxOf(width, height)/maxTileSize)+1,2);	//Determine the nr. of tiles
	//run StarDist and output to the ROI manager (creating a label image works only when not in batch mode, and this is slower)
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'nuclei', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.60000000000001', 'probThresh':'"+probabilityThreshold+"', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'"+starDistTiles+"', 'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	close("nuclei");
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
	labelmap_nuclei = getTitle();

	Ext.CLIJ2_push(image);	//Still necessary?
	Ext.CLIJ2_push(labelmap_nuclei);

	//Filter on area
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei, nucleiStarDist);	//count nuclei detected by StarDist
	run("Clear Results");
	roiManager("reset");
	Ext.CLIJ2_statisticsOfBackgroundAndLabelledPixels(labelmap_nuclei, labelmap_nuclei); //Somehow if you put (image, labelmap) as arguments the pixel count is wrong
	Ext.CLIJ2_pushResultsTableColumn(area, "PIXEL_COUNT");
	Ext.CLIJ2_release(image);	//Still necessary?

	Ext.CLIJ2_excludeLabelsWithValuesOutOfRange(area, labelmap_nuclei, labelmap_nuclei_filtered, minNucleusSize, maxNucleusSize);
	Ext.CLIJ2_release(labelmap_nuclei);

	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei_filtered, nrNuclei);	//get the number of nuclei after filtering
	
	//Create gap
	Ext.CLIJx_extendLabelsWithMaximumRadius(labelmap_nuclei_filtered, extendedlabelmap_gap, gapSize);

	//Expand labels further and subtract to create the rings
	Ext.CLIJx_extendLabelsWithMaximumRadius(extendedlabelmap_gap, extendedlabelmap_gapAndrings, bandSize);
	Ext.CLIJ2_subtractImages(extendedlabelmap_gapAndrings, extendedlabelmap_gap, labelmap_rings);
	Ext.CLIJ2_release(extendedlabelmap_gap);
	Ext.CLIJ2_release(extendedlabelmap_gapAndrings);
	
	print(image + " : " +nucleiStarDist+" nuclei detected by StarDist ; "+nucleiStarDist - nrNuclei+" nuclei with diameter outside ["+d2s(minNucleusSize_setting,0)+" - "+d2s(maxNucleusSize_setting,0)+"] range "+unit+" were removed.");

	//Overlay nuclei and rings with original image
	Ext.CLIJ2_pull(labelmap_nuclei_filtered);
	Ext.CLIJ2_pullBinary(labelmap_nuclei_filtered);
	rename("nuclei_mask");
	//Create a nice blue LUT for the nuclei overlay
    reds = newArray(256); 
    greens = newArray(256); 
    blues = newArray(256);
    for (i=0; i<256; i++) {
        reds[i] = 0;
        greens[i] = 0.5*i;
        blues[i] = i;
    }
    setLut(reds, greens, blues);

	Ext.CLIJ2_pull(labelmap_rings);
	run("glasbey_on_dark");
	selectWindow(image);
	
	run("Add Image...", "image=nuclei_mask x=0 y=0 opacity=25 zero");
	run("Add Image...", "image="+labelmap_rings+" x=0 y=0 opacity=50 zero");

//	Measure intensities on the GPU - CURRENTLY NOT USED; Measurements are performed with MorphoLibJ, because then also the median intensity can be measured.
//	Ext.CLIJ2_statisticsOfLabelledPixels(image, labelmap_nuclei_filtered);
//	nucleiMean = Table.getColumn("MEAN_INTENSITY", "Results");
//	Ext.CLIJ2_statisticsOfLabelledPixels(image, labelmap_rings);
	
	Ext.CLIJ2_clear();

	labelmapNames = newArray(labelmap_nuclei_filtered, labelmap_rings);
	return labelmapNames;
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
	//print("background = "+background);
	resetThreshold();
	return background;
}


function measureROIs(image, nucleiChannel, measurementChannel, labelmapNames, background) {
	//create data arrays
	file_name_image = newArray(nrNuclei);
	cell_nr_image = newArray(nrNuclei);
	DAPI_mean_image = newArray(nrNuclei);
	nuc_area_image = newArray(nrNuclei);
	nuc_mean_image = newArray(nrNuclei);
	nuc_median_image = newArray(nrNuclei);
	cyto_area_image = newArray(nrNuclei);
	cyto_mean_image = newArray(nrNuclei);
	cyto_median_image = newArray(nrNuclei);
	ratio_mean_image = newArray(nrNuclei);
	ratio_median_image = newArray(nrNuclei);
	background_array_image = newArray(nrNuclei);

	for(i=0; i<nrNuclei; i++) {
		file_name_image[i] = image;
		cell_nr_image[i] = i+1;
	}
	//measure nucleus intensity (only background corrected if 'Rolling ball' is selected!)
	selectWindow(image);
	nuclei = "DAPI";
	run("Duplicate...", "title="+nuclei+" duplicate channels="+nucleiChannel);
	run("Intensity Measurements 2D/3D", "input="+nuclei+" labels="+labelmapNames[0]+" mean stddev max min median numberofvoxels volume");	//Measure using MorphoLibJ
	DAPI_mean_image = Table.getColumn("Mean", ""+nuclei+"-intensity-measurements");
	close(nuclei);
	
	//measure the nuclei in the measurementChannel
	selectWindow(image);
	measurement = "Measurement";
	run("Duplicate...", "title="+measurement+" duplicate channels="+measurementChannel);
	run("Intensity Measurements 2D/3D", "input="+measurement+" labels="+labelmapNames[0]+" mean stddev max min median numberofvoxels volume");
	nuc_area_image = Table.getColumn("Volume", ""+measurement+"-intensity-measurements");
	nuc_mean_image = Table.getColumn("Mean", ""+measurement+"-intensity-measurements");
	nuc_mean_image = subtract_scalar_from_array(nuc_mean_image, background);
	nuc_median_image = Table.getColumn("Median", ""+measurement+"-intensity-measurements");
	nuc_median_image = subtract_scalar_from_array(nuc_median_image, background);

	//measure the rings in the measurementChannel
	run("Intensity Measurements 2D/3D", "input="+measurement+" labels="+labelmapNames[1]+" mean stddev max min median numberofvoxels volume");
	cyto_area_image = Table.getColumn("Volume", ""+measurement+"-intensity-measurements");
	cyto_mean_image = Table.getColumn("Mean", ""+measurement+"-intensity-measurements");
	cyto_mean_image = subtract_scalar_from_array(cyto_mean_image, background);
	cyto_median_image = Table.getColumn("Median", ""+measurement+"-intensity-measurements");
	cyto_median_image = subtract_scalar_from_array(cyto_median_image, background);

	close(measurement);

	//Calculate ratios
	for(i=0; i<nrNuclei; i++) {
		ratio_mean_image[i] = nuc_mean_image[i] / cyto_mean_image[i];
		ratio_median_image[i] = nuc_median_image[i] / cyto_median_image[i];
		if(background_subtraction != "Automatic rolling ball") background_array_image[i] = background;
		else background_array_image[i] = "auto";
		if(invertRatio == true) {
			ratio_mean_image[i] = 1/ratio_mean_image[i];
			ratio_median_image[i] = 1/ratio_median_image[i];
		}
	}

	//Handle the results
	if(nrNuclei != 0) {
		selectWindow(resultsTable);
		if(Table.size > 0) {
			//Get all results up to now
			file_name = Table.getColumn("file name", resultsTable);
			cell_nr = Table.getColumn("cell nr", resultsTable);
			DAPI_mean =  Table.getColumn("DAPI intensity", resultsTable);
			nuc_area = Table.getColumn("nuc area", resultsTable);
			cyto_area = Table.getColumn("cyto area", resultsTable);
			nuc_mean = Table.getColumn("nuc mean", resultsTable);
			cyto_mean = Table.getColumn("cyto mean", resultsTable);
			nuc_median = Table.getColumn("nuc median", resultsTable);
			cyto_median = Table.getColumn("cyto median", resultsTable);
			if(invertRatio == false) {
				ratio_mean = Table.getColumn("n/c ratio mean", resultsTable);
				ratio_median = Table.getColumn("n/c ratio median", resultsTable);
			}
			else {
				ratio_mean = Table.getColumn("c/n ratio mean", resultsTable);
				ratio_median = Table.getColumn("c/n ratio median", resultsTable);
			}
			background_array = Table.getColumn("background value", resultsTable);

			//Concatenate the results of the current image
			file_name = Array.concat(file_name, file_name_image);
			cell_nr = Array.concat(cell_nr, cell_nr_image);
			DAPI_mean = Array.concat(DAPI_mean, DAPI_mean_image);	
			nuc_area = Array.concat(nuc_area, nuc_area_image);
			cyto_area = Array.concat(cyto_area, cyto_area_image);
			nuc_mean = Array.concat(nuc_mean,nuc_mean_image);
			cyto_mean = Array.concat(cyto_mean,cyto_mean_image);
			nuc_median = Array.concat(nuc_median,nuc_median_image);
			cyto_median = Array.concat(cyto_median,cyto_median_image);
			ratio_mean = Array.concat(ratio_mean,ratio_mean_image);
			ratio_median = Array.concat(ratio_median,ratio_median_image);
			background_array = Array.concat(background_array, background_array_image);
		}
		else {	//first round
			Table.create("Results");
			file_name = file_name_image;
			cell_nr = cell_nr_image;
			DAPI_mean = DAPI_mean_image;
			nuc_area = nuc_area_image;
			cyto_area = cyto_area_image;
			nuc_mean = nuc_mean_image;
			cyto_mean = cyto_mean_image;
			nuc_median = nuc_median_image;
			cyto_median = cyto_median_image;
			ratio_mean = ratio_mean_image;
			ratio_median = ratio_median_image;
			background_array = background_array_image;
		}
		Table.setColumn("file name", file_name, resultsTable);
		Table.setColumn("cell nr", cell_nr, resultsTable);
		Table.setColumn("DAPI intensity", DAPI_mean, resultsTable);
		Table.setColumn("nuc area", nuc_area, resultsTable);
		Table.setColumn("cyto area", cyto_area, resultsTable);
		Table.setColumn("nuc mean", nuc_mean, resultsTable);
		Table.setColumn("cyto mean", cyto_mean, resultsTable);
		Table.setColumn("nuc median", nuc_median, resultsTable);
		Table.setColumn("cyto median", cyto_median, resultsTable);
		if(invertRatio == false) {
			Table.setColumn("n/c ratio mean", ratio_mean, resultsTable);
			Table.setColumn("n/c ratio median", ratio_median, resultsTable);
		}
		else {
			Table.setColumn("c/n ratio mean", ratio_mean, resultsTable);
			Table.setColumn("c/n ratio median", ratio_median, resultsTable);
		}
		Table.setColumn("background value", background_array, resultsTable);
	}
}


//Adds a scalar to all elements of an array
function subtract_scalar_from_array(array, scalar) {
	added_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		added_array[a]=array[a] - scalar;
	}
	return added_array;
}