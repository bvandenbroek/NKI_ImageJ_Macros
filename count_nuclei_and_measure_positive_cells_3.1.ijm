/* Macro to calculate the amount and fraction of positively stained cells in a (series of) image(s).
 * Requires a multichannel file with channels staining nuclei and positive nuclei/cells as input.
 * Also works on Deltavision, lif etc. files with multiple series
 *
 *	1.8: Added manual illumination correction for André
 *	2.0: Automatic subtraction of 2^15 (32768) if appropriate (Deltavision data sometimes has this).
 *  
 *  2.1: - fixed a bug regarding calibration removal (in Deltavision files)
 *  	 - changed the automatic positive threshold to minimum. Change if needed.
 * 
 * 	3.0: - added possibility to measure in a band around the nucleus
 * 		 - improved image opening (no Bioformats windowless any more)
 * 		 - saving in results folder
 * 		 - Overlay of original channel when thresholding positive cells
 * 		 - Possibility to select a region for analysis (also drawn in final image)
 *  
 *	Bram van den Broek, Netherlands Cancer Institute, 2015-2017
 *	b.vd.broek@nki.nl
 */

var auto_threshold_nuclei = true;
var threshold_method = "auto local threshold";
var threshold_nuc_bias = 0;		//manual absolute bias of the autothreshold method
var min_nuc_threshold = 300;	//minimum threshold for nucleus recognition
var parameter_1 = -5;			//parameter 1 for Auto Local Threshold
var exclude_edges = false;
var watershed = true;
var pause = false;
var analyze_all = false;
var nucleus_size = 15;
var Min_Nucleus_Size = 4;
var Max_Nucleus_Size = 40;
var ch_nuclei = 1;
var ch_positive = 2;
var auto_threshold_positive = true;
var fixed_threshold = false;
var fixed_threshold_number = 0;
var threshold_positive = 0;
var measure_band = false;			//measure the signal in a band around the nucleus
var band_size = 1;					//band around nucleus (in units)
var save_images = true;
var measure_signal = "Mean";
var overlay = true;
var file_name;
var nr_series;
var format;

var crop_borders = false;
var border = 28;

saveSettings();

print("\\Clear");
roiManager("reset");
run("Clear Results");
run("Set Measurements...", "  mean standard median redirect=None decimal=3");
run("Colors...", "foreground=white background=black selection=cyan");
run("Line Width...", "line=1");

if(nImages>0) run("Close All");
path = File.openDialog("Select a File");

run("Bio-Formats Macro Extensions");
run("Bio-Formats Importer", "open=["+path+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
dir = File.getParent(path)+"\\";
savedir= dir+"\\results\\";
if(!File.exists(savedir)) File.makeDirectory(savedir);
Ext.getFormat(path, format);
file_name = File.getName(path);

extension_length=(lengthOf(file_name)- lastIndexOf(file_name, ".")-1);
extension = substring(file_name, (lengthOf(file_name)-extension_length));
file_list = getFileList(dir); //get filenames of directory
getPixelSize(unit, pw, ph, pd);
getDimensions(width, height, channels, slices, frames);
Ext.setId(dir+file_name);
Ext.getSeriesCount(nr_series);
print(nr_series+" series found");

Stack.setDisplayMode("composite");
run("Blue");


//make a list of images with 'extension' as extension.
j=0;
image_list=newArray(file_list.length);	//Dynamic array size doesn't work on some computers, so first make image_list the maximal size and then trim.
for(i=0; i<file_list.length; i++){
	if (endsWith(file_list[i],extension) && !endsWith(file_list[i],"analyzed."+extension)) {
		image_list[j] = file_list[i];
		j++;
	}
}
image_list = Array.trim(image_list, j);	//Trimming the array of images
nuclei_list = newArray(image_list.length*nr_series);
positive_list = newArray(image_list.length*nr_series);

//print("\\Clear");

//---------CONFIG FILE INITIATIONS
tempdir = getDirectory("temp");
print(tempdir);
config_file = tempdir+"\\count_nuclei_config.txt";
if (File.exists(config_file)) {
	print("config file detected: "+config_file);
	config_string = File.openAsString(config_file);
	config_array = split(config_string,"\n");
	if (config_array.length==22) {
		//print("loading config file...");
		ch_nuclei = 		parseInt(config_array[0]);
		ch_positive = 		parseInt(config_array[1]);
		ch_extra = 			parseInt(config_array[2]);
		measure_signal = 	config_array[3];
		measure_band =		parseInt(config_array[4]);
		band_size = 		parseInt(config_array[5]);
		auto_threshold_positive =parseInt(config_array[6]);
		fixed_threshold=	parseInt(config_array[7]);
		fixed_threshold_number = parseInt(config_array[8]);
		nucleus_size =		parseInt(config_array[9]);
		Min_Nucleus_Size =	parseInt(config_array[10]);
		Max_Nucleus_Size =	parseInt(config_array[11]);
		threshold_method =	config_array[12];
		auto_threshold_nuclei =	parseInt(config_array[13]);
		threshold_nuc_bias =	parseFloat(config_array[14]);
		min_nuc_threshold =	parseInt(config_array[15]);
		exclude_edges =		parseInt(config_array[16]);
		overlay = 			parseInt(config_array[17]);
		watershed =			parseInt(config_array[18]);
		pause =				parseInt(config_array[19]);
		analyze_all =		parseInt(config_array[20]);
		save_images =		parseInt(config_array[21]);
	}
}
threshold_method_array = newArray("auto local threshold","Otsu","Li");
measure_signal_method = newArray("Mean","Median","StdDev");
Dialog.create("Options");
	Dialog.addSlider("Nuclei segmentation channel", 1, channels, ch_nuclei);
	Dialog.addSlider("positive index channel", 1, channels, ch_positive);
	if(channels>2) Dialog.addSlider("Channel to quantify intensity", 1, channels, 3);
	Dialog.setInsets(22, 0, 0);
	Dialog.addChoice("Measure the following signal", measure_signal_method, measure_signal);
	Dialog.addCheckbox("Measure in a band around the nucleus, with size", measure_band);
	Dialog.setInsets(-22, 85, 0);
	Dialog.addNumber("",band_size,1,2,unit);
	Dialog.addCheckbox("Automatic thresholding on positive signal?", auto_threshold_positive);
	Dialog.addCheckbox("Fix threshold at", fixed_threshold);
	Dialog.setInsets(-23, 0, 0);
	Dialog.addNumber("", fixed_threshold_number, 0, 3, "");
	Dialog.addNumber("Estimated nucleus diameter",nucleus_size,0,2,unit);
	Dialog.addSlider("Mininum nucleus diameter ("+unit+")", 1, 100, Min_Nucleus_Size);
	Dialog.addSlider("Maximum nucleus diameter ("+unit+")", 1, 100, Max_Nucleus_Size);
	Dialog.addRadioButtonGroup("Nuclei segmentation threshold method", threshold_method_array, 3, 1, threshold_method);
	Dialog.addCheckbox("Automatic segmentation (Otsu & Li only), with threshold bias", auto_threshold_nuclei);
	Dialog.addNumber("", threshold_nuc_bias, 1, 5, "(-"+pow(2,bitDepth)-1+" to "+pow(2,bitDepth)-1+")");
	Dialog.addNumber("Minimum threshold for nuclei detection", min_nuc_threshold, 0, 4, "(0 to "+pow(2,bitDepth)-1+")");
	Dialog.addCheckbox("Exclude nuclei on edges?", exclude_edges);
	Dialog.addCheckbox("Overlay positive channel when manually thresholding nuclei?", overlay);
	Dialog.addCheckbox("Use watershed to separate touching nuclei?", watershed);
	Dialog.addCheckbox("Inspect results for each image?", pause);
	Dialog.addCheckbox("Analyze all "+image_list.length+" ."+extension+" image files in this directory?", analyze_all)
	Dialog.addCheckbox("Save result images?", save_images);
Dialog.show;
ch_nuclei=Dialog.getNumber();
ch_positive=Dialog.getNumber();
if(channels>2) ch_extra=Dialog.getNumber();
else ch_extra=0;
measure_signal=Dialog.getChoice();
measure_band=Dialog.getCheckbox();
band_size=Dialog.getNumber();
auto_threshold_positive=Dialog.getCheckbox();
fixed_threshold = Dialog.getCheckbox();
fixed_threshold_number = Dialog.getNumber();
nucleus_size=Dialog.getNumber();
Min_Nucleus_Size=Dialog.getNumber();
Max_Nucleus_Size=Dialog.getNumber();
threshold_method=Dialog.getRadioButton();
auto_threshold_nuclei=Dialog.getCheckbox();
threshold_nuc_bias = Dialog.getNumber();
min_nuc_threshold = Dialog.getNumber();
exclude_edges=Dialog.getCheckbox();
overlay=Dialog.getCheckbox();
watershed=Dialog.getCheckbox();
pause=Dialog.getCheckbox();
analyze_all=Dialog.getCheckbox();
save_images=Dialog.getCheckbox();

if(threshold_method=="auto local threshold") auto_threshold_nuclei=true;

//Save settings to config file
save_config_file();

current_image_nr = 0;
//Ext.setId(dir+file_name);
//Ext.getSeriesCount(nr_series);
//print(file_name+": "+nr_series+" series");

setBatchMode(true);

do {
	for(s=0;s<nr_series;s++) {	//loop over all series
		if (analyze_all==true) {
			run("Close All");
			file_name = image_list[current_image_nr];	//retrieve file name from image list
			run("Bio-Formats Importer", "open=["+dir+file_name+"] autoscale color_mode=Composite view=Hyperstack stack_order=XYCZT series_"+s+1);
			run("Blue");
		}
		//Remove senseless "calibration" of Deltavision
		info = getImageInfo();
		index = indexOf(info, "Calibration function");
		if (index!=-1) {
			run("Calibrate...", "function=None unit=[Gray Value]");
			run("Subtract...", "value=32768 stack");
				for(c=1;c<=channels;c++) {
			Stack.setChannel(c);
				resetMinAndMax();
				}
			Stack.setChannel(1);
		}
		
		if(crop_borders==true) {
			makeRectangle(border, border, width-2*border, height-2*border);
			run("Crop");
		}
		setBatchMode("show");
		roiManager("Reset");

		if (auto_threshold_nuclei==true) setBatchMode(true);
		roiManager("reset");
		if(analyze_all==true) {
			//run("Close All");
			//file_name = image_list[current_image_nr];		//retrieve file name from image list
			//Ext.openImagePlus(dir+file_name);			//open file using LOCI Bioformats plugin
		}
		else pause=false;
		image=getTitle();

		Stack.setDisplayMode("composite");
		//setTool("freehand");
		//waitForUser("Make a selection marking the region to analyze and click OK. If no selection is made the whole image is analyzed");
		if(selectionType<0 || selectionType>3) run("Select All");
		roiManager("add");
		run("Select None");

		Stack.setChannel(ch_nuclei);
		run("Duplicate...", "title=nuclei");
		run("32-bit");
		run("Cyan Hot");
		roiManager("Select",0);
		roiManager("Set Color", "cyan");
		//roiManager("Set Line Width", 2);
		run("Add Selection...");	//Create temporary overlay of selected region		
		setBatchMode("show");
		run("Clear Outside");
		run("Make Inverse");

		setThreshold(1, 65535);
	//	run("NaN Background");		//Warning: Subtract background cannot handle NANs
		resetThreshold();
		run("Select None");

	//	roiManager("reset");

		selectWindow(image);
		Stack.setChannel(ch_positive);
		run("Duplicate...", "title=positive");
		run("Grays");
		if(channels>2) {
			selectWindow(image);
			Stack.setChannel(ch_extra);
			run("Duplicate...", "title=third_channel");
			run("Grays");
		}

		nuclei = count_nuclei("nuclei");
		positives = get_positives("positive");

		print(file_name+", series "+s+1+": \t"+nuclei+" nuclei, of which "+positives+" with "+measure_signal+" above threshold ("+threshold_positive+" gray values): "+d2s(positives/nuclei*100,1)+"% positive cells.");

		//Prepare final image
		selectWindow("nuclei");
		run("To ROI Manager");	//Move selected region to the ROI manager
		roiManager("Remove Channel Info");
		roiManager("Remove Slice Info");
		roiManager("Remove Frame Info");
		close("nuclei");
		selectWindow(image);
		run("From ROI Manager");//Make overlay of selected region
		Stack.setChannel(ch_nuclei);
		run("Blue");
		run("Enhance Contrast", "saturated=1");
		Stack.setChannel(ch_positive);
		run("Green");
		run("Enhance Contrast", "saturated=0.5");
		Stack.setDisplayMode("color");
		for(i=1;i<=channels;i++) {		//remove unused channels
			if(i!=ch_nuclei && i!=ch_positive) {
				Stack.setChannel(i);
				run("Delete Slice", "delete=channel");
			}
		}
		Stack.setChannel(2);
		run("Add Slice", "add=channel");
		Stack.setChannel(3);
		run("Red");
		run("Paste");		//Paste positive mask in the new channel
		setMinAndMax(0, 255);
		roiManager("Select",0);
		run("Draw", "stack");
		Stack.setDisplayMode("composite");
		setBatchMode("show");
		if(save_images==true) {
			roiManager("Show None");
			saveAs("tiff", savedir+substring(image,0,lengthOf(image)-1-extension_length)+"_series_"+s+1+"_analyzed");
		}
		//if (analyze_all==false || pause==true) setBatchMode(false);

		roiManager("Show All without labels");
		if(pause==true && current_image_nr<image_list.length) {
			Stack.setDisplayMode("composite");
			setBatchMode("show");
			waitForUser("Inspect result and click OK to continue");
			next=getBoolean("Next image? Click 'No' to re-analyze current image.");
		}
		else next=true;
		if(next==true) {
			nuclei_list[current_image_nr*nr_series + s] = nuclei;
			positive_list[current_image_nr*nr_series + s] = positives;
			if(s+1==nr_series) current_image_nr++;
		}
	} //end of for loop over series

	//combine analyzed series into one stack
	if(save_images==true && nr_series>1) {
		run("Close All");
		images_string = "";
		for(s=0;s<nr_series;s++) {
			open(savedir+substring(image,0,lengthOf(image)-1-extension_length)+"_series_"+s+1+"_analyzed.tif");
			images_string += "image"+s+1+"=["+substring(image,0,lengthOf(image)-1-extension_length)+"_series_"+s+1+"_analyzed.tif] ";
			trash = File.delete(savedir+substring(image,0,lengthOf(image)-1-extension_length)+"_series_"+s+1+"_analyzed.tif");
		}
		run("Concatenate...", "  title=[Concatenated Stacks] keep open "+images_string);
		setBatchMode("show");
		saveAs("tiff", savedir+substring(image,0,lengthOf(image)-1-extension_length)+"_analyzed");
		for(s=0;s<nr_series;s++) {
			close(substring(image,0,lengthOf(image)-1-extension_length)+"_series_"+s+1+"_analyzed.tif");
		}
	}	
} while(analyze_all==true && current_image_nr<image_list.length)


run("Set Measurements...", "  mean display redirect=None decimal=3");
run("Measure");	//just measure something to get a Results window if not open
run("Clear Results");
if(analyze_all==true) results_length = image_list.length;
else results_length = 1;
for(i=0;i<results_length;i++) {
	for(s=0;s<nr_series;s++) {
		if(analyze_all==true) setResult("File name", i*nr_series + s, image_list[i]);
		else setResult("File name", i*nr_series + s, file_name);
		setResult("Series", i*nr_series + s, s+1);
		setResult("nuclei count", i*nr_series + s, nuclei_list[i*nr_series + s]);
		setResult("positive count", i*nr_series + s, positive_list[i*nr_series + s]);
		setResult("% positive cells", i*nr_series + s, d2s(positive_list[i*nr_series + s]/nuclei_list[i*nr_series + s]*100,1));
	}
}
selectWindow("Results");
saveAs("text", savedir+"results.txt");

if(analyze_all==true) {
	results_file_totals = File.open(savedir+"\\results_totals.txt");
	print(results_file_totals, "nr\tFile name\tnuclei count\tpositive count\t% positive cells");
	for(i=0;i<results_length;i++) {
		nuclei_count_this_series = Array.slice(nuclei_list,nr_series*i,nr_series*(i+1));
		Array.getStatistics(nuclei_count_this_series, min, max, mean, stdDev);
		nuclei_count_this_file = round(nr_series*mean);		//easiest way of getting the sum of the array
		positive_count_this_series = Array.slice(positive_list,nr_series*i,nr_series*(i+1));
		Array.getStatistics(positive_count_this_series, min, max, mean, stdDev);
		positive_count_this_file = round(nr_series*mean);
		print(results_file_totals, i+"\t"+image_list[i]+"\t"+nuclei_count_this_file+"\t"+positive_count_this_file+"\t"+100*positive_count_this_file/nuclei_count_this_file);
	}
	File.close(results_file_totals);
}

restoreSettings();







function count_nuclei(image1) {
	selectWindow(image1);
	run("Duplicate...", "title=segmented_nuclei duplicate");
	run("Subtract Background...", "rolling="+(3*nucleus_size/pw)+" sliding");
	run("Median...", "radius="+nucleus_size/pw/10);	//Tweaking the median filter
	if(threshold_method=="auto local threshold") {
		run("8-bit");
		setAutoThreshold("Percentile");
		run("Create Selection");
		List.setMeasurements();
		median = List.getValue("Median");
		std = List.getValue("StdDev");
		resetThreshold();
		run("Select None");
		parameter_1=-(4*std);					//empirical choice
		parameter_1=-minOf(-parameter_1,10);	//minimum value
		parameter_1=-maxOf(-parameter_1,2);		//maximum value
		run("Auto Local Threshold", "method=Mean radius="+(nucleus_size/pw/1.5)+" parameter_1="+parameter_1+" parameter_2=0 white");
		print("\\Update:");	//remove statement from auto local threshold plugin
		getThreshold(min,max);
		resetThreshold();
		setThreshold(min,max);
	}
	else if(auto_threshold_nuclei==true) {
		setAutoThreshold(threshold_method+" dark stack");
		getThreshold(min,max);
		resetThreshold();
		setThreshold(maxOf(min+threshold_nuc_bias,min_nuc_threshold),max);
		//print(min+threshold_nuc_bias, min_nuc_threshold, max);
	}
	else {
		setAutoThreshold(threshold_method+" dark stack");
		getThreshold(min,max);
		resetThreshold();
		setThreshold(min+threshold_nuc_bias,max);
		min_old=min;
		selectWindow("segmented_nuclei");
		setBatchMode("show");
		run("Threshold...");
		selectWindow("Threshold");
		waitForUser("Set threshold for segmentation of nuclei and press OK");
		selectWindow("segmented_nuclei");
		getThreshold(min,max);
		threshold_nuc_bias=d2s(min-min_old,1);
		store_bias = getBoolean("Store threshold bias ("+threshold_nuc_bias+") in config file for later use?");
		if (store_bias==true) save_config_file();
		setThreshold(maxOf(min,min_nuc_threshold),max);
	}
	run("Convert to Mask", "  black");
	run("Fill Holes");
	if(watershed==true) run("Watershed");
	setThreshold(127, 255);
	roiManager("Select",0);	//Analyze particles only in selected region
	run("From ROI Manager");
	roiManager("reset");
	if(exclude_edges==true) run("Analyze Particles...", "size="+PI/4*Min_Nucleus_Size*Min_Nucleus_Size+"-"+PI/4*Max_Nucleus_Size*Max_Nucleus_Size+" circularity=0.20-1.00 show=Nothing exclude add");
	else run("Analyze Particles...", "size="+PI/4*Min_Nucleus_Size*Min_Nucleus_Size+"-"+PI/4*Max_Nucleus_Size*Max_Nucleus_Size+" circularity=0.20-1.00 show=Nothing add");
	close("segmented_nuclei");
	//close("nuclei");	//do not close yet, because it contains the selected region as overlay
	
	return roiManager("count");	//minus one because it contains the selected region as well
}



function get_positives(image1) {
	selectWindow(image1);
	run("Subtract Background...", "rolling="+(3*nucleus_size/pw)+" sliding");
	run("Set Measurements...", "  mean standard median redirect=None decimal=3");
	nr_positives=0;
if(auto_threshold_positive == false && fixed_threshold == false) setBatchMode(false);
	if (crop_borders==true) newImage("mask_positive", "16-bit black", width-2*border, height-2*border, 1);
	else newImage("mask_positive", "16-bit black", width, height, 1);
if(auto_threshold_positive == false && fixed_threshold == false) setBatchMode(true);
	signal = newArray(roiManager("count"));
	selectWindow(image1);
	for(i=0;i<roiManager("count");i++) {	//Measure the signals
		roiManager("Select", i);

		//Measure only in band around the nucleus
		if(measure_band==true) run("Make Band...", "band="+band_size);
		//roiManager("update");

		List.setMeasurements();
		signal[i]=List.getValue(measure_signal);
	}
	for(i=0;i<roiManager("count");i++) {	//Drawing the nuclei and filling them with the measured positive signal
		selectWindow("mask_positive");
		roiManager("Select", i);
		changeValues(0, 65535, signal[i]);
	}
	selectWindow("mask_positive");
	setThreshold(1, pow(2,bitDepth-1));
	run("Create Selection");
	setAutoThreshold("Huang dark");		//THRESHOLD: This could be Otsu, Li or something else
	getThreshold(threshold_positive, upper);
	run("Select None");
if(auto_threshold_positive == false && fixed_threshold == false) setBatchMode(false);
	selectWindow("mask_positive");
	run("Threshold...");
	setThreshold(threshold_positive, 255);
	if(auto_threshold_positive == false && fixed_threshold == false) {
		if(overlay==true) {
			zero_pixel = getPixel(0,0);
			setPixel(0,0,255);
			resetMinAndMax();
			run("Add Image...", "image=positive x=0 y=0 opacity=50 zero");
			setThreshold(threshold_positive, 255);
		}
		waitForUser("estimated threshold: "+threshold_positive+". Change it if necessary (upper slider) and click OK.");
		getThreshold(threshold_positive, upper);
		if(overlay==true) setPixel(0,0,zero_pixel);
	}
if(auto_threshold_positive == false && fixed_threshold == false) setBatchMode(true);
if(fixed_threshold == true) threshold_positive = fixed_threshold_number;
	selectWindow("mask_positive");
	for(i=0;i<roiManager("count");i++) {
		//Filling only the nuclei above threshold
		if(signal[i]>threshold_positive) {
			nr_positives++;
			roiManager("Select", i);
			changeValues(0, 65535, 151);
		}
		else {
			roiManager("Select", i);
			changeValues(0, 65535, 0);
		}
	}
	setMinAndMax(0,255);
	//Drawing all nuclei - only necessary when saving the image
	if(save_images==true) {
		for(i=0;i<roiManager("count");i++) {
			roiManager("Select", i);
			setForegroundColor(255,255,255);
			run("Draw");
		}
	}

	run("Select All");
	run("Copy");
	close("mask_positive");
	close("positive");
//	if (auto_threshold_nuclei==false) setBatchMode(false);
	return nr_positives;
}


function save_config_file() {
	config_file = File.open(tempdir+"\\count_nuclei_config.txt");
	print(config_file, ch_nuclei);			//0
	print(config_file, ch_positive);			//1
	print(config_file, ch_extra);			//2
	print(config_file, measure_signal);		//3
	print(config_file, measure_band);		//4
	print(config_file, band_size);			//5
	print(config_file, auto_threshold_positive); //6
	print(config_file, fixed_threshold);		//7
	print(config_file, fixed_threshold_number);//8
	print(config_file, nucleus_size);			//9
	print(config_file, Min_Nucleus_Size);		//10
	print(config_file, Max_Nucleus_Size);		//11
	print(config_file, threshold_method);		//12
	print(config_file, auto_threshold_nuclei);	//13
	print(config_file, threshold_nuc_bias);		//14
	print(config_file, min_nuc_threshold);		//15
	print(config_file, exclude_edges);		//16
	print(config_file, overlay);			//17
	print(config_file, watershed);			//18
	print(config_file, pause);				//19
	print(config_file, analyze_all);		//20
	print(config_file, save_images);		//21
	File.close(config_file);
}