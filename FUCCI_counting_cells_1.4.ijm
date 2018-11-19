#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File extension", value = ".tif") suffix

#@ Integer (label = "Green channel", min = 1, max = 3, style = "scroll bar") ch_green
#@ Integer (label = "Red channel", min = 1, max = 3, style = "scroll bar") ch_red
#@ Boolean (label = "Count the nuclei using a third channel?") use_ch_nuc
#@ Integer (label = "nuclei channel", min = 1, max = 3, style = "scroll bar") ch_nuc

#@ Integer (label = "Pixels to cut off the edge", style = "spinner", value=10) edge
#@ Integer (label = "Minimum nucleus diameter (um) (also used in background subtraction)", style = "spinner", min=0, max=20, value=8) min_diameter

#@ String(label = "Filter choice", choices={"Median filter", "Bilateral filter"}, style="radioButtonHorizontal") filter_choice
#@ Float (label = "Filter radius (pixels)", style = "spinner", min=0, max=10, value=5) filter_radius
#@ Integer (label = "Radius for local thresholding (pixels)", style = "spinner", min=0, max=200, value=50) threshold_radius

#@ Float (label = "Threshold strictness Green", style = "spinner", min=0, max=100, value=10) parameter_1_green
#@ Float (label = "Threshold strictness Red", style = "spinner", min=0, max=100, value=10) parameter_1_red
#@ Float (label = "Threshold strictness Nuclei", style = "spinner", min=0, max=100, value=10) parameter_1_nuc

#@ Boolean (label = "export an RGB image stack (slow)?") create_RGB
#@ Integer (label = "Brightness scale factor for RGB stack green (lower=brighter))", style = "slider", min=1, max=40, value=10) brightness_factor_green
#@ Integer (label = "Brightness scale factor for RGB stack red (lower=brighter))", style = "slider", min=1, max=40, value=10) brightness_factor_red
#@ Integer (label = "Brightness scale factor for RGB stack nuclei (lower=brighter))", style = "slider", min=1, max=40, value=10) brightness_factor_nuc


// Version 1.4:
// - Added possibility to use a third channel for nuclei detection, in stead of the green and red itself
// - Display now also scales before converting to 8-bit when choosing median filtering 


var watershed = true;
var invert_mask = true;

//var edge=10;			//pixels to crop off the edge
//var min_size=80;		//minimum size of the cell nuclei
var max_size=99999;		//maximum size of the cell nuclei
//var filter_radius = 5;	//Median filter before thresholding
//var threshold_radius = 50;		//Radius for local threshold
//var parameter_1 = -15;	//parameter for local threshold. More negative means more strict.
//var brightness_factor = 5;	//Maximum level for displaying and saving the RGB image stack 

min_size = pow(min_diameter/2,2)*PI;

saveSettings();

roiManager("Associate", "true");
run("Colors...", "foreground=white background=black selection=cyan");
setForegroundColor(255,255,255);
setBackgroundColor(0,0,0);
setOption("BlackBackground", true);
//run("Conversions...", " ");

starttime=getTime();

processFolder(input);

endtime=getTime();
showMessage("Finished processing in "+(endtime-starttime)/1000+" seconds!");

restoreSettings();



// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	
	setBatchMode(true);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
//			processFolder(input + File.separator + list[i]);	//Do not process subdirectories
			;
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
	setBatchMode(false);
}

function processFile(input, output, file) {
	if(nImages>0) run("Close All");
	print("\\Clear");
	roiManager("Reset");
	print("Processing: " + input + File.separator + file);
	//open(input + File.separator + file);
	run("Bio-Formats Importer", "open=["+input + File.separator + file+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	name = File.nameWithoutExtension;
	getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pixelWidth, pixelHeight);

//	bilateral_radius = 4/pixelWidth;	//Filter radius set at 4 um.
//	if(filter_choice == "Bilateral filter") {
//		if(unit!="micron") showMessage("Warning: units are in "+unit+", where microns is expexted. This may affect the results. Bilateral filter size is now set at "+bilateral_radius+" pixels.");
//	}
	setBatchMode("show");

	run("Set Measurements...", "mean median limit redirect=None decimal=3");
	Stack.setDisplayMode("composite");
	for(c=1;c<=channels;c++) {
		Stack.setChannel(c);
		run("Grays");
	}
	Stack.setChannel(ch_green);
	run("Green");
	resetMinAndMax();
	Stack.setChannel(ch_red);
	run("Red");
	resetMinAndMax();
	run("Set Measurements...", "  redirect=None decimal=3");

	original = getTitle;
	makeRectangle(edge, edge, width-(2*edge), height-(2*edge));
	run("Duplicate...", "title="+original+"_cropped duplicate");
	cropped = getTitle;

	run("32-bit");	//to preserve noise after background subtraction
	bg_radius_nuc = min_diameter/pixelWidth/2;	//radius for background subtraction of the nucleus channel
	run("Subtract Background...", "rolling="+bg_radius_nuc+" sliding stack");	//subtract changing background value over time (for all channels!)

	run("Duplicate...", "title=dup_"+original+" duplicate");	//duplicate the cropped image stack for processing
	image = getTitle;

	for(c=1;c<=channels;c++) {	//again necessary? Maybe because of conversion to 32-bit
		Stack.setChannel(c);
		run("Grays");
	}
	Stack.setChannel(ch_green);
	run("Green");
	resetMinAndMax();
	Stack.setChannel(ch_red);
	run("Red");
	resetMinAndMax();
	run("Set Measurements...", "  redirect=None decimal=3");

	run("Split Channels");
	if(channels>2) {	//Close the unused channels
		for(c=1;c<channels;c++) {
			if(c!=ch_green && c!=ch_red && use_ch_nuc==false) {
				selectWindow("C"+c+"-"+image);
				run("Close");
			}
			else if (c!=ch_green && c!=ch_red && c!=ch_nuc) {
				selectWindow("C"+c+"-"+image);
				run("Close");
			}
		}
	}


	//GREEN
	selectWindow("C"+ch_green+"-"+image);
	rename("green");
	background_green = get_background();
	signal_green = get_signal();
	stddev_green = get_stddev();

	if(filter_choice == "Median filter") {
		run("Median...", "radius="+filter_radius+" stack");
		setMinAndMax(background_green-stddev_green,signal_green*2);	//scale before converting to 8-bit
		run("8-bit");
	}
	else {
		setMinAndMax(background_green-stddev_green,signal_green*2);	//scale before converting to 8-bit
		run("8-bit");
		stddev_green_8bit = get_stddev();
		run("Bilateral Filter", "spatial="+filter_radius+" range="+2*stddev_green_8bit);
		run("Median...", "radius=1 stack");	//1-pixel median filter to smooth edges
	}
setBatchMode("show");
	run("Auto Local Threshold", "method=Mean radius="+threshold_radius+" parameter_1=-"+parameter_1_green+" parameter_2=0 white stack");
	run("Analyze Particles...", "size="+min_size+"-"+max_size+" show=Masks include stack");
	if(watershed==true) {
		if(invert_mask==true) run("Invert", "stack");
		run("Watershed", "stack");
	}
	rename("mask_green");
	mask_green = getImageID;
	nr_rois = roiManager("Count");
	run("Analyze Particles...", "size="+min_size+"-"+max_size+" show=Nothing add include summarize stack");
	for(i=nr_rois;i<roiManager("Count");i++) {
		roiManager("Select",i);
		Roi.setStrokeColor(0, 255, 0);
	}
	//setBatchMode("show");
	selectWindow("green");
	run("Close");





	//RED
	selectWindow("C"+ch_red+"-"+image);
	rename("red");
	background_red = get_background();
	signal_red = get_signal();
	stddev_red = get_stddev();

	if(filter_choice == "Median filter") {
		run("Median...", "radius="+filter_radius+" stack");
		setMinAndMax(background_red-stddev_red,signal_red*2);	//scale before converting to 8-bit
		run("8-bit");
	}
	else {
		setMinAndMax(background_red-stddev_red,signal_red*2);	//scale before converting to 8-bit
		run("8-bit");
		stddev_red_8bit = get_stddev();
		run("Bilateral Filter", "spatial="+filter_radius+" range="+2*stddev_red_8bit);
		run("Median...", "radius="+filter_radius/2.5+" stack");	//additional median filter to smooth edges
	}
setBatchMode("show");

	run("Auto Local Threshold", "method=Mean radius="+threshold_radius+" parameter_1=-"+parameter_1_red+" parameter_2=0 white stack");
	run("Analyze Particles...", "size="+min_size+"-"+max_size+" show=Masks include stack");
	if(watershed==true) {
		if(invert_mask==true) run("Invert", "stack");
		run("Watershed", "stack");
	}
	rename("mask_red");
	mask_red = getImageID;
	nr_rois = roiManager("Count");
	run("Analyze Particles...", "size="+min_size+"-"+max_size+" show=Nothing add include summarize stack");
	for(i=nr_rois;i<roiManager("Count");i++) {
		roiManager("Select",i);
		Roi.setStrokeColor(255, 0, 0);
	}
	//setBatchMode("show");
	selectWindow("red");
	run("Close");




	//YELLOW
	selectWindow("mask_green");
	run("Select None");
	if(invert_mask==true) run("Invert", "stack");
	selectWindow("mask_red");
	run("Select None");
	if(invert_mask==true) run("Invert", "stack");
	imageCalculator("AND create stack", "mask_red","mask_green");
	rename("mask_yellow");
	mask_yellow = getImageID;
	if(invert_mask==true) run("Invert", "stack");
	nr_rois = roiManager("Count");
	run("Analyze Particles...", "size="+min_size+"-"+max_size+" show=Nothing add include summarize stack");
	for(i=nr_rois;i<roiManager("Count");i++) {
		roiManager("Select",i);
		Roi.setStrokeColor(255, 255, 0);
	}
	run("Invert", "stack");
	//setBatchMode("show");


	//NUCLEI
	if (use_ch_nuc==true) {
		selectWindow("C"+ch_nuc+"-"+image);
		rename("nuclei");
		//run("Median...", "radius=2 stack");	//pre-filtering, because it's noisy
		background_nuc = get_background();
		signal_nuc = get_signal();
		stddev_nuc = get_stddev();

		if(filter_choice == "Median filter") {
			run("Median...", "radius="+filter_radius+" stack");
			setMinAndMax(background_nuc-stddev_nuc,signal_nuc*2);	//scale before converting to 8-bit
			run("8-bit");
		}
		else {
			setMinAndMax(background_nuc-stddev_nuc,signal_nuc*2);	//scale before converting to 8-bit
			run("8-bit");
			stddev_nuc_8bit = get_stddev();
			run("Bilateral Filter", "spatial="+filter_radius+" range="+2*stddev_nuc_8bit);
			run("Median...", "radius=1 stack");	//1-pixel median filter to smooth edges
		}
	setBatchMode("show");
		run("Auto Local Threshold", "method=Mean radius="+threshold_radius+" parameter_1=-"+parameter_1_nuc+" parameter_2=0 white stack");
		run("Analyze Particles...", "size="+min_size+"-"+max_size+" show=Masks include stack");
		if(watershed==true) {
			if(invert_mask==true) run("Invert", "stack");
			run("Watershed", "stack");
		}
		rename("mask_nuc");
		mask_nuc = getImageID;
		nr_rois = roiManager("Count");
		run("Analyze Particles...", "size="+min_size+"-"+max_size+" show=Nothing add include summarize stack");
		for(i=nr_rois;i<roiManager("Count");i++) {
			roiManager("Select",i);
			Roi.setStrokeColor(128, 128, 128);
		}
		selectWindow("nuclei");
		run("Close");
	}



	count_green = newArray(frames);
	count_red = newArray(frames);
	count_yellow = newArray(frames);
	if (use_ch_nuc==true) count_nuc = newArray(frames);

	run("Clear Results");
	if(isOpen("Results")) {
		selectWindow("Results");
		run("Close");
	}
	selectWindow("Summary of mask_green");
	rename("Results");
	for(i=0;i<frames;i++) count_green[i] = getResult("Count",i);
	run("Close");
	
	selectWindow("Summary of mask_red");
	rename("Results");
	for(i=0;i<frames;i++) count_red[i] = getResult("Count",i);
	run("Close");

	selectWindow("Summary of mask_yellow");
	rename("Results");
	for(i=0;i<frames;i++) count_yellow[i] = getResult("Count",i);
	run("Close");

	if (use_ch_nuc==true) {
		selectWindow("Summary of mask_nuc");
		rename("Results");
		for(i=0;i<frames;i++) count_nuc[i] = getResult("Count",i);
		run("Close");
	}
		
	//For some reason one of these windows gets renamed to "Results" as well. This is to revert that change.
	selectImage(mask_green);
	rename("mask_green");
	selectImage(mask_red);
	rename("mask_red");
	selectImage(mask_yellow);
	rename("mask_yellow");
	if (use_ch_nuc==true) {
		selectImage(mask_nuc);
		rename("mask_nuc");
	}

	for(i=0;i<frames;i++) {
		//subtract yellow counts to prevent double counting of cells
		count_green[i] = count_green[i] - count_yellow[i];
		count_red[i]   = count_red[i] - count_yellow[i];

		setResult("Green", i, count_green[i]);
		setResult("Red", i, count_red[i]);
		setResult("Yellow", i, count_yellow[i]);
		setResult("G+R+Y)", i, count_green[i]+count_red[i]+count_yellow[i]);
		if (use_ch_nuc==true) setResult("Total", i, count_nuc[i]);

		setResult("Red/Green", i, count_red[i]/count_green[i]);
		setResult("Red/Yellow", i, count_red[i]/count_yellow[i]);
		
		if (use_ch_nuc==true){
			setResult("Green/Total", i, count_green[i]/count_nuc[i]);
			setResult("Red/Total", i, count_red[i]/count_nuc[i]);
			setResult("Yellow/Total", i, count_yellow[i]/count_nuc[i]);
			setResult("(G+R+Y)/Total", i, (count_green[i]+count_red[i]+count_yellow[i])/count_nuc[i]);
		}
	}
	updateResults;

	selectWindow("mask_green");
	run("Close");
	selectWindow("mask_red");
	run("Close");
	selectWindow("mask_yellow");
	run("Close");
	if (use_ch_nuc==true) {
		selectWindow("mask_nuc");
		run("Close");	
	}
	close("mask_green");
	close("mask_red");
	close("mask_yellow");
	if (use_ch_nuc==true) close("mask_nuc");

	if(create_RGB==true) {
		selectWindow(cropped);
		Stack.setChannel(ch_green);
		run("Green");
		setMinAndMax(background_green,brightness_factor_green*signal_green);
		Stack.setChannel(ch_red);
		run("Red");
		setMinAndMax(background_red,brightness_factor_red*signal_red);
		Stack.setChannel(ch_nuc);
		run("Grays");
		setMinAndMax(background_nuc,brightness_factor_nuc*signal_nuc);
		run("RGB Color", "frames");
		rename(cropped+"_RGB");
		setBatchMode("show");

		selectWindow(cropped+"_RGB");
		roiManager("Show all Without labels");
		if(roiManager("count")>0) {
			showStatus("Creating Overlay...");
			run("From ROI Manager");
			run("Flatten", "stack");
		}
		saveAs("tif",output + File.separator + name + "_processed");
	}
	selectWindow("Results");
	saveAs("Results",output + File.separator + name + "_results.xls");


setBatchMode("show");
}




function get_background() {
	run("Set Measurements...", "area mean standard median limit redirect=None decimal=3");
	List.setMeasurements("limit");
	background = List.getValue("Median");
	return background;	
}

function get_signal() {
	run("Set Measurements...", "area mean standard median limit redirect=None decimal=3");
	setAutoThreshold("Triangle dark stack");
	List.setMeasurements("limit");
	signal = List.getValue("Median");
	resetThreshold;
	return signal;	
}

function get_stddev() {
	run("Set Measurements...", "area mean standard median limit redirect=None decimal=3");
	setAutoThreshold("Triangle dark stack");
	List.setMeasurements("limit");
	stddev = List.getValue("StdDev");
	//print(stddev);
	resetThreshold;
	return stddev;	
}
