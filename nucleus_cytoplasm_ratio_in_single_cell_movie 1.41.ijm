/*
 * Macro to quantify cytoplasm/nucleus ratio
 * Required input: an image (stack).
 * 
 * Multi-file support: the user should first select a file in the to-be-analyzed directory.
 * and set the correct settings in a dialog. Then the user is asked whether to analyze all the
 * files (with the same extension) in the directory and whether the macro should pause after each image.
 * 
 * Works on time series.
 *  
 * Bram van den Broek, Netherlands Cancer Institute, 06-16-2014
 * 
 * Version 1.1, February 2015
 * - many bugfixes, speed improvements, and more. Just never use an older version again.
 * 
 * Version 1.2, February 2015
 * - fixed mistake: sum slices now really does sum slices
 * - option to not subtract the background
 * - automatic subtraction of 2^15 (32768) (Deltavision data), if appropriate
 * 
 * Version 1.3, February 2015 * 
 * - Option to subtract a fixed value as background
 * 
 */

saveSettings();
run("Conversions...", " ");
setOption("BlackBackground", true);
roiManager("Associate", "true");

var pause=false;
var exclude_edges = false;		//exclude nuclei on edge of image
var cytoplasm_width = 2;		//width of the cytoplasm ring around the nucleus
var spacer = 1;					//width of the spacing between nucleus and cytoplasm
var SNR = 2;					//nuclei-wide SNR cytoplasm intensity threshold
var manual_background = 0;		//offset value to subtract
var clean=false;				//cleanup after processing
var large_dots=false;			//display cytoplasm as large dots in stead of small dots
var current_slice = 1;
var median_radius_nuclei = 1;	//size of median filter in nuclei detection procedure		
var analyze_all = false;		//analyze all images in the directory
var pause_after_file = false;
var auto_threshold_nuclei = true;	//find nuclei based on automatic thresholding
var watershed = true;			//separate nuclei using watershed algorithm

var nuclei_found = true;
var current_image_nr = 0;
var start_time = 0;
var run_time = 0;
var th_nuc_min = 0;				//currently not used
var th_nuc_max = 255;			//currently not used
var outliers_threshold=10;		//threshold for cytoplasm outlier removal before retreiving median and stddev of nucleus
var outlier_radius=2;			//outlier removal for nucleus segmentation
var min_cytoplasm_size = 0;		//minimum cytoplasm size in pixels (?)
var median_background;

var avg_median;				//average of outlier-removed median of all nuclei (ok, a bit dubble but you can never be too sure) 
var avg_stddev;				//average of outlier-removed stddev of all nuclei
var avg_background;			//average of background (everything but nuclei)

//default settings
var ch_nuclei = 1;
var ch_cytoplasm = 2;
var Min_Nucleus_Size = 0;
var Max_Nucleus_Size = 100000;
var calculate_nuc_threshold = true;
var bgsubtr = true;
var verbose = false;
var slices=0;				//default value
nuclei_method_array = newArray("automatically detect the best slice", "manually select slice", "maximum intensity projection");
var nuclei_method = "automatically detect the best slice";
cytoplasm_projection_method_array = newArray("sum slices", "maximum intensity projection", "manually select slice", "Extended depth of field projection");
var cytoplasm_projection_method = "sum slices";


if(nImages>0) run("Close All");
//path = File.openDialog("Select a File");
run("Bio-Formats Windowless Importer");
run("Bio-Formats Macro Extensions");
dir = getDirectory("image");

//Subtract 2^15 (32768) if appropriate
getDimensions(width, height, channels, slices, frames);
run("Set Measurements...", "min redirect=None decimal=3");
List.setMeasurements();
min = List.getValue("Min");
if(min>32768) {
	run("Subtract...", "value=32768 stack");
	for(c=1;c<=channels;c++) {
		Stack.setChannel(c);
		resetMinAndMax();
	}
	Stack.setChannel(1);
}


file_name = getInfo("image.filename");
results_file = dir+"\\"+File.nameWithoutExtension+"_results";
merged_image_file = dir+"\\"+File.nameWithoutExtension+"_analyzed";

extension_length=(lengthOf(file_name)- lengthOf(File.nameWithoutExtension)-1);
extension = substring(file_name, (lengthOf(file_name)-extension_length));
file_list = getFileList(dir); //get filenames of directory


//make a list of images with 'extension' as extension.
j=0;
image_list=newArray(file_list.length);	//Dynamic array size doesn't work on some computers, so first make image_list the maximal size and then trim.
for(i=0; i<file_list.length; i++){
	if (endsWith(file_list[i],extension)) {
		image_list[j] = file_list[i];
		j++;
	}
}
image_list = Array.trim(image_list, j);	//Trimming the array of images
print("\\Clear");
print("Directory contains "+file_list.length+" files, of which "+image_list.length+" "+extension+" files.");


//---------CONFIG FILE INITIATIONS
tempdir = getDirectory("temp");
config_file = tempdir+"\\nucleus_cytoplasm_ratio_macro_config.txt";
if (File.exists(config_file)) {
	config_string = File.openAsString(config_file);
	config_array = split(config_string,"\n");
	if (config_array.length==15 || config_array.length==13) {
		//numbers
		ch_nuclei = parseInt(config_array[0]);
		ch_cytoplasm = parseInt(config_array[1]);
		Min_Nucleus_Size = parseInt(config_array[2]);
		Max_Nucleus_Size = parseFloat(config_array[3]);
		spacer = parseFloat(config_array[4]);
		cytoplasm_width = parseFloat(config_array[5]);
		SNR = parseInt(config_array[6]);
		bgsubtr = parseInt(config_array[7]);
		manual_background = parseInt(config_array[8]);
		auto_threshold_nuclei = parseInt(config_array[9]);
		exclude_edges = parseInt(config_array[10]);
		calculate_nuc_threshold = parseInt(config_array[11]);
		verbose = parseInt(config_array[12]);
		//text (saved and read but not displayed as default in the dialog!)
		if(config_array.length==14) nuclei_method = config_array[13];
		if(config_array.length==14) cytoplasm_projection_method = config_array[14];
	}
}

Stack.getDimensions(width, height, channels, slices, frames);
if (bitDepth()!=24 && channels==1) //showMessage("Only one channel detected.");
if (bitDepth()==24) run("Make Composite");
if (bitDepth()==32) run("16-bit");
Stack.getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pw, ph, pd);
setLocation(0, 0);
resetMinAndMax();
if(ch_nuclei>channels) ch_nuclei=channels;
if(ch_cytoplasm>channels) ch_cytoplasm=channels;


//---------OPEN DIALOG
Dialog.createNonBlocking("Options");
	Dialog.addSlider("nuclei channel nr", 1, channels, ch_nuclei);
	Dialog.addSlider("cytoplasm channel nr", 1, channels, ch_cytoplasm);
	if(slices>1) Dialog.addChoice("Select method for nuclei z-projection", nuclei_method_array, nuclei_method);
	if(slices>1) Dialog.addChoice("Select projection method for cytoplasm z-stack", cytoplasm_projection_method_array, cytoplasm_projection_method);
	Dialog.addSlider("Mininum nucleus size ("+unit+"^2)", 1, 10000, Min_Nucleus_Size);
	Dialog.addSlider("Maximum nucleus size ("+unit+"^2)", 1, 10000, Max_Nucleus_Size);
	Dialog.addNumber("Space between nucleus and cytoplasm", spacer, 2, 4, unit);
	Dialog.addNumber("Width of the cytoplasm ring around the nucleus", cytoplasm_width, 2, 4, unit);
	Dialog.addSlider("Signal-to-Noise Ratio for cytoplasm detection", 0.1, 5.0, SNR);	//nr of nucleus-wide stddevs of peak height
	Dialog.addCheckbox("Automatically subtract background?", bgsubtr);
	Dialog.addNumber("Manually subtract background value (uncheck auto-background)", manual_background, 0, 4, "gray values");
	Dialog.addCheckbox("Automatic threshold for nuclei segmentation", auto_threshold_nuclei);
	Dialog.addCheckbox("Exclude nuclei on edge of image", exclude_edges);
	Dialog.addCheckbox("Calculate nucleus threshold for every frame?", calculate_nuc_threshold);	//otherwise stack histogram is used
	Dialog.addCheckbox("Verbose mode (show intermediate screens during processing)", verbose);
	//Dialog.addHelp(url) For later use
Dialog.show;
	ch_nuclei=Dialog.getNumber();
	ch_cytoplasm=Dialog.getNumber();
	if(slices>1) nuclei_method=Dialog.getChoice();
	if(slices>1) cytoplasm_projection_method=Dialog.getChoice();
	Min_Nucleus_Size = Dialog.getNumber();
	Max_Nucleus_Size = Dialog.getNumber();	
	spacer = Dialog.getNumber();
	cytoplasm_width = Dialog.getNumber();
	SNR = Dialog.getNumber();
	bgsubtr = Dialog.getCheckbox();
	manual_background = Dialog.getNumber();
	auto_threshold_nuclei=Dialog.getCheckbox();
	exclude_edges=Dialog.getCheckbox();
	calculate_nuc_threshold=Dialog.getCheckbox();
	verbose=Dialog.getCheckbox();

//---------SAVE SETTINGS IN CONFIG FILE
save_config_file();

//enquire if all files in the directory should be analyzed if the directory contains more than one files with the same extension.
if(image_list.length>1){
	analyze_all=getBoolean("Analyze all "+image_list.length+" "+extension+" files in this directory with these settings?");
	if(analyze_all==true) pause_after_file=getBoolean("Pause after each file?");
}

start_time=getTime();

//START OF DO...WHILE LOOP FOR ANALYZING ALL IMAGES IN A DIRECTORY
do{

/*if(frames>1) {
	rename("original_stack");
	//run("Reduce Dimensionality...", "channels slices keep");
	rename(file_name);
}
*/
roiManager("reset");
run("Clear Results");

if(verbose==false) setBatchMode(true);

if(analyze_all==true) {
	run("Close All");
	file_name = image_list[current_image_nr];	//retrieve file name from image list
	current_image_nr++;
	Ext.openImagePlus(dir+file_name);		//open file using LOCI Bioformats plugin
	//run("32-bit");
	print("current image: "+file_name);
	extension_index = indexOf(file_name, extension)-1;	//index of file extension
	results_file = dir+"\\"+substring(file_name,0,extension_index)+"_results";	//name of results file
	merged_image_file = dir+"\\"+substring(file_name,0,extension_index)+"_analyzed";//name of analyzed image
	Stack.getDimensions(width, height, channels, slices, frames);
}
if (channels>1 && ch_cytoplasm!=ch_nuclei) {
	run("Split Channels");
	selectWindow("C"+ch_cytoplasm+"-"+file_name);
	rename("cytoplasm_original");
	selectWindow("C"+ch_nuclei+"-"+file_name);
	rename("nuclei_original");
}
else if (channels>1) {	//if nuclei and cytoplasm channel are the same
	Stack.setChannel(ch_cytoplasm);
	run("Reduce Dimensionality...", "slices frames keep");
	rename("cytoplasm_original");
	run("Duplicate...", "title=nuclei_original duplicate");
}
else {
	rename("cytoplasm_original");
	run("Duplicate...", "title=nuclei_original duplicate");
}

//---------METHODS FOR NUCLEI DETECTION
if (slices==1) {};
else if (nuclei_method=="automatically detect the best slice") {	//select slice with largest standard deviation
	run("Set Measurements...", "area mean standard median skewness stack limit redirect=None decimal=3");
	setAutoThreshold("Li dark stack");
	run("Select All");
	roiManager("Add");
	roiManager("Multi Measure");
	resetThreshold();
	roiManager("reset");
	var stddev_array = newArray(slices);
	for(i=0;i<slices;i++) {
		stddev_array[i] = getResult("StdDev1",i);
	}
	Array.getStatistics(stddev_array, min, max);
	index_max = newArray();
	index_max = indexOfArray(stddev_array, max);
	setSlice(index_max[0]+1);	//select slice with largest standard deviation
	print("Slice "+index_max[0]+1+" used for nuclei detection");
	run("Duplicate...", "title=nuclei duplicate slices="+index_max[0]+1+" title=nuclei");
	message("Slice "+index_max[0]+1+" selected for nuclei detection.");
}
else if (nuclei_method=="manually select slice") {
	setBatchMode("show");
	waitForUser("Select slice for analysis of nuclei");
	current_slice = getSliceNumber();
	setSlice(current_slice);
	run("Duplicate...", "title=nuclei duplicate slices="+current_slice+" title=nuclei");
}
else if (nuclei_method=="maximum intensity projection") {
	run("Z Project...", "projection=[Max Intensity] all");
}
rename("nuclei");
run("32-bit");






//---------MAIN PART


segment_nuclei();

if(nuclei_found==true) {
	selectWindow("cytoplasm_original");
	message("cytoplasm before background subtraction");
	if(bgsubtr==false && manual_background!=0) {
		run("32-bit");
		run("Subtract...", "value="+manual_background+" stack");
	}
	if (slices>1) z_project();
	rename("cytoplasm");
	run("32-bit");
//	if(verbose==true) setBatchMode(true);	//enable batch mode anyway for this part

	keep_single_nucleus_per_frame();	//select single nucleus and delete extra ROIs

	var nr_of_frames = roiManager("count");
	var background = newArray(nr_of_frames);
	var nucleus_area = newArray(nr_of_frames);
	var nucleus_mean = newArray(nr_of_frames);
	var nucleus_median = newArray(nr_of_frames);
	var nucleus_stddev = newArray(nr_of_frames);
	var cytoplasm_area = newArray(nr_of_frames);
	var cytoplasm_mean = newArray(nr_of_frames);
	var cytoplasm_median = newArray(nr_of_frames);
	var cytoplasm_stddev = newArray(nr_of_frames);
	
	var cyto_nuc_median_ratio = newArray(nr_of_frames);
	var cyto_nuc_mean_ratio = newArray(nr_of_frames);
	var error_cyto_nuc_ratio = newArray(nr_of_frames);

	if(bgsubtr==true) subtract_background("cytoplasm");	//automatically subtract the background outside the cell in the cytoplasm channel (offset)
	measure_cell();

	plot_results("median of nucleus", nucleus_median, nucleus_stddev, "red");
	plot_results("median of cytoplasm", cytoplasm_median, cytoplasm_stddev, "#008800");
	plot_results("cytoplasm/nucleus ratio", cyto_nuc_median_ratio, error_cyto_nuc_ratio, "black");	//to do: calculate error bars


	cleanup();
	handle_results();
	run("Merge Channels...", "c1=cytoplasm c2=nuclei create");
	Stack.setChannel(2);
	run("Red");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(1);
	run("Green");
	run("Enhance Contrast", "saturated=0.35");
	roiManager("Show All without labels");
	Stack.setFrame(1);
	setBatchMode("show");
	run("Set... ", "zoom=400");
	saveAs("tiff", merged_image_file);
}
if(analyze_all==true && current_image_nr<image_list.length) {
	if(pause_after_file==true) waitForUser("Inspect results and click OK to continue with the next file");
	if(pause_after_file==true) run("Close All");
	nuclei_found = true; //reset for next loop
}
if(analyze_all==true && current_image_nr==image_list.length) {
	run_time=round((getTime()-start_time)/1000);
	showMessage("Finished analyzing "+image_list.length+" files in "+run_time+" seconds!");
}
if(analyze_all==false) {
	run_time=round((getTime()-start_time)/1000);
	//showMessage("Finished in "+run_time+" seconds!");
}

//END OF DO...WHILE LOOP
} while(analyze_all==true && current_image_nr<image_list.length);

restoreSettings;












//---------FUNCTIONS

function change_zeros_to_NaN(image) {
	selectWindow(image);
	getDimensions(image_width, image_height, image_channels, image_slices, image_frames);
	for(i=1;i<=image_frames;i++) {
		if(image_frames>1) Stack.setFrame(i);
		changeValues(0, 0, -10000);
		setThreshold(-9999, 65536);
		run("NaN Background", "slice");
	}
}


function segment_nuclei() {
	run("Duplicate...", "title=segmented_nuclei duplicate");
	run("Properties...", "unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
	selectWindow("segmented_nuclei");
	//setBatchMode("show");
	run("Remove Outliers...", "radius="+2*outlier_radius+" threshold="+outliers_threshold+" which=Bright stack");	//to produce a flatter nucleus for segmentation
	run("Log", "stack");
	run("Median...", "radius="+median_radius_nuclei+" stack");
	run("Duplicate...", "title=nuclei_before_segmentation duplicate");
	resetMinAndMax();
	selectWindow("segmented_nuclei");
	change_zeros_to_NaN("segmented_nuclei");
	resetMinAndMax();
	if(auto_threshold_nuclei==true)	{
		if (calculate_nuc_threshold==true) {
			//setAutoThreshold("Otsu dark");						//no stack histogram used
			run("Convert to Mask", "method=Otsu background=Dark calculate black");	//threshold calculated for each frame
		}
		else {
			//setAutoThreshold("Otsu dark stack");
			run("Convert to Mask", "method=Otsu background=Dark black");		//threshold calculated from first frame
		}
	}
	else {
		setAutoThreshold("Otsu dark stack");
		run("Threshold...");
		selectWindow("Threshold");
		waitForUser("Set threshold for segmentation of nuclei");
		run("Convert to Mask", "stack");
//		getThreshold(th_nuc_min, th_nuc_max);	//Currently not used
	}
	run("Fill Holes", "stack");
	if(watershed==true) run("Watershed", "stack");
	setThreshold(127, 255);
	run("Set Measurements...", "area mean standard centroid stack redirect=None decimal=3");
	if(exclude_edges==true) run("Analyze Particles...", "size=Min_Nucleus_Size-Max_Nucleus_Size circularity=0.20-1.00 show=Nothing display exclude add stack");
	else run("Analyze Particles...", "size=Min_Nucleus_Size-Max_Nucleus_Size circularity=0.2-1.00 show=Nothing display add stack");
//setBatchMode(false);
	if(nr_of_frames==0) {
		nuclei_found=false;
		print("No nuclei found!");
	}
	else print(nr_of_frames+" putative nuclei detected in "+frames+" frames");
}


function z_project() {
	if (cytoplasm_projection_method=="sum slices") {
		run("Z Project...", " projection=[Sum Slices] all");
	}
	else if (cytoplasm_projection_method=="manually select slice") {
		setBatchMode("show");
		waitForUser("Select slice for analysis of cytoplasm");
		current_slice = getSliceNumber();
		setSlice(current_slice);
		run("Duplicate...", "title=cytoplasm duplicate slices="+current_slice);
	}
	else if (cytoplasm_projection_method=="maximum intensity projection") {
		run("Z Project...", " projection=[Max Intensity] all");
	}
	if (cytoplasm_projection_method=="Extended depth of field projection") {
		Extended_Depth_of_Field();
		selectWindow("cytoplasm");
	}
}


function keep_single_nucleus_per_frame() {
	showStatus("Calculating central nucleus for each frame...");	//This is very fast, so you don't see it anyway.
	X = newArray(nResults);	//Array with centroids of nuclei
	Y = newArray(nResults);
	distance_to_center = newArray(nResults);
	X0 = getWidth/2;	//center coordinate of image
	Y0 = getHeight/2;
	
	ROI_indices = newArray(nResults);
	total=0;		//running total number of nuclei
	n=0;			//'ROI to delete' counter
	for(f=1;f<=frames;f++) {
		j=0;		//running number of nuclei per frame
		Array.fill(X, 0);	//reset arrays every frame
		Array.fill(Y, 0);
		Array.fill(distance_to_center, 9999);
		for(i=0;i<nResults;i++) {
			if(getResult("Slice", i)==f) {	//Check if this nucleus belongs to slice (frame) f
				X[j]=getResult("X", i);	//get centroid coordinates from measurements in segment_nuclei()
				Y[j]=getResult("Y", i);
				distance_to_center[j]=sqrt( (X[j]-X0)*(X[j]-X0) + (Y[j]-Y0)*(Y[j]-Y0));	//square of the distance
				//print(f, distance_to_center[j]);
				j++;
				total++;
			}
		}
		rankPos = Array.rankPositions(distance_to_center);
		k=0;	//counter for this frame
		//print("j="+j);
		//print("total="+total);
		//print("position="+rankPos[0]);
		for(i=(total-j);i<total;i++) {	//run this loop only j times on the right ROIs
			if(k!=rankPos[0]) {
				ROI_indices[n]=i;	//ROIs to be deleted
				//print("deleting ROI "+i);
				n++;
			}
		k++;
		}
	}
	ROIs_to_delete = Array.trim(ROI_indices, n);	//trim array because it was created too long 
	print("Deleting "+ROIs_to_delete.length+" ROIs");
	roiManager("select", ROIs_to_delete);
	roiManager("Delete");
}


function subtract_background(image) {
	selectWindow(image);
	run("Clear Results");
	run("Select None");
	run("Set Measurements...", "mean median limit redirect=None decimal=3");
	for(i=0;i<nr_of_frames;i++) {
		Stack.setFrame(i+1);
		setAutoThreshold("Triangle");
		List.setMeasurements("limit");
		background[i] = List.getValue("Median");
	}
	resetThreshold();
	run("Select None");
	Array.print(background);
	Array.sort(background);
	median_background = background[floor(lengthOf(background)/2)];	//take the median of median_background
	print("median background value (outside cell): "+median_background);
	run("Subtract...", "value="+median_background+" stack");
}


function measure_cell() {
	selectWindow("cytoplasm");
	run("Clear Results");
	run("Set Measurements...", "area mean standard median redirect=None decimal=3");
	for(i=0;i<nr_of_frames;i++) {
		roiManager("select",i);
		run("Enlarge...", "enlarge=-"+spacer/2);	//shrink the nucleus
		roiManager("Set Color", "yellow");
		roiManager("Update");
		roiManager("Rename", "nucleus_"+i);

		//measure nucleus
		List.setMeasurements();
		nucleus_area[i] = List.getValue("Area");
		nucleus_mean[i] = List.getValue("Mean");
		nucleus_median[i] = List.getValue("Median");
		nucleus_stddev[i] = List.getValue("StdDev");

		//measure cytoplasm
		run("Enlarge...", "enlarge="+spacer);	//grow to avoid nucleus edge
		run("Make Band...", "band="+cytoplasm_width);
		roiManager("Add");		//add cytoplasm band to ROI manager
		roiManager("select",nr_of_frames+i);
		roiManager("Set Color", "cyan");
		roiManager("Rename", "cytoplasm_"+i)

		List.setMeasurements();
		cytoplasm_area[i] = List.getValue("Area");
		cytoplasm_mean[i] = List.getValue("Mean");
		cytoplasm_median[i] = List.getValue("Median");
		cytoplasm_stddev[i] = List.getValue("StdDev");

		cyto_nuc_median_ratio[i] = cytoplasm_median[i]/nucleus_median[i];
		cyto_nuc_mean_ratio[i]   = cytoplasm_mean[i]/nucleus_mean[i];
		error_cyto_nuc_ratio[i]  = sqrt(pow(nucleus_stddev[i]/(nucleus_area[i]/(pw*pw))/nucleus_median[i],2)+pow(cytoplasm_stddev[i]/(cytoplasm_area[i]/(pw*pw))/cytoplasm_median[i],2));
	}
}




function plot_results(name, array, errors, color) {
	//TO DO: get number of arguments and display error bars if errors are passed.
	Array.getStatistics(array, min, max);
	Plot.create(name, "frames", name, array);
	Plot.setLimits(0, frames-1, 0, max);
	Plot.setColor(color);
	Plot.add("boxes", array);
	Plot.add("error bars", errors);
	Plot.show();
	setBatchMode("show");
}


function cleanup() {
	if(clean==true) {
		for(i=1;i<=channels;i++) {
			if(isOpen("C"+i+"-"+file_name)) {
				selectWindow("C"+i+"-"+file_name);
				run("Close");	
			}
		}
		selectWindow("nuclei");
		run("Close");
		selectWindow("segmented_nuclei");
		run("Close");
		selectWindow("cytoplasm");
		run("Close");
	}
}

function handle_results() {
	//Put measured data in Results table
	run("Clear Results");
	for(i=0;i<nr_of_frames;i++) {
		setResult("nuc_area ("+unit+"^2)", i, nucleus_area[i]);
		setResult("nuc_mean", i, nucleus_mean[i]);
		setResult("nuc_median", i, nucleus_median[i]);
		setResult("nuc_stddev", i, nucleus_stddev[i]);
		
		setResult("cyt_area ("+unit+"^2)", i, cytoplasm_area[i]);
		setResult("cyt_mean", i, cytoplasm_mean[i]);
		setResult("cyt_median", i, cytoplasm_median[i]);		
		setResult("cyt_stddev", i, cytoplasm_stddev[i]);		
		
		setResult("cyt/nuc_ratio_mean", i, cyto_nuc_mean_ratio[i]);
		setResult("cyt/nuc_ratio_median", i, cyto_nuc_median_ratio[i]);
		setResult("error_cyt/nuc_ratio", i, error_cyto_nuc_ratio[i]);
	}
	updateResults();
	//Save results and merged image
	selectWindow("Results");
	saveAs("text", results_file);
}


function save_config_file() {
	config_file = File.open(tempdir+"\\nucleus_cytoplasm_ratio_macro_config.txt");
	//numbers
	print(config_file, ch_nuclei);
	print(config_file, ch_cytoplasm);
	print(config_file, Min_Nucleus_Size);
	print(config_file, Max_Nucleus_Size);
	print(config_file, spacer);
	print(config_file, cytoplasm_width);
	print(config_file, SNR);
	print(config_file, bgsubtr);
	print(config_file, manual_background);
	print(config_file, auto_threshold_nuclei);
	print(config_file, exclude_edges);
	print(config_file, calculate_nuc_threshold);
	print(config_file, verbose);
	//text
	if(slices>1) print(config_file, nuclei_method);
	if(slices>1) print(config_file, cytoplasm_projection_method);
	File.close(config_file);	
}	
	
function indexOfArray(array, value) {
	count=0;
	for (a=0; a<lengthOf(array); a++) {
		if (d2s(array[a],3)==d2s(value,3)) {
			count++;
		}
	}
	if (count>0) {
		indices=newArray(count);
		count=0;
		for (a=0; a<lengthOf(array); a++) {
			if (array[a]==value) {
				indices[count]=a;
				count++;
			}
		}
		return indices;
	}
}


function message(txt) {
if (pause==true){
	waitForUser(txt);
	}
}


function Extended_Depth_of_Field() {
	if(verbose==true) setBatchMode(true);	//always use BatchMode here
	radius=2;
	//Get start image properties
	w=getWidth();
	h=getHeight();
	stack_before_focussing = getTitle();

	for(i=1;i<=frames;i++) {
		//print("focusing frame "+i);
		selectWindow(stack_before_focussing);
		Stack.setFrame(i);
		run("Reduce Dimensionality...", "slices keep");
		rename("cytoplasm_frame_"+i);
		d=nSlices();
		source=getImageID();
		origtitle=getTitle();
		rename("tempnameforprocessing");
		sourcetitle=getTitle();

		//Generate edge-detected image for detecting focus
		run("Duplicate...", "title=["+sourcetitle+"_Heightmap] duplicate range=1-"+d);
		heightmap=getImageID();
		heightmaptitle=getTitle();
		run("Find Edges", "stack");
		run("Maximum...", "radius="+radius+" stack");
		//Alter edge detected image to desired structure
		run("32-bit");
		for (x=0; x<w; x++) {
			showStatus("Creating focused image from stack...");
			showProgress(x/w);
			for (y=0; y<h; y++) {
				slice=0;
				max=0;
				for (z=0; z<d; z++) {
					setZCoordinate(z);
					v=getPixel(x,y);
					if (v>=max) {
						max=v;
						slice=z;
					}
				}
				for (z=0; z<d; z++) {
					setZCoordinate(z);
					if (z==slice) {
						setPixel(x,y,1);
					} else {
						setPixel(x,y,0);
					}
				}
			}
		}
		run("Gaussian Blur...", "sigma="+radius+" stack");

		//Generation of the final image

		//Multiply modified edge detect (the depth map) with the source image
		run("Image Calculator...", "image1="+sourcetitle+" operation=Multiply image2="+heightmaptitle+" create 32-bit stack");multiplication=getImageID();
		//Z project the multiplication result
		run("Z Project...", "start=1 stop="+d+" projection=[Sum Slices]");
		//Some tidying

		rename(origtitle+"_focused");
		selectImage(heightmap);
		close();
		selectImage(multiplication);
		close();
		selectImage(source);
		run("Close");
		//rename(origtitle);

		if(i>1) run("Concatenate...", "stack1=cytoplasm_frame_1_focused stack2=cytoplasm_frame_"+i+"_focused title=cytoplasm_frame_1_focused");	//concatenate focused frames
	}
	rename("cytoplasm");
	if(verbose==true) {
		setBatchMode("show");
		setBatchMode(false);
	}
}