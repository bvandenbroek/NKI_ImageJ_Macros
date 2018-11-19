
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = "tif") suffix

#@ Integer (label = "Median filter radius", value=2, style="spinner", min=0, max=20) filter_radius
#@ String(choices={"Global Threshold", "Local Threshold"}, style="listBox", value="Local Threshold", label = "Threshold type") threshold_type
#@ String(label = "Global threshold method (if applicable)", value="Li", choices={"Default", "Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"}, style="listbox") global_th_method
#@ String(label = "Local threshold method (if applicable)", value="Niblack", choices={"Bernsen", "Contrast", "Mean", "Median", "MidGrey", "Niblack", "Otsu", "Phansalkar", "Sauvola"}, style="listbox") local_th_method
#@ Integer (label = "Local threshold radius", value=10, style="spinner", min=0) radius
#@ Float (label = "Parameter 1 value", style = "spinner", value=0,5, min=-10000, max=10000) par1
#@ Float (label = "Parameter 2 value", style = "spinner", value=-5, min=-10000, max=10000) par2


#@ Integer (label = "Intensity threshold (to prevent false counts in empty wells)", style = "spinner", min=0, max=65535, value=2000) min_intensity

#@ Integer (label = "Lower diameter limit", style = "spinner", min=0, max=1000, value=10) lower_diameter_limit
#@ Integer (label = "Upper diameter limit", style = "spinner", min=0, max=1000, value=100) upper_diameter_limit
#@ Boolean (label = "Exclude nuclei on edges", value=true) exclude_edges

#@ Boolean (label = "Show images during processing", value=false) show_images
#@ Boolean (label = "Save RGB images with overlay", value=true) save_images
#@ String(choices={"TIF", "PNG"}, style="listBox", value="TIF", label = "as which file type") save_file_type

lower_size_limit = lower_diameter_limit*lower_diameter_limit/4*PI;
upper_size_limit = upper_diameter_limit*upper_diameter_limit/4*PI;

var n=0;
var current_image_nr=0;
var processtime=0;

saveSettings;

run("Conversions...", "scale");
run("Close All");
run("Set Measurements...", "area mean redirect=None decimal=3");
setBatchMode(true);

if(!File.exists(output)) {
	create = getBoolean("The specified output folder "+output+" does not exist. Create?");
	if(create==true) File.makeDirectory(output);		//create the output folder if it doesn't exist
	else exit;
}
if(isOpen("Summary")) close("Summary");

start = getTime();
print("\\Clear");
print("\n");

scanFolder(input);
processFolder(input);

selectWindow("Summary");
saveAs("results",output + File.separator + "Results.txt");
end = getTime();
print("-------------------------------------------------------------------");
print("Finished processing "+n+" images in "+d2s((end-start)/1000,1)+" seconds. ("+d2s((end-start)/1000/n,1)+" seconds per image)");

restoreSettings;


// function to scan folders/subfolders/files to count files with correct suffix
function scanFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			scanFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			n++;
	}
}


// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix)) {
			current_image_nr++;
			showProgress(current_image_nr/n);
			processFile(input, output, list[i]);
		}
	}
}

function processFile(input, output, file) {
	starttime = getTime();
	print("\\Update1:Processing image "+current_image_nr+"/"+n+": " + input + file);
	print("\\Update2:Average speed: "+d2s(current_image_nr/processtime,1)+" images per minute.");
	time_to_run = (n/(current_image_nr/processtime)-processtime);
	if(time_to_run<5) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes.");
	else if(time_to_run<60) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. You'd better get some coffee.");
	else if(time_to_run<480) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes ("+d2s(time_to_run/60,1)+" hours). You'd better go and do something useful.");
	else if(time_to_run<1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. ("+d2s(time_to_run/60,1)+" hours). You'd better come back tomorrow.");
	else if(time_to_run>1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. This is never going to work. Give it up!");

	roiManager("reset");
	open(input + File.separator + file);
	original=getTitle;
	run("Enhance Contrast", "saturated=0.35");
	if(show_images==true) setBatchMode("show");
	run("Properties...", "unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
//	well_id = File.directory;
	well_id = substring(file,lengthOf(file)-lengthOf(suffix)-4,lengthOf(file)-lengthOf(suffix)-1);
//	well_id = substring(well_id,lastIndexOf(well_id,File.separator)+1,lengthOf(well_id));

	run("Duplicate...", "title=["+well_id+"]");
	if(filter_radius>0) run("Median...", "radius="+filter_radius);
	
	resetMinAndMax;

	//Get rid of zero pixel edges (due to e.g. stitching)
	List.setMeasurements();
	mean = List.getValue("Mean");
	changeValues(0,0,mean);
	//Prevent that scaling when converting in empty images causes false measurements
	getMinAndMax(min,max);
	if(max<min_intensity) setMinAndMax(min,min_intensity);

	if(threshold_type == "Global Threshold") {
		setAutoThreshold(global_th_method+" dark");
		run("Convert to Mask");
	} else {
		run("8-bit");
		run("Auto Local Threshold", "method="+local_th_method+" radius="+radius+" parameter_1="+par1+" parameter_2="+par2+" white");
	}
	run("Watershed");
	if(exclude_edges) run("Analyze Particles...", "size="+lower_size_limit+"-"+upper_size_limit+" circularity=0.33-1.00 display exclude summarize add");
	else run("Analyze Particles...", "size="+lower_size_limit+"-"+upper_size_limit+" circularity=0.33-1.00 display summarize add");
	selectWindow(original);
	roiManager("Show all without labels");
	if(roiManager("count")>0) run("From ROI Manager");	//make overlay
	run("Flatten");
//	setBatchMode("show");
	if(save_images==true) saveAs(save_file_type, output + File.separator + well_id);
	run("Close All");
	endtime = getTime();
	processtime = processtime+(endtime-starttime)/60000;
}
