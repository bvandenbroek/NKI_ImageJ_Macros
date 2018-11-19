 
#@ File (label = "Input intensity file", style = "file") input_intensity_file
#@ File (label = "Input lifetime file", style = "file") input_lifetime_file

#@ File (label = "Output directory", style = "directory") output

#@ Boolean (label = "Segment cells (otherwise use ROIs from roiManager)", value=true) segment_cells
#@ Integer (label = "Nuclei channel", value=2, style="spinner", min=0, max=5) ch_nuc
#@ Integer (label = "Cells channel", value=1, style="spinner", min=0, max=5) ch_cells
//IGNORED #@ Integer (label = "Number of positions", value=100, style="spinner", min=0, max=10000) positions

#@ Integer (label = "Median filter radius before nuclei segmentation", value=1, style="spinner", min=0, max=10) filter_radius
#@ String(label = "Local threshold method", value="Niblack", choices={"Bernsen", "Contrast", "Mean", "Median", "MidGrey", "Niblack", "Otsu", "Phansalkar", "Sauvola"}, style="listbox") method
#@ Integer (label = "Local threshold radius", value=50, style="spinner", min=0) radius
#@ Float (label = "Parameter 1 value", style = "spinner", value=0.2, min=-10000, max=10000) par1
#@ Float (label = "Parameter 2 value", style = "spinner", value=-3, min=-10000, max=10000) par2
#@ Integer (label = "Lower diameter limit", style = "spinner", min=0, max=1000, value=10) lower_diameter_limit
#@ Integer (label = "Upper diameter limit", style = "spinner", min=0, max=1000, value=100) upper_diameter_limit
#@ Float (label = "Minimum circularity", style = "spinner", min=0, max=1, value=0.33) min_circularity

#@ Boolean (label = "Delete rois with low intensity/small cells?", value=true) delete_bad_rois
#@ Float (label = "Minimum cell intensity", style = "spinner", min=0, max=1E9, value=2.0) min_cell_intensity
#@ Float (label = "Minimum cell area (um^2)", style = "spinner", min=0, max=1E9, value=100.0) min_cell_area

#@ Boolean (label = "Exclude nuclei on edges", value=true) exclude_edges

#@ Boolean (label = "Refine cell ROIs by thresholding on intensity (slow!)", value=true) threshold_cells

#@ Boolean (label = "Save images with overay as PNG file?") save_images
#@ Boolean (label = "Display images on screen while analyzing?", value=false) display_images



lower_size_limit = lower_diameter_limit*lower_diameter_limit/4*PI;
upper_size_limit = upper_diameter_limit*upper_diameter_limit/4*PI;

lifetime_string = "FlimMeanDecayTime 1 ch1";
//intensity_string = "FlimMeanPhotonArivalTime 1 ch1"	//This is not lifetime, but (some measure of) intensity; actually the one with the least amount of noise
intensity_string = "";

var median_radius = 1;	//median radius when refining ROIs
var n=0;
var current_image_nr=0;

saveSettings;

run("Bio-Formats Macro Extensions");
run("Conversions...", "scale");
setBackgroundColor(0, 0, 0);
run("Colors...", "foreground=white background=black selection=#007777");
setOption("BlackBackground", true);	//This is the important one

print("\\Clear");
if(nImages>0) run("Close All");
run("Set Measurements...", "area mean median standard min integrated limit redirect=None decimal=3");
if(display_images==false) setBatchMode(true);

if(!File.exists(output)) {
	create = getBoolean("The specified output folder "+output+" does not exist. Create?");
	if(create==true) File.makeDirectory(output);		//create the output folder if it doesn't exist
	else exit;
}
if(isOpen("Summary")) close("Summary");

start = getTime();

processFile(input_intensity_file, input_lifetime_file, output);

//selectWindow("Summary");
//saveAs("results",output + File.separator + "Results.txt");
end = getTime();
print("-------------------------------------------------------------------");
print("Finished processing "+roiManager("count")+" cells in "+d2s((end-start)/1000,1)+" seconds. ("+d2s(roiManager("count")/((end-start)/1000),1)+" cells per second)");

restoreSettings;




function processFile(input_intensity_file, input_lifetime_file, output) {
//	filename = File.getName(path);
//	Ext.setId(path);
//	Ext.getSeriesCount(nr_series);
//	open_FALCON_images(path);

open(input_intensity_file);
rename("Intensity_stack");
getPixelSize(unit, pixelWidth, pixelHeight);

open(input_lifetime_file);
rename("Lifetime_stack");
run("physics black");
run("Multiply...", "value=1E9 stack");
setMinAndMax(1, 4);

//	getDimensions(width, height, channels, slices, frames);

	bits = bitDepth;

	selectWindow("Lifetime_stack");
	getDimensions(width_lifetime, height_lifetime, channels_lifetime, slices, frames);
	run("Properties...", "channels="+channels_lifetime+" slices="+slices+" frames="+frames+" unit="+unit+" pixel_width="+pixelWidth+" pixel_height="+pixelWidth+" voxel_depth=1.0000");

	selectWindow("Intensity_stack");
	getDimensions(width_int, height_int, channels_int, slices, frames);
	
	run("Split Channels");
	nuclei = "C"+ch_nuc+"-Intensity_stack";
	cells = "C"+ch_cells+"-Intensity_stack";
	
	if(segment_cells==true) {
		//segment nuclei using Voronoi
		roiManager("reset");
		selectWindow(nuclei);
		run("Duplicate...", "title=["+nuclei+"_filtered] duplicate");
		if(filter_radius>0) run("Median...", "radius="+filter_radius);
		//	run("32-bit");
		//	run("ROF Denoise", "theta=10");	//alternative filtering
		resetMinAndMax;
		run("8-bit");
		run("Auto Local Threshold", "method="+method+" radius="+radius+" parameter_1="+par1+" parameter_2="+par2+" white stack");
		run("Watershed", "stack");

		if(exclude_edges) run("Analyze Particles...", "size="+lower_size_limit+"-"+upper_size_limit+" circularity="+min_circularity+"-1.00 show=Masks exclude stack");
		else run("Analyze Particles...", "size="+lower_size_limit+"-"+upper_size_limit+" circularity="+min_circularity+"-1.00 show=Masks stack");

		run("Invert", "stack");
		run("Voronoi", "stack");
		setThreshold(0,0.001);	//threshold on all voronoi distances
		run("Analyze Particles...", "size="+lower_size_limit+"-Infinity show=Nothing clear add stack");
	}
	
//Refine cell ROIs based on thresholding on the intensity
	if(threshold_cells==true) {
		selectWindow(cells);
		run("Duplicate...", "title="+cells+"_filtered duplicate");
		run("Median...", "radius="+median_radius+" stack");
		cells_filtered = getTitle();
		
		//create list of ROI coordinates
		var ROI_x = newArray(roiManager("count"));	//containers for selection locations
		var ROI_y = newArray(roiManager("count"));
		for(i=0;i<roiManager("count");i++) {
			roiManager("Select",i);
			Roi.getBounds(x, y, ROI_width, ROI_height);
			//getSelectionBounds(x, y, ROI_width, ROI_height);
			ROI_x[i]=x;
			ROI_y[i]=y;
		}
		roiManager("Show None");
		//run("Clear Results");
		nr_cells = roiManager("count");


	//selectWindow(cells);
	showStatus("Refining "+nr_cells+" ROIs...");
		if(display_images==true) setBatchMode("hide");	//Always do this in Batch Mode
		n=0;
		for(i=0;i<nr_cells;i++) {
		//	if(i%1000==0) {
				showProgress(i/nr_cells);
				showStatus("Refining ROI "+i+"//"+nr_cells+"...");
		//	}
			selectWindow(cells_filtered);
			roiManager("Select", i);
			run("Duplicate...", "title=cell_"+i+1);
			//roiManager("Update");

			run("Clear Outside");	//Not necessary because measurements are taken inside the selection (?)
showStatus("Refining ROI "+i+"//"+nr_cells+"...");
			//setAutoThreshold("MaxEntropy dark");
			setThreshold(maxOf(min_cell_intensity/2,2), 255);	//set the threshold to half the minimum intensity, but at least 2

			List.setMeasurements("limit");
			mean = List.getValue("Mean");
			area = List.getValue("Area");
			if(area > min_cell_area && mean > min_cell_intensity) {
				run("Convert to Mask");
				run("Fill Holes");
				run("Create Selection");
				getSelectionBounds(x_cell, y_cell, ROI_width_cell, ROI_height_cell);
				close("Mask");

				//re-set ROI locations
				selectWindow(cells_filtered);
				run("Restore Selection");
				roiManager("Select", i);
				run("Restore Selection");
				setSelectionLocation(ROI_x[n]+x_cell, ROI_y[n]+y_cell);
				//Roi.move(x_cell, y_cell);
				roiManager("Update");
			}
			else {
				//print("ROI "+i+1+" is invalid. Area="+area+", Mean="+mean);
				roiManager("Select", i);
				getSelectionBounds(x_cell, y_cell, ROI_width_cell, ROI_height_cell);
				Roi.setStrokeColor("red");
				setSelectionLocation(ROI_x[n]+x_cell, ROI_y[n]+y_cell);
				roiManager("update");
				roiManager("delete");
				i--;
				nr_cells--;	//remove this cell from the list
			}
			n++;	
			close("cell_"+i+1);
		}
		if(display_images==true) {
			selectWindow(cells_filtered);
			setBatchMode("show");
			setBatchMode(false);	//Return from Batch Mode
		}
		print(i+" \/ "+n+" detected cells have area > "+min_cell_area+" and intensity > "+min_cell_intensity);
	}
	
	//scale lifetime stack to match the size of the intensity stack
	selectWindow("Lifetime_stack");
//	run("Multiply...", "value=1E9 stack");
	run("Select None");
	scale = width_int/width_lifetime;
	run("Scale...", "x="+scale+" y="+scale+" z=1.0 width="+width_int+" height="+height_int+" depth="+frames+" interpolation=None average process create");
	rename("Lifetime_stack_scaled");
	for(i=1;i<=frames;i++) {
		setSlice(i);
		changeValues(0,0,NaN);
	}
	setBatchMode("show");
	showStatus("Measuring cells...");
	nr_cells = roiManager("count");	//Number of valid cells
	mean_lifetime_ = newArray(nr_cells);
	cell_area_ = newArray(nr_cells);
	cell_intensity_ = newArray(nr_cells);

	//multiply lifetime with intensity (for normalization)
	imageCalculator("Multiply create 32-bit stack", "Lifetime_stack_scaled",cells);
	rename("Lifetime_times_intensity");

	//Set zeroes to NaN
	selectWindow(cells);
	run("32-bit");
	roiManager("show all without labels");
	for(i=1;i<=frames;i++) {
		setSlice(i);
		changeValues(0,0,NaN);	//Set zeroes to NaN, also in intensity stack
	}
	run("Grays");
	setBatchMode("show");
	run("Enhance Contrast", "saturated=5");
	print(nr_cells+" cells detected");
	run("Set Measurements...", "area mean redirect=None decimal=3");

	//measure the intensity and area of each cell
	showStatus("Measuring cell intensity and area...");
	j=0;
	ROIs_to_delete = newArray(nr_cells);
	for(i=0;i<nr_cells;i++) {
		if(i%1000==0) showProgress(i/nr_cells);
		roiManager("select",i);
		List.setMeasurements();
		cell_intensity_[i] = List.getValue("Mean");			//intensity in ROI
		cell_area_[i] = List.getValue("Area");				//Area in ROI
		//Delete ROIs with dim and/or small cells
		if(delete_bad_rois==true) {
			if(cell_area_[i] < min_cell_area || cell_intensity_[i] < min_cell_intensity) {	//then delete the ROI
				roiManager("delete");
				i--;
				nr_cells--;	//remove this cell from the list. Next i will overwrite the data. To do: check for the last one
//				ROIs_to_delete[j]=i;
//				j++;
//				deleteFromArray(cell_intensity_, i);
//				deleteFromArray(cell_area_, i);
			}
		}
	}
	showProgress(0);
	showStatus("Deleting ROIs (this may take a while...");
	ROIs_to_delete = Array.trim(ROIs_to_delete, j);
	if(ROIs_to_delete.length>0) {
		roiManager("Select",ROIs_to_delete);
		roiManager("Delete");
	}
	nr_cells = roiManager("count");
	print("Deleting "+ROIs_to_delete.length+" cells with intensity < "+min_cell_intensity+" and area > "+min_cell_area+" um^2.\n"+nr_cells+" cells remain.");

	//measure the lifetime of each cell
	selectWindow("Lifetime_times_intensity");
	for(i=0;i<nr_cells;i++) {
		roiManager("select",i);
		List.setMeasurements();
		mean_lifetime_[i] = List.getValue("Mean")/cell_intensity_[i];	//Lifetime normalized on intensity
	}
	
	//Array.getStatistics(mean_lifetime_, min, max, mean, stdDev);
	//print("Mean lifetime: "+mean+" ± "+stdDev);	//Doesn't work because of some NaN values
/*
	selectWindow();
	run("32-bit");	//for merging
	run("Merge Channels...", "c1=["+DAPI+"] c2=["+GFP+"_corrected] create");
	roiManager("Show All without Labels");
	Stack.setChannel(2);
	setMinAndMax(10,150);
*/
	run("Clear Results");
	for(i=27590;i<nr_cells;i++) {
		setResult("cell_area", i-27590, cell_area_[i]);
		setResult("cell intensity", i-27590, cell_intensity_[i]);
		setResult("lifetime", i-27590, mean_lifetime_[i]);
	}
	updateResults();
	run("Distribution...", "parameter=lifetime or=256 and=[1 - 4]");
	setBatchMode("show");
	selectWindow("Lifetime_stack_scaled");
	setMinAndMax(1, 4);
	roiManager("show all without labels");
//	if(save_images==true) saveAs("Tif", output + File.separator + file + "_analyzed");
//	saveAs("Results", output + File.separator + file + "_results.txt");
//	run("Close All");
}

restoreSettings();


function open_FALCON_images(input_lif_file) {
	//Open only the metadata
	run("Bio-Formats Importer", "open=["+input_lif_file+"] autoscale color_mode=Default display_metadata rois_import=[ROI manager] view=[Metadata only] stack_order=Default");
	//put metadata into arrays
	metadata = split(getInfo("window.contents"),'\n');
	run("Close");
	keys = newArray(metadata.length);
	values = newArray(metadata.length);
	for(i=0;i<metadata.length;i++) {
		temp = split(metadata[i],'\t');
		keys[i] = temp[0];
		values[i] = temp[1];
	}

/*
	//print metadata
	for(i=0;i<metadata.length;i++) {
		print(i+ ": "+ keys[i] + " = \t" + values[i]);
	}
*/

	//Note: BioFormats messes up the numbers. Series in the metadata start with 'series 0', so everything is shifted!
	//Therefore: the arrays have to be sorted correctly, because the series names are not padded with zeros and are sorted incorrectly by BioFormats. Could cause sequence problems.
	
	//Create a string to feed into Bioformats with the series to open
	series_to_open_string = "";
	n=0;
	for(i=0;i<values.length;i++) {
		if(startsWith(keys[i]," Series") && endsWith(values[i], lifetime_string)) {	//Find the name
			print(i+": "+keys[i]+" = "+values[i]);
			series_to_open_string += "series_"+parseInt(substring(keys[i], 8, indexOf(keys[i],"Name")))+1 + " ";	//+1 because BioFormats starts counting at Series0 (in macros)
			n++;
		}
	}

	for(i=0;i<values.length;i++) {
//		if(endsWith(values[i], intensity_string)) {	//Find the name
		if(startsWith(keys[i]," Series") && !matches(values[i],".*FLIM.*") && !matches(values[i],".*Overview.*")) {	//Find the name
			print(i+": "+keys[i]+" = "+values[i]);
			series_to_open_string += "series_"+parseInt(substring(keys[i], 8, indexOf(keys[i],"Name")))+1 + " ";
		}
	}

	//Open the lifetime and intensity series
	run("Bio-Formats Importer", "open=["+input_lif_file+"] autoscale color_mode=Default concatenate_series view=Hyperstack stack_order=XYCZT "+series_to_open_string);
	rename("Lifetime_stack");
	if(nImages==1) {	//Then lifetime and intensity have the same dimensions and are concatenated 
		run("Make Substack...", "delete slices="+n+1+"-"+2*n);
	}
	else {
		selectImage(1);	//Select the first image in the Window menu (=intensity_stack)
	}
	rename("Intensity_stack");	
	run("Enhance Contrast", "saturated=0.35");
	selectWindow("Lifetime_stack");
	run("Multiply...", "value=1E9 stack");
	run("physics black");
	setMinAndMax(1, 5);
}


//Delete the specified entry, shifting the other positions 
function deleteFromArray(array, position) {
	if (position<lengthOf(array)) {
//		temparray_a = Array.slice(array,0,i-1);
//		temparray_b = Array.slice(array,i+1,array.length);
//		Array.concat(temparray_a,temparray_b);
		Array.rotate(array, -position);
		array = Array.slice(array,1,array.length);
		Array.rotate(array, position);
	}
}