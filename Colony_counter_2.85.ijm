
var dpi = 800;
var automatic_detection = true;	//Automatic well detection
//var well_diameter = 25000;	//size of the well in micrometers

var min_diameter = 10;		//minimum diameter of a colony (in micrometers)
var min_circ = 0.15;			//minimum circularity of the detected objects (after watershedding if activated)
var max_diameter = 10000;

var exclude = false;				//Exclude colonies on the edge of the well
var watershed = true;			//Watershed detected objects
var analyze_all = false;		//analyze all files in the folder
var manual_adjustment = true;	//Manually adjust the selected circle in the well
var save_images = false;
var pause = false;				//Pause after each image
var	scale = 0.25;				//Scaling for well detection
var	size_outliers = 10;			//For outlier removal in well detection procedure
var	peak_radius = 5;			//radius of cross correlation peak in pixels
var	background_correct = true;	//Correct for background (if this yields bright edges)
var	timeout = 20000;			//timeout for searching the correct number of wells
var rolling_ball_radius= 500;	//Background subtraction in ROIs, just before auto local threshold
var radius = 5000;				//radius of the auto_local_threshold (in units)

var smooth_factor = 5;			//smoothing is required for large fluffy colonies
 
var parameter1 = 10;			//threshold level (higher means more strict)
var	error = false;
var roi_color = "red";
var template = "template_normalized";
var dim = 1;				//start value of the dimension;
var RGB = false;			//Is set to true in case of RGB images

var nr_wells = 6;			//Number of wells
var nr_wells_x=0;
var nr_wells_y=0;
var	well_diameter = 34000;
var exclude_edge_wells = true;
var realtime_colonies = true;

var make_template_selection=false;
var template_path = "D:\\DATA\\Lisa Koob\\template_6well_Lisa.tif";
var pixel_size=dpi/25.4;		//pixel size in micrometers

saveSettings();
//To do: first get color, then put back at the end
//getInfo("selection.color"); //doesn't work
//run("Colors...", "foreground=white background=black selection=red");

setBackgroundColor(0,0,0);
roiManager("Reset");
roiManager("Show None");
run("Clear Results");
print("\\Clear");
run("Set Measurements...", "area mean standard min centroid redirect=None decimal=3");

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
//print(nr_series+" series found");


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

//#@ File (label = "Input file", style = "file") input
Dialog.create("Options");
	Dialog.addNumber("Resolution of the scan",dpi,0,5, "");
	Dialog.addNumber("Number of wells",nr_wells,0,5, "");
	Dialog.addNumber("Well diameter",well_diameter,0,5,unit);
	Dialog.addCheckbox("Automatic well detection?",automatic_detection);
	Dialog.addNumber("smooth factor (scales with colony size)",smooth_factor,0,4,"");
	Dialog.addNumber("radius for background subtraction",rolling_ball_radius,0,4,unit);
	Dialog.addNumber("radius of the local_thresholding",radius,0,4,unit);
	Dialog.addNumber("Threshold offset (higher means more strict)", parameter1,0,2,"");
	Dialog.addNumber("Minimum colony diameter",min_diameter,0,5,unit);
	Dialog.addNumber("Maximum colony diameter",max_diameter,0,5,unit);
	Dialog.addSlider("Minimum circularity", 0, 1, min_circ);
	Dialog.addCheckbox("Use watershed to separate touching colonies?", watershed);
	Dialog.addCheckbox("Exclude edge colonies?", exclude);
//	Dialog.addCheckbox("Create a new well template in the image?", make_template_selection);
	items = newArray("Use the pre-defined template", "Create a new template in the image");
	Dialog.addRadioButtonGroup("Template options", items, 2, 1, "Use the pre-defined template");
	Dialog.addCheckbox("Save images?", save_images);
	Dialog.addCheckbox("Analyze all "+image_list.length+" ."+extension+" image files in this directory?", analyze_all);
	Dialog.addCheckbox("Pause after each image", pause);

Dialog.show;
dpi= Dialog.getNumber;

nr_wells = Dialog.getNumber();
well_diameter = Dialog.getNumber();
automatic_detection = Dialog.getCheckbox();
smooth_factor = Dialog.getNumber();
rolling_ball_radius = Dialog.getNumber();
radius = Dialog.getNumber();
parameter1 = Dialog.getNumber();
min_diameter = Dialog.getNumber();
max_diameter = Dialog.getNumber();
min_circ = Dialog.getNumber();
watershed = Dialog.getCheckbox();
exclude = Dialog.getCheckbox();
//make_template_selection = Dialog.getCheckbox();
if(Dialog.getRadioButton()=="Create a new template in the image") make_template_selection = true;
else make_template_selection = false;
save_images = Dialog.getCheckbox();
analyze_all = Dialog.getCheckbox();
pause = Dialog.getCheckbox();

pixel_size=dpi/25.4;
min_area = pow(min_diameter/2,2)*PI;
max_area = pow(max_diameter/2,2)*PI;
current_image_nr = 0;

setBatchMode(true);


do {
	for(s=0;s<nr_series;s++) {	//loop over all series
		if (analyze_all==true) {
			run("Close All");
			file_name = image_list[current_image_nr];	//retrieve file name from image list
			run("Bio-Formats Importer", "open=["+dir+file_name+"] autoscale color_mode=Composite view=Hyperstack stack_order=XYCZT series_"+s+1);
			roiManager("Reset");
			roiManager("Show None");
			setBatchMode("show");
		}
		original = getTitle();
		getDimensions(width, height, channels, slices, frames);
		getPixelSize(unit, pixel_width, pixel_height);
		if(pixel_width != pixel_height) exit("Shame on you! Go and buy a camera with square pixels!");
//		pixel_size = pixel_width;	//For now: Use the hardcoded value defined above.
		run("Properties...", "unit=microns pixel_width="+pixel_size+" pixel_height="+pixel_size+" voxel_depth=0");

		if(bitDepth==24 || channels==3) {
			RGB==true;
			split_RGB_image(original);
		}
		else {
			run("8-bit");
			selectWindow(original);
			//rename(original+"_original");
			original = getTitle();
			run("Duplicate...", "title="+original+"_colonies");
			run("Duplicate...", "title="+original+"_wells");	//duplicate images - not very elegant, but easier to cope with both RGB images that can be split and grayscale images that cannot be split.
		}
		getDimensions(width, height, channels, slices, frames);

		well_diameter = well_diameter/pixel_size;

		width = getWidth;
		height = getHeight;
		//well_size_px = well_diameter/pixel_size;	//currently already in pixels

		if(automatic_detection == true) detect_wells(original+"_wells");
		else {
			//adapt from earlier versions
		}

		if(automatic_detection == true) { 
			run("Labels...", "color=white font=16 show use bold");
			for(roi_nr=0;roi_nr<roiManager("count");roi_nr++) {
				roiManager("Show All");
				selectWindow(original+"_colonies");
				roiManager("Select",roi_nr);

				count_colonies(original+"_colonies");
			}
		}
		else if(realtime_colonies == true) {
			showMessage("Left-click on the image to activate realtime quantification(TM). End with Space bar.");
			shift=1;
			ctrl=2; 
			rightButton=4;
			alt=8;
			leftButton=16;
			insideROI = 32; // requires 1.42i or later

			x2=-1; y2=-1; z2=-1; flags2=-1;
			setOption("DisablePopupMenu", true);
			while (!isKeyDown("space")) {
				selectWindow(original);
				getCursorLoc(x, y, z, flags);
				if (x!=x2 || y!=y2 || z!=z2 || flags!=flags2) {
					s = " ";
					if (flags&leftButton!=0) {
						run("Remove Overlay");
						roiManager("Show None");
						if(roiManager("count")>0) {
							roiManager("select",roiManager("count")-1);
							roiManager("delete");
						}
						makeOval(x-well_diameter/2, y-well_diameter/2, well_diameter, well_diameter);
						roiManager("add");
						roiManager("select",roiManager("count")-1);
						count_colonies();
						run("Add Selection...");
						roiManager("select",roiManager("count")-1);
						roiManager("rename",""+nResults);
						run("Labels...", "color=white font=16 show use bold");
						roiManager("Show All with labels");
						//roi_name = Roi.getName;
						//Roi.setName(roi_name+" - "+nResults);	//add the number of colonies to the roi name
						//Roi.setStrokeColor(roi_color);	//reset the active ROI color
					}
					if (flags&rightButton!=0) s = s + "<right>";
 					if (flags&shift!=0) s = s + "<shift>";
					if (flags&ctrl!=0) s = s + "<ctrl> ";
					if (flags&alt!=0) s = s + "<alt>";
					if (flags&insideROI!=0) s = s + "<inside>";
					// print(x+" "+y+" "+z+" "+flags + "" + s);
				}
				x2=x; y2=y; z2=z; flags2=flags;
				wait(10);
			}
			setOption("DisablePopupMenu", false);
		}

//		if(automatic_detection == true) {
			run("From ROI Manager");
			roiManager("Show All with labels");
			run("Labels...", "color=white font=16 show use bold");
		
			selectWindow("Summary");
			saveAs("Results", savedir+"\\"+substring(file_name,0,lastIndexOf(file_name,".")-1)+"_analyzed.xls");
			run("Close");
			if(save_images==true) {
				//run("Hide Overlay");
				//saveAs("png", savedir+"\\"+File.nameWithoutExtension);
				//run("Show Overlay");
				run("Labels...", "color=white font=64 show use bold");	//somehow the fonts are really small in the flattened image
				run("Flatten");
				saveAs("tif", savedir+"\\"+substring(file_name,0,lastIndexOf(file_name,".")-1)+"_analyzed");
				selectWindow(original);
				run("Labels...", "color=white font=16 show use bold");
			}
//		}
		if(analyze_all==true && pause==true  && current_image_nr<image_list.length) waitForUser("Click OK to continue");
		else print("Finished");
		if(current_image_nr<image_list.length && analyze_all==true) {
			if(s+1==nr_series) current_image_nr++;
		}
		if(analyze_all==true && current_image_nr<image_list.length) close();

	} //end of for loop over series
	
} while(analyze_all==true && current_image_nr<image_list.length)

restoreSettings();







function split_RGB_image(image) {
	selectWindow(image);
	rename("temp_image");
	run("RGB Color");
	rename(image);
	close("temp_image");
	run("Colour Deconvolution", "vectors=[H DAB] hide");

	imageCalculator("Average", image+"-(Colour_2)", image+"-(Colour_3)");
	rename(image+"_wells");
	run("8-bit");
	setBatchMode("show");
	imageCalculator("Difference", image+"-(Colour_1)", image+"_wells");
	rename(image+"_colonies");
	run("8-bit");
	run("Invert");
	setBatchMode("show");
	close(image+"-(Colour_3)");
	selectWindow(image);
//	rename(image+"_original");
	setBatchMode("show");
}





function detect_wells(well_image) {

	selectWindow(well_image);
	run("Set Measurements...", "area mean standard min limit redirect=None decimal=3");
	roiManager("Reset");

	original = getImageID();
	run("Select None");

//	setBatchMode(true);

	run("Scale...", "x="+scale+" y="+scale+" interpolation=Bilinear average create");
	image = getImageID();
	run("Remove Outliers...", "radius="+size_outliers+" threshold=0 which=Dark");

	if(background_correct==true) background_value = replace_background(image);
//waitForUser("background removed");
	roiManager("Reset");

	selectImage(image);
//run("Southwest");		Test: adding shadows...
	run("Find Edges");
	image_edges = getImageID();
	dim = pad_image_edges(image_edges, dim, "Top-Left");
	normalized_image = normalize_image(image_edges);

	//Template creation or loading (TO DO)
	if (make_template_selection == true) {
		selectImage(original);
		template = create_template();
		pad_image_edges(template, dim, "Center");	//send dimension retreived from the original image 
		template_title = normalize_image(template);
	}
	else {
		open(template_path);
		rename("template_normalized");
		template_title = getTitle();
	}


	//Cross-correlate the image with a template
	//TO DO: Or make binary first - check which one works better
	run("FD Math...", "image1=["+normalized_image+"] operation=Correlate image2=["+template_title+"] result=CC do");

	run("Variance...", "radius="+peak_radius);
	run("Gaussian Blur...", "sigma="+2*peak_radius);
	setAutoThreshold("Otsu dark");
	run("Find Maxima...", "noise=0 output=[Point Selection] above");
	getSelectionCoordinates(x_, y_);
	start = getTime();
	//print(x_.length+" wells detected");
	maxima_to_delete = get_faulty_maxima(x_,y_);
	x_ = Array_delete_elements(x_, maxima_to_delete);
	y_ = Array_delete_elements(y_, maxima_to_delete);
	//print(x_.length+" wells detected after correction for faulty maxima.");
	//Iteratively change the threshold until the correct number of maxima are found.
	while (x_.length > nr_wells && getTime-start<timeout) {
		getThreshold(th_min,th_max);
		setThreshold(th_min + 2000000*(abs(x_.length-nr_wells)),th_max);
		run("Find Maxima...", "noise=0 output=[Point Selection] above");
		getSelectionCoordinates(x_, y_);
		//print(x_.length+" wells detected; threshold = "+th_min);
		maxima_to_delete = get_faulty_maxima(x_,y_);
		x_ = Array_delete_elements(x_, maxima_to_delete);
		y_ = Array_delete_elements(y_, maxima_to_delete);
		//print(x_.length+" wells detected after correction for faulty maxima.");
	}
	while(x_.length < nr_wells && getTime-start<timeout) {
		showStatus("Optimizing threshold for well detection");
		getThreshold(th_min,th_max);
		setThreshold(th_min - 2000000*(abs(x_.length-nr_wells)),th_max);	//change threshold slower if the number of detected wells is close to the expected number
		run("Find Maxima...", "noise=0 output=[Point Selection] above");
		getSelectionCoordinates(x_, y_);
		//print(x_.length+" wells detected; threshold = "+th_min);
		maxima_to_delete = get_faulty_maxima(x_,y_);
		x_ = Array_delete_elements(x_, maxima_to_delete);
		y_ = Array_delete_elements(y_, maxima_to_delete);
		//print(x_.length+" wells detected after correction for faulty maxima.");
	}
	if(getTime-start>timeout) print("timeout (>"+timeout/1000+" s) when optimizing thresholds.");
	if(x_.length != nr_wells) {
		print("Warning: "+x_.length+" wells detected, but "+nr_wells+" expected.");
		error=true;
	}

	get_well_layout(x_,y_);
	
	//error handling
	if(error == true) {
		print("Warning: Unable to faithfully detect the well layout!");
		print("Well layout found (x-y): "+nr_wells_x+" by "+nr_wells_y);
	}
	if(error == false && nr_wells_y * nr_wells_x == nr_wells) {	// no error case
		print("Well layout found (x-y): "+nr_wells_x+" by "+nr_wells_y);
	}
	else if (x_.length == nr_wells && nr_wells_y * nr_wells_x == nr_wells && error == true) {
		setTool("oval");
		selectImage(original);
		roiManager("Show All");
		setBatchMode("show");
		waitForUser("Automatic well detection failed. Please move wells manually ok to continue.");
		//TO DO: Sorting doesn't work
	}
	else if (x_.length != nr_wells && nr_wells_y * nr_wells_x == nr_wells) {
		setTool("oval");
		selectImage(original);
		roiManager("Show All");
		setBatchMode("show");
		waitForUser("Automatic well detection failed. Please add or remove wells manually to or from the ROI manager and click OK.");
		//TO DO: Sorting doesn't work
	}
	else {
		setTool("oval");
		roiManager("Show All");
		selectImage(original);
		showMessage("Automatic well detection REALLY failed. Please manually add well selections by clicking in the middle of each well.");
		manual_well_selection();
		getSelectionCoordinates(x_, y_);
		x_ = scale_array(x_, scale);	//scaling necessary because the manual selection is done on the original image
		y_ = scale_array(y_, scale);
		get_well_layout(x_,y_);
	}
	//TO DO: If well detection fails, but the layout is ok, search for maxima in the CC window at the expected empty location ± well_diameter/2
	//Also, prevent 2 maxima closer together than well_diameter/2
	sort_ROIs(x_,y_);

	selectImage(original);
	roiManager("Show All");
}


function get_faulty_maxima (x,y) {
	faulty_maxima = newArray(9999);	//allow a maximum of 9999 detections
	m=0;
	for(i=0;i<x.length;i++) {
		distances = newArray(x.length); 
		for(j=i;j<x.length;j++) {
			distances[j]=sqrt( ((x[j]-x[i])*(x[j]-x[i])) + ((y[j]-y[i])*(y[j]-y[i])) );
		}
		//Array.print(distances);
		for(j=i;j<x.length;j++) {	//start at i, because the first maxima are the highest, and we want to throw out the lower maximum that is close
			if (distances[j]!=0 && distances[j]<200) {
				//print(i+": BINGO! Maximum "+j+" should be kicked out! (distance: "+distances[j]+")");
				faulty_maxima[m] = j;
				m++;
			}
		}
		//Do not allow wells touching the edges
		if(exclude_edge_wells == true) {
			if(x[i]<well_diameter/2*scale || x[i]>((width-well_diameter/2)*scale) || y[i]<well_diameter/2*scale || y[i]>((height-well_diameter/2)*scale)) {
				//print("Edge maximum found at "+x[i]+","+y[i]);
				faulty_maxima[m] = i;
				m++;
				//setPixel(x[i], y[i], 65536);
			}
		}
	}
	faulty_maxima = Array.trim(faulty_maxima,m);
	faulty_maxima = Array.sort(faulty_maxima);

	return faulty_maxima;
}

function Array_delete_elements(array, to_delete) {
	n=0;
	for(i=0;i<to_delete.length;i++) {
		if(i==0) {
			//Array.print(array);
			Array.rotate(array,-to_delete[i]+n);
			//Array.print(array);
			array = Array.slice(array,1,array.length);
			//Array.print(array);
			Array.rotate(array,to_delete[i]-n);
			//Array.print(array);
			n++;
		}
		else if(to_delete[i]!=to_delete[i-1]) {	//Check for double occurences
			//Array.print(array);
			Array.rotate(array,-to_delete[i]+n);
			//Array.print(array);
			array = Array.slice(array,1,array.length);
			//Array.print(array);
			Array.rotate(array,to_delete[i]-n);
			//Array.print(array);
			n++;
		}
	}
	return array;
}



function count_colonies(colonies_image) {
	selectWindow(colonies_image);
	Roi.setStrokeColor("green");	//color the active well green
	well_name = Roi.getName;
	run("Duplicate...", "title="+well_name);
//showimage("before subtraction");
	run("Subtract Background...", "rolling="+rolling_ball_radius/pixel_size+" light sliding");
//showimage("after subtraction");
	if(smooth_factor>0) run("Median...", "radius="+smooth_factor);
	run("Auto Local Threshold", "method=Mean radius="+radius/pixel_size+" parameter_1="+parameter1+" parameter_2=0");
	if(watershed==true) run("Watershed");

	run("Restore Selection");
//showimage("before clear");
	run("Clear Outside");
//showimage("before analyze particles");
	if(exclude==true) run("Analyze Particles...", "size="+min_area+"-"+max_area+" circularity="+min_circ+"-1.00 display clear summarize exclude show=[Bare Outlines]");
	else run("Analyze Particles...", "size="+min_area+"-"+max_area+" circularity="+min_circ+"-1.00 display clear summarize show=[Bare Outlines]");
	rename("outlines_"+well_name);
//showimage("outlined colonies");
	run("Invert");
	run("Cyan");

//	selectWindow("Results");
	Area = newArray(nResults);
	for(i=0;i<nResults;i++) {
			Area[i] = getResult("Area",i);
	}
	Array.getStatistics(Area, min, max, mean, stdDev);
	print("\\Update:\n");	//remove line printed by auto local threshold routine
	print("\\Update:"+original+":"  +nResults+" colonies with mean diameter "+d2s(2*sqrt(mean/PI),0)+" "+unit);

	selectWindow(colonies_image);
	getSelectionBounds(x,y,width,height);
	roi_name = Roi.getName;
	Roi.setName(roi_name+" - "+nResults);	//add the number of colonies to the roi name
	Roi.setStrokeColor(roi_color);	//reset the active ROI color
	selectWindow(original);
	run("Add Image...", "image=outlines_"+well_name+" x="+x+" y="+y+" opacity=100 zero");
	close(well_name);
	close("outlines_"+well_name);
}


//TO DO: Make this into a macro with a dedicated button in the ImageJ window
function create_template() {
	setTool("Oval");
	setKeyDown("shift");
	waitForUser("Create a selection around a typical well (hold Shift while dragging) and press OK.");
	setKeyDown("none");
	run("Duplicate...", "Title=template_unscaled");
	run("Remove Outliers...", "radius=30 threshold=50 which=Dark");
	run("Find Edges");
	run("Clear Outside");
	run("Scale...", "x="+scale+" y="+scale+" interpolation=Bilinear average create title=template");
	template = getImageID();
	close("template_unscaled");
	return template;
}


function manual_well_selection() {
	roiManager("Reset");
	setTool("multipoint");
	waitForUser("Please mark the center of each well, "+nr_wells+" in total, and click OK.");
	getSelectionCoordinates(x_manual, y_manual);
	while(x_manual.length != nr_wells) {
		waitForUser("Select exactly "+nr_wells+" wells. Use Alt-leftclick to remove a point.");
	}
}



function get_well_layout(x_, y_) {
	nr_wells_y = get_steps(y_);
	nr_wells_x = get_steps(x_);

	selectImage(image);
	for(i=0;i<x_.length;i++) {
		makeOval((x_[i]/scale)-well_diameter/2, (y_[i]/scale)-well_diameter/2, well_diameter, well_diameter);
		Roi.setStrokeColor(roi_color);
		roiManager("add");
		run("Select None");
	}
}


function replace_background(image) {	//replace bright background value by the mean in order to minimize edges
	selectImage(image);
	roiManager("Reset");
	setAutoThreshold("Minimum");	//Was Default, but too lenient in some cases.
	List.setMeasurements();
	mean = List.getValue("Mean");
	setAutoThreshold("Minimum dark");
	run("Analyze Particles...", "size=50000000-Infinity show=Nothing include add");
	if(roiManager("count")>1) {
		roiManager("Select All");
		roiManager("Combine");
		roiManager("Reset");
		roiManager("add");
		roiManager("Select",0);
		roiManager("rename", "background");
		changeValues(0,65535,mean);
	}
	run("Select None");
	resetThreshold();
	roiManager("Show None");
	return mean;
}


function normalize_image(image) {	//subtract mean and divide by stddev to normalize the cross correlation
	selectImage(image);
	run("Duplicate...", "title="+image+"_normalized");
	normalized_image = getTitle();
	run("32-bit");
	List.setMeasurements();
	mean = List.getValue("Mean");
	stddev = List.getValue("StdDev");
	run("Subtract...", "value="+mean);
	run("Divide...", "value="+stddev);
	resetMinAndMax();
	run("Enhance Contrast", "saturated=0.35");
	return normalized_image;
}


function pad_image_edges(image, size, location) {
	selectImage(image);
	width = getWidth();
	height = getHeight();
	w=1;
	h=1;
	while(width>pow(2,w) || size>pow(2,w)) w++;
	while(height>pow(2,h) || size>pow(2,w)) h++;
	dim = pow(2,maxOf(w,h));
	//print("dimension of FFT: "+dim);
	run("Canvas Size...", "width="+dim+" height="+dim+" position="+location+" zero");
	return dim;
}


function get_steps(array) {
	array_sort = Array.copy(array);
	array_sort = Array.sort(array_sort);
	array_diff = newArray(array.length-1);
	for(i=0;i<array_diff.length;i++) {
		array_diff[i] = array_sort[i+1] - array_sort[i];
	}
	Array.getStatistics(array_diff, min, max, mean, stdDev);
	array_maxima = Array.findMaxima(array_diff, stdDev, 1);	//check for steps larger than stdDev
	Array.sort(array_maxima);

	//check the location of the steps
	nr_occurences = array_maxima.length+1;
	//print(array.length);
	//print(array_maxima.length);
	//print(nr_occurences);
	step = array.length/nr_occurences;
	j=0;
	for(i=step-1 ; i<array.length-1 ; i=i+step) {
		//print(i);
		//print(array_maxima[j]);
		if(i != array_maxima[j]) error=true;
		j++;
	}
	return nr_occurences;
}


function sort_ROIs(x_ , y_) {
	rank_x_ = Array.rankPositions(x_);	//First x, then y.
	//Array.print(y_);
	//Array.print(x_);
	//Array.print(rank_x_);
	x_sort = Array.copy(x_);
	x_sort = Array.sort(x_sort);
	y_sort = Array.copy(y_);
	y_sort = Array.sort(y_sort);

	y_coords_ = newArray(nr_wells_y);
	for(i=0;i<x_.length;i++) {
		roiManager("select",rank_x_[i]);
		roiManager("rename",fromCharCode(65+floor(i/nr_wells_y)));	//add letters to the well name
	}

	for(i=0;i<x_.length;i++) {
		k = i%nr_wells_y;
		y_coords_[k] = y_[rank_x_[i]];	//fill the temp array y_coords_ with the correct values for each column
	
		if(k==nr_wells_y-1) {
			//Array.print(y_coords_);
			rank_y_coords_ = Array.rankPositions(y_coords_);
			//Array.print(rank_y_coords_);
			for(j=0;j<y_coords_.length;j++) {
				roiManager("select",rank_x_[i+j-nr_wells_y+1]);
				//print(rank_x_[i+j-nr_wells_y+1]+" -> "+rank_y_coords_[j]);
				roiManager("rename",call("ij.plugin.frame.RoiManager.getName", rank_x_[i+j-nr_wells_y+1])+IJ.pad(rank_y_coords_[j]+1, 2));	//numbers
			}
		}
	}
	roiManager("sort");
}


function scale_array(array, scale) {
	for(i=0;i<array.length;i++) {
		array[i] = array[i] * scale;
	}
	return array;
}


//Returns the number of times the value occurs within the array
function occurencesInArray(array, value) {
	count=0;
	for (a=0; a<lengthOf(array); a++) {
		if (array[a]==value) {
			count++;
		}
	}
	return count;
}

function showimage(string) {
	setBatchMode("show");
	waitForUser(string);
}