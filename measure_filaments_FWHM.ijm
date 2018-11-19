// Macro to measure the average thickness of filaments in an image.
// Bram van den Broek, The Netherlands Cancer Institute, 2013-2018
// b.vd.broek@nki.nl
//

#@ Integer (label = "line width", style = "spinner", value=50, min=1, max=1000) line_width
#@ Integer (label = "fitted segments length", style = "spinner", value=100, min=10, max=10000) block_size
#@ Float (label = "Minimum R-squared value to accept fit", style = "spinner", value=0.90, min=0, max=1) min_R_sq
#@ Float (label = "Maximum filament FWHM to accept fit (um)", style = "spinner", value=0.2, min=0, max=1000) max_width
#@ Boolean (label = "Show montage of straigtened filaments") show_montage

var filament_channel = 1		//In case of mutichannel image
var interpolation_interval = 10;//interval for interpolating the lines - necessary for getting the profile.
var zoom = 400;					//Final zoom of segment images
var overlay_color = "#7fff0000"	//alpha,red,green,blue

//Arrays to hold the data
var filament_array = newArray(999);
var width_array = newArray(999);
var R_squared_array = newArray(999);
var amplitude_array = newArray(999);

var segment_nr=0;
var nr_fits=0;
var first_fit=true;
var continue_drawing = true;
var lines=1;
var fit_stack;

var verbose=false;

saveSettings();

setOption("Bicubic", true);
roiManager("reset");
print("\\Clear");
run("Remove Overlay");
setTool("polyline");
run("Overlay Options...", "stroke=red width=1 fill=#88880000 set");

run("Close All");
path=File.openDialog("Select file...");
open(path);

//////////////////
//CONVERTING TO 8 BIT - ONLY FOR RGB IMAGES - REMOVE FOR RAW DATA
//if(bitDepth!=8) run("8-bit");
//////////////////


image = getTitle();
getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pw, ph);
dir = getDirectory("image");
file_name_without_extension = File.nameWithoutExtension;
if (File.exists(dir+file_name_without_extension+"_ROIs.zip")) {		//load ROIs of previously lines in this file
	roiManager("Open", dir+file_name_without_extension+"_ROIs.zip");
	lines=roiManager("count")+1;	//+1 because this number is the next line
	load_ROIs = getBoolean(""+lines-1+" lines have already been drawn for this image.\nDo you want to keep them?");
	if (load_ROIs==false) {
		roiManager("reset");
		lines=1;
	}
}

run("Line Width...", "line="+line_width);

if(verbose==false) setBatchMode(true);

run("Duplicate...", "duplicate channels="+filament_channel);
rename("f*i*l*a*m*e*n*t*s_"+image);
filaments = getTitle();
selectWindow(image);

//Draw lines to be analyzed
while(continue_drawing==true) {
	waitForUser("Draw line "+lines+" and press OK");
	if(selectionType()!=-1) {
		setSelectionName("line_"+lines);
		run("Interpolate", "interval="+interpolation_interval);		//smoothing the lines
		roiManager("Add");
		roiManager("Show All with labels");
		continue_drawing = getBoolean("continue to draw lines?");
		lines++;
	}
	else {
		if (roiManager("count")>0) continue_drawing = !getBoolean("Proceed calculation with "+roiManager("count")+" lines?");
		else exit("No lines drawn. Exiting macro.");
	}
}
getMinAndMax(min_level,max_level);
getLut(reds, greens, blues);		//retreive LUT

lines=roiManager("count");
x_points = create_ramp_array(line_width-1,pw);	//create x-data for the fit

newImage("segment_stack", "32-bit black", line_width, block_size, 1);
run("Properties...", "unit=nm pixel_width="+pw+" pixel_height="+ph);

for(i=0;i<lines;i++) {
	selectWindow(filaments);
	roiManager("select",i);
	run("Straighten...");
	rename("filament_"+i);
	run("Rotate 90 Degrees Right");

	fit_and_align_segments(i);
	

}
selectWindow("segment_stack");
run("Select None");
setSlice(1);
setBatchMode("show");
setLut(reds, greens, blues);
setMinAndMax(min_level,max_level);
run("Set... ", "zoom="+zoom);

//handle and print statistics
filament_array = Array.trim(filament_array, segment_nr);
width_array = Array.trim(width_array, segment_nr);
amplitude_array = Array.trim(amplitude_array, segment_nr);
R_squared_array = Array.trim(R_squared_array, segment_nr);
Array.getStatistics(width_array, min, max, mean, stdDev);
print("\n"+segment_nr+" of "+nr_fits+" segments were accepted, with R^2 > "+min_R_sq+" and FWHM < "+max_width);
print("Average width: "+mean+" +- "+stdDev+" (stddev)");
selectWindow(image);
run("Select None");

//Get the average profile and fit it
selectWindow("segment_stack");
setBatchMode("show");
run("Z Project...", "projection=[Sum Slices]");
rename("summed_segment");
setBatchMode("show");
run("Select None");
run("Set... ", "zoom="+zoom);
setLut(reds, greens, blues);
//setMinAndMax(min_level,max_level);
run("Select All");
y_points_filament = getProfile();
x_points = create_ramp_array(line_width,pw);
Fit.doFit("Gaussian", x_points, y_points_filament);
Fit.plot;
rename("Fitted_average_segment");
//rename("Fit of the averaged filament");
FWHM_avg = Fit.p(3)*2.35;			//FWHM =~ 2.35 * stddev
print("\nFWHM of the averaged filament = "+FWHM_avg+" "+unit+", R-squared="+Fit.rSquared);
setBatchMode("show");
selectWindow("averaged_segment");
run("Select None");
run("Out [-]");	//Otherwise not displayed correctly. Reason unknown.
run("In [+]");

//Print to result window
for(i=0;i<width_array.length;i++) {
	setResult("filament_nr",i,filament_array[i]);
	setResult("width",i,width_array[i]);
	setResult("amplitude",i,amplitude_array[i]);
	setResult("R^2",i,R_squared_array[i]);
}


if (show_montage==true) {
	//make montage of filaments
	run("Images to Stack", "method=[Copy (top-left)] name=all_filaments title=filament use");
	if(lines>1) run("Make Montage...", "columns="+lines+" rows=1 scale=1 increment=1 border=1 font=12");
	rename("montage_filaments");
	close("all_filaments");

	rename("montage");
	setBatchMode("show");
//	run("Rotate 90 Degrees Left");
}
setBatchMode("show");


//save results and ROIs to files
if (roiManager("count")>0) roiManager("Save", dir+file_name_without_extension+"_ROIs.zip");
saveAs("Results", dir+file_name_without_extension+"_results.txt");

setTool("hand");



restoreSettings();

////////////////////////////////////////////


function create_ramp_array(length,step) {
	array = newArray(length);
	for(i=0;i<length;i++) {
		array[i] = i*step;
	}
	return array;
}


function fit_and_align_segments(line_nr) {
	selectWindow("filament_"+line_nr);
	height=getHeight();
	y=0;
	do {
		selectWindow("filament_"+line_nr);
		makeRectangle(0, y, line_width-1, block_size);
		y_points = getProfile();
		initialGuesses = newArray(0,1,line_width/2,10*pw);
		Fit.doFit("Gaussian", x_points, y_points);
		Fit.plot;
		if (verbose==false) setBatchMode(true);
		if (first_fit==true) {
			fit_stack = getImageID;		//Get the ID of the plot window
			rename("Stack of fitted segments");
			setBatchMode("show");
			first_fit=false;
		}
		else {
			run("Copy");
			close();
			selectImage(fit_stack);
			run("Add Slice");
			run("Paste");
			setBatchMode("show");
		}
		nr_fits++;
		FWHM = Fit.p(3)*2.35;			//FWHM =~ 2.35 * stddev
		print("FWHM="+FWHM+" "+unit+", R-squared="+Fit.rSquared);
		if(Fit.rSquared>min_R_sq && FWHM<=max_width) {
			//Add fit results to data arrays
			filament_array[segment_nr] = line_nr+1;
			width_array[segment_nr] = FWHM;
			amplitude_array[segment_nr] = Fit.p(1);
			R_squared_array[segment_nr] = Fit.rSquared;
			segment_nr+=1;				//Only increment for accepted segments
			selectWindow("filament_"+line_nr);
			run("Copy");				//Copy-paste the approved line segment into the segment_stack
			selectWindow("segment_stack");
			if(segment_nr>1) {
				run("Add Slice");
				run("Paste");
			}
			else run("Paste");
			translation = (line_width-2)/2 - Fit.p(2)/pw;
			run("Translate...", "x="+translation+" y=0 interpolation=Bicubic slice");

		}
		else {
			selectWindow("filament_"+line_nr);
			   run("Add Selection...", "fill="+overlay_color);

			selectImage(fit_stack);
			getDimensions(w,h,ch,slices,frames);
			if(slices==1) {
				run("Add Slice");			//add slice and remove later to prevent the whole stack becoming red
				Stack.setSlice(1);
				run("Select All");
				   run("Add Selection...", "fill="+overlay_color);	//add red overlay to this segment in the fit window
				Stack.setSlice(2);
				run("Delete Slice");
			}
			else {
				run("Select All");
				   run("Add Selection...", "fill="+overlay_color);	//add red overlay to this segment in the fit window
			}

		}
		y+=block_size;
	} while (y<height);
	selectWindow("filament_"+line_nr);
	run("Select None");
	setLut(reds, greens, blues);
	setMinAndMax(min_level,max_level);
	run("Flatten");							//Flatten red overlays to an RGB image
	close("filament_"+line_nr);
	rename("filament_"+line_nr);
}

