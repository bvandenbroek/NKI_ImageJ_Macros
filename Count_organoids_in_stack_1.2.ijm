

var nr_series;
var file_name;
var format;

var manual_threshold = true;
var median_radius = 2;
var kernel = 11;
var autothreshold = false;
choice_projection_array = newArray("Stack focuser plugin", "Minimum Intensity projection");
var choice_projection = "Stack focuser plugin";	//default
var remove_combine_ROIs = true;
var process_all_files = false;
var pause = false;
var draw = false;
var scale = 0.5;
var background_radius = 50	//radius for background subtraction (in pixels)
var min_size = 50;
var edge = 80;				//Determines the search circle/ellipse size from the edge
unit="um";

saveSettings();

print("\\Clear");
run("Clear Results");

if(nImages>0) run("Close All");
path = File.openDialog("Select any file in the folder to be processed");

setBatchMode(true);

run("Bio-Formats Macro Extensions");
//run("Bio-Formats Importer", "open=["+path+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_1");
Ext.setId(path);
Ext.getFormat(path, format);
print(format);

//getDimensions(width, height, channels, slices, frames);

file_name = File.getName(path);
dir = File.getParent(path)+"\\";
savedir = dir+"\\results\\";
if(!File.exists(savedir)) File.makeDirectory(savedir);
extension = substring(file_name, lastIndexOf(file_name,".")+1);

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

//Dialog
Dialog.create("Background subtraction: options");
	Dialog.addSlider("Stack focuser kernel size (pixels)", 3, 40, kernel);
	Dialog.addNumber("Median radius for filtering organoids",median_radius,1,3,"pixels");
	Dialog.addChoice("Projection method", choice_projection_array, choice_projection);
	Dialog.addNumber("Minimum organoid diameter",min_size,0,5,unit);
	Dialog.addCheckbox("Automatic thresholding?", autothreshold);
	Dialog.addCheckbox("Manually remove/combine found organoids?", remove_combine_ROIs);
	Dialog.addCheckbox("Process all "+image_list.length+" ."+extension+" ("+format+") files.?", process_all_files);
	Dialog.addCheckbox("Pause after each image? ", pause);
	Dialog.addCheckbox("Draw the selected region in the image? ", draw);
Dialog.show;
kernel=Dialog.getNumber();
median_radius = Dialog.getNumber();
choice_projection = Dialog.getChoice();
min_size = Dialog.getNumber();
autothreshold=Dialog.getCheckbox();
remove_combine_ROIs=Dialog.getCheckbox();
process_all_files=Dialog.getCheckbox();
pause=Dialog.getCheckbox();
draw=Dialog.getCheckbox();


print("\\Clear");
print("Directory contains "+file_list.length+" files, of which "+image_list.length+" ."+extension+" ("+format+") files.");
run("Set Measurements...", "area mean integrated median stack limit redirect=None decimal=3");
current_file_nr=0;
do {
	if(process_all_files==true) {
		run("Close All");
		file_name = image_list[current_file_nr];		//retrieve file name from image list
	}
	else file_name = File.getName(path);
	Ext.setId(dir+file_name);
	Ext.getSeriesCount(nr_series);

	//run("Clear Results");
	roiManager("Reset");
	
	for(i=0;i<nr_series;i++) {
	print("Processing file "+current_file_nr+1+"/"+image_list.length+", series "+i+1+"/"+nr_series+": "+dir+file_name+"...");
	run("Bio-Formats Importer", "open=["+dir+file_name+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+i+1);
	
	name = getTitle();
	name = replace(name,"\\/","-");	//replace slashes by dashes in the name
	name = substring(name, 0,lengthOf(name)-lengthOf(extension)-1);	//remove extension
	rename(name);

	run("Scale...", "x="+scale+" y="+scale+" z=1.0 interpolation=Bilinear average process create");
	close(name);
	
	frame_interval = Stack.getFrameInterval();
	getDimensions(width, height, channels, slices, frames);
	
	getVoxelSize(width, height, depth, unit);
	if(choice_projection == "Stack focuser plugin") run("Stack Focuser ", "enter="+kernel);
	else run("Z Project...", "projection=[Min Intensity]");
	setVoxelSize(width, height, depth, unit);

	focussed = getTitle();
	run("Duplicate...", "title=filtered duplicate");
	run("Remove Outliers...", "radius=6 threshold=200 which=Dark");	//remove crap
	run("Subtract Background...", "rolling="+background_radius+" light sliding");
	run("Median...", "radius="+median_radius);

	setBatchMode("show");

	getDimensions(width, height, channels, slices, frames);
	makeOval(edge,edge,width-(2*edge),height-(2*edge));

	run("Enhance Contrast", "saturated=0.35");
	run("8-bit");
//	run("Auto Local Threshold", "method=Mean radius=50 parameter_1=25 parameter_2=0");
	run("Threshold...");
	setAutoThreshold("Li");
	if(autothreshold==false) {
		setBatchMode("show");
		waitForUser("Change threshold if necessary");
	}
	wait(25);

	run("Analyze Particles...", "size="+scale*scale*PI*min_size*min_size/4+"-Infinity exclude include add");

	if(remove_combine_ROIs==true) combine_ROIs(focussed);
	
	selectWindow(focussed);
	if(roiManager("count")>0) {
		for(j=0;j<roiManager("count");j++) {
			roiManager("Select",j);
			run("Measure");
			setResult("File name", nResults-1, file_name);
			updateResults();
			roiManager("Set Line Width", 2);
			//	roiManager("Set Color", "White");
			if(draw==true) run("Draw", "slice");
		}
	}
	else {	//if no ROIs are found
		run("Measure");
		setResult("Area", current_file_nr, 0);
		setResult("Mean", current_file_nr, 0);
		setResult("X", current_file_nr, 0);
		setResult("Y", current_file_nr, 0);
		setResult("Major", current_file_nr, 0);
		setResult("Minor", current_file_nr, 0);
		setResult("Angle", current_file_nr, 0);
		setResult("Circ.", current_file_nr, 0);
		setResult("IntDen", current_file_nr, 0);
		setResult("RawIntDen", current_file_nr, 0);
		setResult("AR", current_file_nr, 0);
		setResult("Round", current_file_nr, 0);
		setResult("Solidity", current_file_nr, 0);
		setResult("File name", current_file_nr, file_name);
	}
	saveAs("Tiff",savedir+name+"_focused");
	
	if(pause==true) {
		setBatchMode("show");
		if(roiManager("count")>0) waitForUser("Click OK to continue.");
		else waitForUser("No organoid found! Click OK to continue.");
	}		

		saveAs("results", savedir+name+".xls");
		close();
	}
	current_file_nr++;
} while (process_all_files==true && current_file_nr<image_list.length);


selectWindow("Results");
updateResults();
if(process_all_files==true) saveAs("results",savedir+"Results_"+image_list.length+"_files.txt");

if(process_all_files==true) {
	run("Close All");
	print("Finished analyzing "+image_list.length+" images.");
}

//setBatchMode(false);

restoreSettings();










function combine_ROIs(image1) {

shift=1;
ctrl=2; 
rightButton=4;
alt=8;
leftButton=16;
insideROI = 32;

flags=-1;
//x2=-1; y2=-1; z2=-1; flags2=-1;

selectWindow(image1);
roiManager("Show All without labels");
setOption("DisablePopupMenu", true);
setBatchMode(true);
resetMinAndMax();
setBatchMode("show");
showMessage("Combine ROIs of split nuclei. Select ROIs with shift-left mouse button and right-click to merge them.\nLeft clicking while holding CTRL deletes a ROI.\nPress CTRL-shift-leftclick bar when finished editing. This information will be printed to the log window.");
print("Combine ROIs of split nuclei. Select ROIs with shift-left mouse button and right-click to merge them.\nLeft clicking while holding CTRL deletes a ROI.\nPress CTRL-shift-leftclick bar when finished editing.");

for(i=0;i<roiManager("Count");i++) {	//unselect all ROIs
	roiManager("Select",i);
	roiManager("Rename", "ROI "+i);
	Roi.setProperty("selected",false);
	Roi.setStrokeColor("cyan");
}

//while(!isKeyDown("space")) {
run("Remove Overlay");
setTool("freehand");
roiManager("Show All Without Labels");
setOption("DisablePopupMenu", true);
while(flags!=19) {						//exit loop by pressing ctrl-shift-leftclick
	getCursorLoc(x, y, z, flags);
//	if (x!=x2 || y!=y2 || z!=z2 || flags!=flags2) {

//Check if a freehand selection is present. Do something with it.
//	if(selectionType>0 && selectionType<=3) {
//		print("selection found");
//	}

	if(flags==17 || flags==18)	{
		for(i=0;i<roiManager("Count");i++) {
			roiManager("Select",i);
			if(Roi.contains(x, y)==true) {
			selected = Roi.getProperty("selected");
				//click to select a single ROI
				if(flags==17 && selected==false) {		//select ROI
					//print("selecting ROI "+i);
					Roi.setStrokeColor("red");
					Roi.setProperty("selected",true);
				}
				else if(flags==17 && selected==true) {	//deselect ROI
					//print("deselecting ROI "+i);
					Roi.setStrokeColor("cyan");
					Roi.setProperty("selected",false);
				}
				else if(flags==18) {	//delete ROI
					roiManager("Delete");
					for(j=0;j<roiManager("Count");j++) {	//deselect all ROIs and rename
						roiManager("Select",j);
						roiManager("Rename", "ROI "+j);
					}
				}
			}
		}
	roiManager("Deselect");
	run("Select None");
	updateDisplay();
	}
	if(flags==4) {
		selected_ROI_array = newArray(roiManager("Count"));	//create array with indices of selected ROIs
		j=0;
		for(i=0;i<roiManager("Count");i++) {
			roiManager("select",i);
			selected = Roi.getProperty("selected");
			if(selected==true) {
				selected_ROI_array[j] = i;
				j++;
				//print(j);
			}
		}
		//check if more than one ROI is selected. If yes, combine the selected ROIs and update the list
		selected_ROI_array = Array.trim(selected_ROI_array,j);
		//print(selected_ROI_array.length + " ROIs selected");
		if(selected_ROI_array.length > 1) {
			//print("combining ROIs:");
			Array.print(selected_ROI_array);
			roiManager("Select",selected_ROI_array);
			roiManager("Combine");
			roiManager("Update");
//			for(i=1;i<selected_ROI_array.length;i++) {	
			to_delete_array = Array.copy(selected_ROI_array);														//selecting and deleting redundant ROIs
			to_delete_array = Array.slice(selected_ROI_array,1,selected_ROI_array.length);	//create array without the first element
				roiManager("Deselect");
				//print("deleting ROIs:");
				Array.print(to_delete_array);
				roiManager("select", to_delete_array);
				roiManager("Delete");
			roiManager("Select",selected_ROI_array[0]);
			//print("repairing ROI "+selected_ROI_array[0]);
			run("Enlarge...", "enlarge=1 pixel");			//remove wall between ROIs by enlarging and shrinking with 1 pixel
			run("Enlarge...", "enlarge=-1 pixel");
			roiManager("Update");
			
			setKeyDown("none");
			
			for(i=0;i<roiManager("Count");i++) {	//deselect all ROIs and rename
			roiManager("Select",i);
			roiManager("Rename", "ROI "+i);
			Roi.setProperty("selected",false);
			Roi.setStrokeColor("cyan");
			}
		}
	roiManager("Deselect");
	run("Select None");
	updateDisplay();	//doesn't work...?
	}
//	x2=x; y2=y; z2=z; flags2=flags;

	wait(50);
}	//end of while loop

}
