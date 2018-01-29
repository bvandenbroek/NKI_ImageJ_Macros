roiManager("Reset");
run("Set Measurements...", "area mean standard min median limit redirect=None decimal=3");

setBatchMode(true);

original = getTitle();
getDimensions(width, height, channels, slices, frames);

backgrounds = newArray(channels);
true_values = newArray(channels);
false_values = newArray(channels);

print("------------------------");
for(i=1;i<=channels;i++) {
	selectWindow(original);
	Stack.setChannel(i);
	//run("k-means Clustering ...", "number_of_clusters=3 cluster_center_tolerance=0.0001000 enable_randomization_seed randomization_seed=48");
	
	//get background value (outside chromosomes) on the original image
	setAutoThreshold("Mean");
	List.setMeasurements("limit");
	backgrounds[i-1] = List.getValue("Median");
print("bg: "+backgrounds[i-1]);
	resetThreshold();
	run("Select None");
	
	//get median true value (positive chromosomes)
	run("Duplicate...", "duplicate title=ch"+i+" channels="+i);
	run("Median...", "radius=2");

	setAutoThreshold("Mean dark");
	getThreshold(th_leakthrough,th_max);
	run("Create Selection");
	setAutoThreshold("Otsu dark");
	getThreshold(th_positive,th_max);
		run("Create Selection");
		roiManager("Add");
		selectWindow(original);
		roiManager("Select", roiManager("count")-1);
		roiManager("Rename", "ch"+i+"_positive");
		roiManager("Update");
		List.setMeasurements("limit");
		true_values[i-1] = List.getValue("Median");
		run("Select None");
print("true: "+true_values[i-1]);

	//get median false value (leakthrough & cross-excitation in the other chromosomes)
	selectWindow("ch"+i);
	setThreshold(th_leakthrough,th_positive);
		run("Create Selection");
		roiManager("Add");
		selectWindow(original);
		roiManager("Select", roiManager("count")-1);
		roiManager("Rename", "ch"+i+"_negative");
		roiManager("Update");
		List.setMeasurements("limit");
		false_values[i-1] = List.getValue("Median");
		run("Select None");
print("false: "+false_values[i-1]);
	
	resetThreshold();
	close("ch"+i);
}
run("Select None");

setBatchMode(false);

//Normalize the image
//	selectWindow(original);
	waitForUser("Select window to normalize");
	run("Duplicate...", "duplicate title=["+original+"_normalized]");
	run("32-bit");
	for(i=1;i<=channels;i++) {
		Stack.setChannel(i);
		run("Macro...", "code=v=((v-("+false_values[i-1]+backgrounds[i-1]+"))/("+true_values[i-1]-false_values[i-1]-backgrounds[i-1]+")) slice");
	}
	setBatchMode("show");


/////// Steps to incorporate (earlier)
/*
 * Make mask of DAPI (now sum) (use a bit of median filtering) 
 * Multiply original with mask
 * Divide original by DAPI (now sum): result is a normalized image
 *  
 * 
 */

image = getTitle;
getDimensions(width, height, channels, slices, frames);
Inf = 1/0;
string = "";
for(i=1;i<=channels;i++) {
	selectWindow(image);
	Stack.setChannel(i);
	setThreshold(0.3,Inf);
	run("Create Mask");
	rename("color_"+i);
	string = string + "image"+i+"=color_"+i+" ";
}
run("Concatenate...", "  title=masks "+string);
run("Divide...", "value=255.000 stack");
setMinAndMax(0, 2);