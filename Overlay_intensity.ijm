//Overlay intensity image with ratio/FLIM image - the fast way
//Images need to be opened before

smooth=false;
radius=2;
//setBatchMode(true);

waitForUser("Select Ratio/FLIM image and adjust LUT and Min/Max");
ratio = getTitle();
if(smooth==true) {
	run("Mean 3D...", "x=1 y=1 z="+radius); //smooth in time
}
//setMinAndMax(2000, 3000);
//run("Physics Black");
run("Duplicate...", "title=splits duplicate");
run("RGB Color");
run("Split Channels");

waitForUser("Select intensity image or stack and adjust Min/Max");
//	run("Divide...", "value=64 stack");
//run("8-bit");
//setMinAndMax(0,50);	//brightness of intensity channel
intensity = getTitle();

run("Apply LUT", "stack");

imageCalculator("Multiply 32-bit stack", "splits (red)", intensity);
rename("Red");
imageCalculator("Multiply 32-bit stack", "splits (green)", intensity);
rename("Green");
imageCalculator("Multiply 32-bit stack", "splits (blue)", intensity);
rename("Blue");
run("Merge Channels...", "c1=Red c2=Green c3=Blue");
rename(ratio);
run("Set... ", "zoom=200");
	
close(ratio +"(red)");
close(ratio +" (green)");
close(ratio +" (blue)");
