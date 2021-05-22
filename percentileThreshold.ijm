args = getArgument();
if(args != "") percentile = parseFloat(args);
else percentile = 0.1;

getRawStatistics(nPixels, mean, min, max, std, histogram);
total = 0;
bin=0;
while (total < nPixels*percentile) {
	total += histogram[bin];
	bin++;
} 
setThreshold(0,bin-1);
