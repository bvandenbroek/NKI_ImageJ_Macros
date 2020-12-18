lowerPercentile = 0.2;
upperPercentile = 0.8;
getRawStatistics(nPixels, mean, min, max, std, histogram);

lowerTotal = 0;
upperTotal = nPixels;
i=0;
while (lowerTotal < nPixels*lowerPercentile) {
	lowerTotal += histogram[i];
	//print("lower: "+lowerTotal);
	i++;
}
j=histogram.length-1;
while (upperTotal > nPixels*upperPercentile) {
	upperTotal -= histogram[j];
	//print("upper: "+upperTotal);
	j--;
}
setThreshold(i,j);