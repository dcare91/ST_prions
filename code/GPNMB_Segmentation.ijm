/* ImageJ macro using .tif stacks containing four channels from MacNets plates (i.e. stained cell images),
 * creates cell ROIs and stores them as masks for use in the second macro, GPNMB_ROIAnalysis.
 * After usage, check detection quality (e.g. via checking resulting thumbnails), if unsatisfactory, adjust setMinAndMax first.
 * Tested with ImageJ/Fiji version 2.16.0/1.54p.
 * Written by Lisa K. Polzer (0009-0007-5279-0883) for the GPNMB project.
 */

//set up for macro: define paths, path variables, ImageJ settings
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec); //‘month’ and ‘dayOfWeek’ are zero-based indices
dateAndTime = "Date_"+year+"-"+month+"-"+dayOfMonth+"_Time_" +hour+"-"+minute+"-"+second;
inputDirectory = getDirectory("Choose input directory.");  
fileList = getFileList(inputDirectory);
run("Input/Output...", "jpeg=85 gif=-1 file=.csv use_file copy_row save_column save_row"); //configure I/O settings: TIFFs are saved with MSB (big endian) by default

//dialog for other parameter variables
channels = newArray("1", "2", "3", "4");
Dialog.create("Segmentation Parameters");
Dialog.addChoice("Channel for ROI segmentation:", channels, channels[0]); //zero-based index
Dialog.addNumber("Maximum for setMinAndMax:",400); //GPNMB project: 400 for Iba1, 10000 for Map2
Dialog.addNumber("Foreground smoothing sigma:",2); //GPNMB: 2
Dialog.show();
chosenChannel = Dialog.getChoice();
channelNumber = parseInt(chosenChannel); //convert to integer
maxNum = Dialog.getNumber();
sigma = Dialog.getNumber();
outputDirectory = inputDirectory+"GPNMB_Segmentation"+dateAndTime+"_ch"+channelNumber+"_max"+maxNum+File.separator; //define output directory path
File.makeDirectory(outputDirectory); //generate output directory

//set up of macro's main function
setBatchMode(false); //batchMode(true) lets the script run without displaying images in the GUI
for (i = 0; i < fileList.length; i++){
	if  (matches(fileList[i], ".*w"+chosenChannel+".*")) { //limit the function to channel chosen via dialog
		main (inputDirectory, outputDirectory);
	}
}
setBatchMode(false);

function main (inputDirectory, outputDirectory) {
	//open .tif of chosen channel
	open(inputDirectory+fileList[i]);
	filename = getTitle();
	print("Segmenting "+filename+"...");
	
	//foreground smoothing on duplicate to remove Gaussian noise
	run("Duplicate...", "title=Smoothing duplicate");
	selectWindow("Smoothing");
	run("Gaussian Blur...", "sigma="+sigma);
	
	//BG subtraction
	run("Subtract Background...", "rolling=20 sliding"); //adjust as needed
	
	//threshold to include cell candidates, create mask with cells
	setMinAndMax(0, maxNum); 
	setAutoThreshold("Default dark no-reset");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Default background=Default black");
		
	//store segmented cell masks
	segmented_path = outputDirectory + filename + "_segmented.tif" ;
	saveAs("TIFF", segmented_path);
		
	//filter via particle size and circularity
		// Run Analyze Particles with different settings depending on channel
		if (channelNumber == 3) { //filtering for Iba1 channel, adjust as needed
		    run("Analyze Particles...", "size=350-Infinity pixel circularity=0.00-0.50 show=Masks clear include");
		} else if (channelNumber == 2) { //filtering for Map2 channel, adjust as needed
		    run("Analyze Particles...", "size=1000-Infinity pixel circularity=0.00-0.70 show=Masks clear ");
		} else {
		    print("No ROI filtering defined for channel " + channelNumber);
		}	
		
	//convert segmented image bit-depth	(from 8 to 16 bits)
	setOption("ScaleConversions", true); 
	run("16-bit"); 
		
	//store filtered cell masks
	filtered_path = outputDirectory + filename + "_" +  maxNum + "_" + sigma + "_filtered.tif";
	saveAs("TIFF", filtered_path);
	//for troubleshooting: check on segmentation & filtering here
	
	close("*");
}
print("\nFinished running macro 'GPNMB_Segmentation'!\n");
