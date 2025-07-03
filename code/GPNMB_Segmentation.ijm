/* ImageJ macro using .tif stacks containing 4 channels from MacNets plates (i.e. stained cell images)
 * >>> creates cell ROIs and stores them as masks to then use in macros GPNMB_ROIAnalysis
 * After usage, run macro, check if detection is good (e.g. via checking resulting thumbnails), 
 * if not, check setMinAndMax first.
 * Written by Lisa K. Polzer for the GPNMB project.
 */

//set up for macro: defines variables, generates necessary components, etc.
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec); //Note that ‘month’ and ‘dayOfWeek’ are zero-based indices.
dateAndTime = "Date_"+year+"-"+month+"-"+dayOfMonth+"_Time_" +hour+"-"+minute+"-"+second;
inputDirectory = getDirectory("Choose input directory."); 
fileList = getFileList(inputDirectory);
run("Input/Output...", "jpeg=85 gif=-1 file=.csv use_file copy_row save_column save_row"); //determines .tifs to be saved with MSB (big endian) by default

//dialog for variable other parameters
channels = newArray("1", "2", "3", "4"); //"GPNMB: 1", "GFP: 2", "IBA1: 3", "MAP2: 4"
Dialog.create("z-Cropping");
Dialog.addChoice("Channel:", channels, channels[1]); //numbers here start from 0
Dialog.addNumber("Maximum for setMinAndMax:",400); //GPNMB projetc: 400 for Iba1, 10000 for Map2
Dialog.addNumber("Foreground smoothing sigma:",2);
Dialog.show();
chosenChannel = Dialog.getChoice();
channelNumber = parseInt(chosenChannel); // Convert to integer
maxNum = Dialog.getNumber();
smallSigma = Dialog.getNumber();
outputDirectory = inputDirectory+"GPNMB_Segmentation"+dateAndTime+channelNumber+maxNum+File.separator; //defines output folder path
File.makeDirectory(outputDirectory); //generatres output folder

//set up of macro's main function
setBatchMode(false); //batchMode lets script run without calling images to GUI
for (i = 0; i < fileList.length; i++){
	if  (matches(fileList[i], ".*w"+chosenChannel+".*")) { //limits the function to tif files  //to limit for date also, for prev. channel numbers: (matches(fileList[i], ".*20241206.*") && matches(fileList[i], ".*ch0"+threshChannel_number+".*"))
		main (inputDirectory, outputDirectory);
	}
}
setBatchMode(true);

function main (inputDirectory, outputDirectory) {
	//opens, smoothes .tif stack for thresholding
	open(inputDirectory+fileList[i]);
	filename = getTitle();
	
	//foreground smoothing on duplicate to remove Gaussian noise
	run("Duplicate...", "title=smallBlur duplicate");
	smallBlur = getTitle();
	selectWindow("smallBlur");
	run("Gaussian Blur...", "sigma="+smallSigma);
	
	//BG subtraction
	run("Subtract Background...", "rolling=20 sliding");
	
	//threshold to include cell candidates, create mask with cells, filter based on specific metrics
	setMinAndMax(0, maxNum); 
	setAutoThreshold("Default dark no-reset");
	setOption("BlackBackground", true);
	run("Convert to Mask", "method=Default background=Default black");
		
	//store segmented cells
	segmented_path = outputDirectory + filename + "_segmented.nrrd" ;
	saveAs("TIFF", segmented_path);
		
	//filter via particle size and circularity, change bit-depth
		// Run Analyze Particles with different settings depending on channel
		if (channelNumber == 3) {
		    run("Analyze Particles...", "size=350-Infinity pixel circularity=0.00-0.50 show=Masks clear include");
		} else if (channelNumber == 4) {
		    run("Analyze Particles...", "size=1000-Infinity pixel circularity=0.00-0.70 show=Masks clear ");
		} else {
		    print("No Analyze Particles action defined for channel " + channelNumber);
		}	

	setOption("ScaleConversions", true); //converts segmented image from 8 to 16 bits to match other files' types
	run("16-bit"); 
		
	//store segmented cells
	filtered_path = outputDirectory + filename + "_filtered.nrrd" ;
	saveAs("TIFF", filtered_path);
	//waitForUser; //for troubleshooting, to check on segmentation & filtering
	
	close("*");
}
print("\nFinished running macro 'GPNMB_Segmentation'!\n");
