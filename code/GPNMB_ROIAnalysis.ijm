/* ImageJ macro using .tif stacks containing 4 channels from MacNets plates (i.e. stained cell images)
 * >>> uses ROIs from GPNMB_Segmentation mecro and gets intensity averages, other measurements from them per chosen channel
 * WARNING: setBatchMode(true) seems to break this macro...
 * written by Lisa K. Polzer
 */

//set up for macro: defines variables, generates necessary components, etc.
parentDir = getDirectory("Choose parent directory.");
inputDirectory = getDirectory("Choose input directory.");
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec); //Note that ‘month’ and ‘dayOfWeek’ are zero-based indexes.
dateAndTime = "Date_"+year+"-"+month+"-"+dayOfMonth+"_Time_" +hour+"-"+minute+"-"+second;
fileList = getFileList(inputDirectory);
chosen_channel = getString("Chose channel: GPNMB:1 GFP:2 Iba1:3 MAP2:4", 1);
run("Input/Output...", "jpeg=85 gif=-1 file=.csv use_file copy_row save_column save_row"); //determines .tifs to be saved with MSB (big endian) by default
outputFile = inputDirectory + chosen_channel+"_allResults"+dateAndTime+".csv";  // Save results in the same folder

//set up of macro's main function
setBatchMode(false); //batchMode lets script run without calling images to GUI
for (i = 0; i < fileList.length; i++){
	if (matches(fileList[i], ".*filtered.*")) { //previously, to filter for dates: (matches(fileList[i], ".*20241206.*") && matches(fileList[i], ".*filtered.*"))
		main (inputDirectory, parentDir);
	}
}
setBatchMode(false);

function main (inputDirectory, parentDir) {
	//open Iba1/MAP2 filtered masks, open chosen channel raw images
	open(inputDirectory+fileList[i]);
	filename = getTitle();
	filename_base = substring(filename,0,lengthOf(filename)-18); 
	filename_chosen = filename_base+chosen_channel+".TIF";
	print(filename_base);
	print(filename_chosen);
	open(parentDir+filename_chosen);
	
	//get ROIs from specific channel-based ROIs
	selectImage(filename);
	setAutoThreshold("Default no-reset");
	run("Analyze Particles...", "size=0-Infinity pixel circularity=0-1 show=Nothing clear include add"); //get Iba1/MAP2-based ROIs via thresholding mask, note: edges included 
	
	//select chosen channel, add ROIs (according to if/else) to ROI manager
	selectImage(filename_chosen);
	run("Set Measurements...", "area mean min shape redirect=None decimal=3");
	roiManager("measure");
	
	// Save results table in shared and individual files
	File.append(filename_base, outputFile);
	ROI_number = Table.size();
	
	for (r = 0; r < ROI_number; r++) {
	    roiIndex = r + 1;  // ROI index (1-based)
	    ROI_area = Table.get("Area", r);  // Get area
	    ROI_mean = Table.get("Mean", r);  // Get mean intensity
	    ROI_min  = Table.get("Min", r);  // Get min intensity
	    ROI_max  = Table.get("Max", r);  // Get max intensity
	    ROI_circ = Table.get("Circ.", r);  // Get circularity
	    File.append(roiIndex + "," + ROI_mean + "," + ROI_min + "," + ROI_max + "," + ROI_area + "," + ROI_circ, outputFile);
	}
	
	//waitForUser;
	results_path = inputDirectory + File.separator + filename_base + chosen_channel + "_Results"+dateAndTime+".csv";
	Table.save(results_path);
	
	close("*");
}

print("\nFinished running macro 'GPNMB_ROIAnalysis'!\n");
