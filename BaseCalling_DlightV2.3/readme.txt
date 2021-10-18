Copyright (c) 2021, Xiaonan Fu, Biochemistry department,
University of Washington, Seattle WA98103, USA
	
Dlight User Guide
Dlight is a Matlab script to extract sub-pixel intensity for DNA sequencing base calling. Please note that the generation of final fastq file requires 3DEC.

Please contact Xiaonan Fu (xnfu@uw.edu) for any problems, bugs or suggestions.
	
	Pre-requisite:
	Ubuntu 14.01 LTS	
	Matlabv2019 or new version
	3Dec

	To fully run the pipeline:
	1, set the parameters in DlightP2_3.m:
	    totalCycle, the total sequencing cycles included in your folder. start from 1.
	    imageSize, the pixel width and height of collected sequencing images.
	    Fov, select folder name.
	    Image_dir, the main folder name.
            nthread, multiple thread running.

	2. The other parameters were optimalized for current platform, which might need to be adjusted if necessary.

	3. Running the script in Matlab will generate "intensity_file.txt"
	python3 AYB-style_cif_pos_writer.py ./intensity_file.txt 4 1
	~/biosoft/3dec-master/bin/3Dec -L -o ./output -q s_4_1

	4. A fastq file of s_4_1.fastq will be produced.
			
License
	3Dec is subject to  Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
	International Public License. A copy of the licence is attached with the software. You can 
	also obtain one at http://creativecommons.org/licenses/by-nc-sa/4.0/.
	
	Please notice that the source codes in the "include" folder are subject to different licences 
	such as MPL, MIT or BSD and the author of 3Dec does not have their copyright. 
	Licences for them can be found within or along with their files.
	
	
