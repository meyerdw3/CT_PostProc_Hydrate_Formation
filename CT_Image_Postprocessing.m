%---------------------------------------------------------------------------------------------------------
%INTRODUCTION
%Code was written by Dylan Whitney Meyer
%The University of Texas at Austin, Jackson School of Geosciences
%Last Updated: August 1, 2017

%This code performs the semi-autonomous full procedure for post-processing CT images 
%collected from the modified CT scanner in the Petroleum and Geosystems Engineering 
%Department at the University of Texas at Austin. The code begins by extracting the 
%attenuation data from the raw data file output by the CT scanner. The data is then masked, 
%filtered, and corrected to produce a reliable dataset for further analysis. The porosity, 
%bulk density, and formation front parameters are calculated from the finalized dataset. 
%The code has the ability to produce plots to visualize the data produced during post-
%processing, but that option can be turned off if desired. At no point during the process is
%any data overwritten. All iterations of the data from raw to corrected are saved in 
%separate folders to ensure continuity and for error analysis.

%Before running this script, the user is required to initialize the following folder 
%structure:

%   ./HVTXXXX - Folder for experiment being processed
%   ./HVTXXXX/scan_info/ - Folder with information on the scan parameters
%   ./HVTXXXX/pile_files/ - Folder with PILE files

%The scan_info folder should contain a .csv document that is a list of the elapsed time for 
%each scan taken during the experiment. The pile_files folder should contain the raw PILE 
%files output by the CT scanner. The user is also required to provide additional information
%about the scan parameters below.

%After the post-processing completes, this script will have produced and populated the 
%following folder structure:

%   ./HVTXXXX/mat_files - Folder with MAT files of the raw and processed CT data
%   ./HVTXXXX/mat_files/HVTXX_YY/ - Folder with MAT data from one scan
%   ./HVTXXXX/mat_files/HVTXX_YY/raw - Folder with raw MAT data
%   ./HVTXXXX/mat_files/HVTXX_YY/masked - Folder with masked MAT data
%   ./HVTXXXX/mat_files/HVTXX_YY/filtered - Folder with filtered MAT data
%   ./HVTXXXX/mat_files/HVTXX_YY/corrected - Folder with corrected MAT data
%   ./HVTXXXX/mat_files/HVTXX_YY/density - Folder with density data
%   ./HVTXXXX/mat_files/HVTXX_YY/front - Folder with front data

%   ./HVTXXXX/tif_files - Folder with TIFF images of every slice from every scan
%   ./HVTXXXX/tif_files/HVTXX_YY/ - Folder with TIFF images from one scan

%   ./HVTXXXX/results - Folder with finalized results from the post-processing
%   ./HVTXXXX/results/front_data - Folder with front parameter results
%   ./HVTXXXX/results/ortho_data - Folder with orthogonal density results
%   ./HVTXXXX/results/ortho_plots - Folder with orthogonal density plots
%   ./HVTXXXX/results/porosity - Folder with porosity results

%Several of the subroutines used in this code were written by Kehua You at the University 
%of Texas at Austin, Jackson School of Geosciences:

%   brine_density.m
%   gas_density.m
%   methane_solubility.m
%   calc_Vr.m
%   duan.m

%Several of the subroutines used in this code were obtained with permission from David 
%DiCarlo at the  University of Texas at Austin, Petroleum and Geosystem Engineering 
%Department:

%   circle_eq.m
%   circlepoints.m
%   imcircle.m
%   Pile_to_tif_mat.m

%The following subroutine was download from the MathWorks File Exchange system in order to
%create contoured colorbars when plotting the final data:

%   cbarf.m

%---------------------------------------------------------------------------------------------------------
%INITIALIZATION

fprintf('--- Initializing Program --- ')

%Add the path to the script subroutines to the search path
addpath('./SUBROUTINES/')

%---------------------------------------------------------------------------------------------------------
%EXPERIMENTAL PARAMETER DEFINITIONS
%These values change between experiments and should be updated to reflect the 
%experiment currently being processed

%SUBROUTINES:
%   methane_solubility.m
%   brine_density.m
%   gas_density.m

P = 12.24; %Experimental pore pressure, [MPa]
T = 1.05; %Experimental temperature, [degC]
phi = 0.405; %Gravimetric sample porosity, []
L = 10.2; %Sample length, [cm]
r = 2.54; %Sample radius, [cm]
cl = 0.07; %Initial brine salinity, [wt% NaCl]
m = (methane_solubility(P, T, ((cl / 0.05844) / ...
                         (1 - cl))) * 0.016); %Methane solubility, [wt%]
rhow = (brine_density(P, T, ((cl / 0.05844) / ...
                         (1 - cl)), (m / 0.016)) / 1000); %Initial brine density, [g/cm3]
rhog = (gas_density(P, T) / 1000); %Methane density, [g/cm3]
rhos = 2.65; %Grain density, [g/cm3]
Vtot = (pi() * (r ^ 2) * L); %Total sample volume, [mL]
Vpore = (Vtot * phi); %Sample pore volume, [mL]

%---------------------------------------------------------------------------------------------------------
%SCAN PARAMETER DEFINITIONS
%These values defined the parameters used during each scan and construct the data storage 
%structure for the post-processing procedure

%Define the names of the experiment being processed, the saturated and unsaturated scans, 
%and the scan prefix and the location of the scan information and PILE files
exp_name = 'HVT0023'; %Name of experiment being processed
wet_scan = 'HVT23_Wet_130'; %Name of the saturated scan
dry_scan = 'HVT23_Dry_130'; %Name of the unsaturated scan
scan_pref = 'HVT23_'; %Scan name prefix
scan_sufx = '_130'; %Scan name suffix
pile_files = strcat('./', exp_name, '/pile_files/'); %Location of PILE files
scan_info = strcat('./', exp_name, '/scan_info/'); %Location of the scan information

%Define the scan parameters
n_scans = 36; %Total number of scans to be processed
n_slices = 35; %Total number of slices per scan
slice_space = 0.3; %Distance between slices, [cm]
slice_dist = [L:(-1 * slice_space):0]; %Distance of each slice from the outlet, [cm]

%Define the directory names for the post-processing output and create the folders
results = strcat('./', exp_name, '/results/'); %Results folder location
mat_files = strcat('./', exp_name, '/mat_files/'); %Matlab files folder location
tif_files = strcat('./', exp_name, '/tif_files/'); %TIFF files folder location
mkdir(results);
mkdir(mat_files);
mkdir(tif_files);

%Define a list of scans that should be skipped for any reason (incorrect  data collection 
%or recording, forgot to export it, corrupt PILE file). Hopefully this list is empty, but 
%mistakes happen during scanning and this provides flexibility in the post-processing 
%procedure.
%EXAMPLE USAGE: skip_list = [1 2 3];
skip_list = [1 26];

%Load the scan times and the experimental pump volume and mass balance results
scan_times = csvread(strcat(scan_info, 'scan_times.csv'));

%---------------------------------------------------------------------------------------------------------
%DEFINE ESSENTIAL FILE NAMES

%Define directories of the saturated and unsaturated scan PILE files
wet_pile_file = strcat(pile_files, 'pile_', wet_scan);
dry_pile_file = strcat(pile_files, 'pile_', dry_scan);

%Initialize cell arrays to hold the directories and names of each PILE file and the names 
%of each scan
file_names = cell(n_scans, 1);
pile_names = cell(n_scans, 1);
scan_names = cell(n_scans, 1);

%Populate the cell arrays with the appropriate data
for i = 1:n_scans
    %Skip the scan if it's found in the skip list
    if (sum(i == skip_list) == 1)
        continue
    end

    %If the scan number is less than 10, add the buffer '0' to the names
    if (i < 10)
        %PILE file directory
        file_names{i} = [strcat(pile_files, 'pile_', scan_pref, '0', num2str(i), ...
                                scan_sufx)];
        
        %PILE file name
        pile_names{i} = [strcat('pile_', scan_pref, '0', num2str(i), scan_sufx)];
        
        %Scan name
        scan_names{i} = [strcat(scan_pref, '0', num2str(i), scan_sufx)];
    end

    %If the scan number is greater than 10, no buffer '0' is required
    if (i >= 10)
        %PILE file directory
        file_names{i} = [strcat(pile_files, 'pile_', scan_pref, num2str(i), scan_sufx)];
        
        %PILE file name
        pile_names{i} = [strcat('pile_', scan_pref, num2str(i), scan_sufx)];
        
        %Scan name
        scan_names{i} = [strcat(scan_pref, num2str(i), scan_sufx)];
    end
end

%---------------------------------------------------------------------------------------------------------
%POST-PROCESSING CHECKPOINT DEFINITIIONS
%Adding checkpoints to the procedure improves functionality by keeping track of how far 
%the post-processing has progressed. At the beginning all of the checkpoints are initialized
%to 'n', indicated that none of them have been completed. As the procedure completes 
%sections, it will update these checkpoints to 'y', indicating the task has been completed. 
%If there is an issue during post-processing, the user can correct the error and then 
%restart the process from where the error occurred, rather than performing it from the 
%beginning. Additionally, if the user wants to only perform one portion of the process, 
%they can change the default settings below to indicate which piece of the code they wish 
%to execute.

%Each time the script is run, it will ask the user if they want to reinitialize the
%checkpoints to the defaults. If the post-processing is being run for the first time or the
%user only wants to run a certain portion, they should select 'Yes'. If they want to restart
%the procedure from where they left off, they should select 'No'.
reply = questdlg('Would you like to reinitialize the post-processing checkpoints?',...
                 'Reinitialize Checkpoints?', 'Yes', 'No', 'Yes');

if strcmp(reply, 'Yes')    
    %CHECKPOINT 1: Extract the scan data from the PILE files and save them as TIFF 
    %and MAT files
    data_extract = 'n';
    
    %CHECKPOINT 2: Select the position of the sample in each slice
    pos_pick = 'n';
    
    %CHECKPOINT 3: Mask every slice in every scan and center the sample in the matrix
    scan_mask = 'n'; 
    
    %CHECKPOINT 4: Apply the 10-pixel median filter to every slice in every scan
    scan_filt = 'n';
    
    %CHECKPOINT 5: Determine and apply the tube heat correction to every slice in every 
    %scan
    heat_correct = 'n';
    
    %CHECKPOINT 6: Calculate the 3D sample porosity from the saturated and 
    %unsaturated scans
    porosity = 'n';
    
    %CHECKPOINT 7: Calculate the 3D sample bulk density for every scan
    density = 'n';
    
    %CHECKPOINT 8: Determine the essential front parameters for every scan
    v_aff = 'n';
    
    %CHECKPOINT 9: Record the 2D distribution of the porosity and bulk density 
    %orthogonal to the flow direction
    ortho_data = 'n';
end

%Define how the data will the plotted and recorded ('n' means that particular type
%of plot will not be produced)
axial_plots = 'n';
ortho_plots = 'n';
colorbars = 'n';

fprintf('DONE\n')

%---------------------------------------------------------------------------------------------------------
%CHECKPOINT 1: DATA EXTRACTION
%This section of code loads each PILE file exported from the CT scanner and extracts the 
%CT attenuation data in the form of TIFF and MAT files.

%SUBROUTINES:
%   Pile_to_tif_mat.m

fprintf('--- Extracting Data From PILE Files --- \n')

while (data_extract == 'n')
    %Extract the data from the saturated and unsaturated PILE files
    Pile_to_tif_mat(exp_name, wet_pile_file, wet_scan, n_slices);
    Pile_to_tif_mat(exp_name, dry_pile_file, dry_scan, n_slices);
    
    %Iterate through all the successful scans and extract the data from the PILE files
    for i = 1:n_scans
        %Skip the scan if it's found in the skip list
        if (sum(i == skip_list) == 1)
            continue
        end

        %Run the current scan through the conversion script
        Pile_to_tif_mat(exp_name, strjoin(file_names(i, 1)), ...
                                  strjoin(scan_names(i, 1)), n_slices);
        
        %Display the completion progress of the current procedure
        fprintf('Percent Complete: %d\n', round(((i - 1) / n_scans) * 100))
    end
    
    %Display the process as complete
    fprintf('Percent Complete: %d\n', round(100))
    
    %Close any open figures
    close all

    %Update the checkpoint status
    data_extract = 'y';
end

fprintf('DONE\n')

%---------------------------------------------------------------------------------------------------------
%CHECKPOINT 2: POSITION PICKING
%This section of code was modified from a script written by David DiCarlo that allows 
%the user to pick the sample radius and then the centerpoint of the sample in every slice 
%in the saturated scan. This information is used later to remove all the data exterior to 
%the sample and the center the sample within the data matrix.

%SUBROUTINES:
%   circle_eq.m
%   circlepoints.m

fprintf('--- Picking Sample Positions --- ')

while(pos_pick == 'n')
    %STEP 1: Determine the radius of the sample
    %The program displays a slice from the middle of the sample, to reduce any artifacts 
    %produced by the endcaps, and then the user pick three points on the edge of the 
    %sample to define the circle. The user is allowed to correct the center position 
    %and radius if necessary.
    
    %Define the name of the slice in the middle of the sample and load the slice
    name = strcat(mat_files, wet_scan, '/raw/slice', num2str(floor(n_slices / 2)), '.mat');
    load(name);

    %Open a new figure window and allow it to plot multiple objects
    figure
    hold on

    %This next section displays the image and allows the user to pick three points on 
    %the perimeter. This while loop runs as long as the array, pts, has less than three 
    %points in it.
    pts = [];
    while (size(pts,1) ~= 3)
        clear x_r y_r

        %Display the image. The numbers in brackets control the contrast and brightness.
        imshow(slice,[800 2200]);
        title('First Slice');

        %Prompt the user to pick a point and record it
        questdlg('Choose three points on the edge of the circle and hit ENTER',...
                 'Circle Selection','OK', 'Cancel', 'OK');
        [x_r(:,1),y_r(:,1),P] = impixel;
        pts = [x_r y_r];

        %When pts contains three points, this statements runs and determines the radius 
        %using the script, circle_eq.
        if (size(pts, 1) == 3)
            
            %Record the center position (x, y) and radius into three separate variables
            [circ1,circ2,circ3] = circle_eq(pts(1,:),pts(2,:),pts(3,:));

            %Transfer and round the center position and radius into an array, circ
            circ(1) = round(circ1);
            circ(2) = round(circ2);
            circ(3) = round(circ3);

            %Transfer the center position into a vector, I_cen
            I_cen = [circ(1) circ(2)];

            %Transfer the radius into a variable, rad
            rad = circ(3);

            %Determines the circle pixels
            circ_pixels = circlepoints(I_cen(1), I_cen(2), rad);
            hold on

            %Display the circle on the image, so that the user can tell if it is in the 
            %right place
            plot(circ_pixels(:,1), circ_pixels(:,2), 'g.', 'markersize', 2);
            
        else
            
            %If the user picks more or less than three points, the user is told to retry
            printf('Need three and only three points. Retry.\n');
            
            %In order to retry the array holding the selected points needs to be cleared 
            %and the figure needs to be closed
            pts = [];
            close all;
            
        end
    end
    
    %STEP 2: Check the sample radius and allow for adjustments
    %Ask the user if the center position and radius are correct
    reply = questdlg('Circle Check', 'Circle OK?', 'Yes', 'No', 'No');

    %If these values are not correct the user is asked to adjust them. A dialog box will 
    %open and the user will be able to enter new values for the x-y position and radius. 
    %The drawn circle will be adjusted and the program will ask the user again if the 
    %circle is correct.
    while strcmp(reply, 'No')
        
        %Display the dialog box that prompts the user for changes
        prompt = {'Adjust x-center?', 'Adjust y-center?', 'Adjust radius?'};
        dlg_title = 'Adjust circle?';
        num_lines = 1;
        def = {num2str(I_cen(1)), num2str(I_cen(2)), num2str(rad)};
        answer = inputdlg(prompt, dlg_title, num_lines, def);

        %Overwrite the old values with the new ones
        I_cen(1) = str2num(char(answer(1)));
        I_cen(2) = str2num(char(answer(2)));
        rad = str2num(char(answer(3)));
        
        %Close the old figure before plotting the new one
        close all

        %Redisplay the figure without the circle
        imshow(slice, [800 2200]);
        hold on

        %Calculate the points that make up the circle and display them onto the image.
        circ_pixels = circlepoints(I_cen(1), I_cen(2), rad);
        plot(circ_pixels(:,1), circ_pixels(:,2), 'g.', 'markersize', 2);

        %Ask the user if the new circle parameters are OK, if not then it returns to 
        %the beginning of the while loop to ask the user to adjust the parameters again.
        reply = questdlg('Circle Check','Circle OK?','Yes','No','No');
    end
    
    %Close all open figures before next step
    close all
    
    %STEP 3: Pick the center positions for the sample in every slice
    %Iterate through the slice and allow the user to correct the sample position. The 
    %sample radius is kept constant.
    for i = 1:n_slices
        %Define the file name of the current slice and load it
        name = strcat(mat_files, wet_scan, '/raw/slice', num2str(i), '.mat');
        load(name);

        %Display the image and label it with the slice number
        title_slice = strcat('slice', num2str(i));
        imshow(slice,[750 1500]);
        title(title_slice);
        hold on

        %Calculate the points that make up the circle from the previous picked center 
        %point and display them onto the image.
        plot(circ_pixels(:,1),circ_pixels(:,2),'g.','markersize',2);

        %Ask the user if the center position is correct
        reply = questdlg(' Check circle position', 'Circle OK?', 'Yes', 'No', 'No');

        %If the centerpoint is not correct, the user is asked to adjust the centerpoint.
        while strcmp(reply,'No')
            
            %Set up the dialog box and and prompt the user for changes
            prompt = {'Adjust x-center?', 'Adjust y-center?'};
            dlg_title = 'Adjust circle?';
            num_lines = 1;
            def = {num2str(I_cen(1)), num2str(I_cen(2))};
            answer = inputdlg(prompt, dlg_title, num_lines, def);

            %Overwrite the old values with the new ones
            I_cen(1) = str2num(char(answer(1)));
            I_cen(2) = str2num(char(answer(2)));
            
            %Close the old figure before plotting the new one
            close all

            %Redisplay the image without the circle
            imshow(slice, [750 1500]);
            title(title_slice);
            hold on

            %Calculate the points that make up the circle and display them onto the image
            circ_pixels = circlepoints(I_cen(1), I_cen(2), rad);
            plot(circ_pixels(:,1), circ_pixels(:,2), 'g.', 'markersize', 2);

            %Ask the user if the center position is correct. If not, then the program 
            %loops back to the beginning of the while statements and tries again.
            reply = questdlg('Circle Check', 'Circle OK?', 'Yes', 'No', 'No');
        end
        
        %Close all open figures before next step
        close all

        %Save the center position of each slice into an array
        I_cen_vec(i,:) = I_cen;

        %Save the center position and radius into a mat file. As one is an array and the 
        %other is a single value, the information will be saved into a matlab struct with 
        %two parts (rad and I_cen_vec) that can be separately called. Look below for how 
        %to reference either piece of information.
        outfile_name = strcat(strcat(scan_info, 'position', '.mat'));
        save(outfile_name, 'I_cen_vec', 'rad');
    end

    %Update the checkpoint status
    pos_pick = 'y';
end

fprintf('DONE\n')

%Load sample location data. It's important to do this externally to the while statement 
%above, because it assures the position data is loaded into the program for later use, 
%regardless of whether or not that while statement was executed.
pos = load(strcat(scan_info, 'position.mat'));
rad = pos.rad;
centers = pos.I_cen_vec;

%---------------------------------------------------------------------------------------------------------
%CHECKPOINT 3: DATA MASKING
%This section crops the raw data down to a square matrix with the sample circumscribed 
%in the center and then set all the attenuation values outside the sample to 0.

%SUBROUTINES:
%   imcircle.m

fprintf('--- Applying Slice Masks --- \n')

%Initialize the sample mask. This command will create a 2D matrix with dimensions of 
%((2 * rad) + 1) with a circle circumscribed within the matrix boundaries. Inside the 
%circle, the value of each cell is equal to 1. Outside the circle, the value of each cell 
%is equal to 0.
mask = int16(imcircle(2 * rad + 1));

%Sum up the elements in the sample mask to calculate the total number of pixels within 
%the sample.
mask_sum = sum(sum(mask));

while (scan_mask == 'n')

    %Create new directories to hold the masked MAT files of the WET and DRY scans.
    wet_path = strcat(mat_files, wet_scan, '/masked/');
    dry_path = strcat(mat_files, dry_scan, '/masked/');
    mkdir(wet_path);
    mkdir(dry_path);
    
    %Iterate through the slices in the unsaturated and saturated scans and apply the mask
    for i = 1:n_slices
        %Load a slice from both the WET and DRY scans
        slice_wet = load(strcat(mat_files, wet_scan, '/raw/slice', num2str(i), '.mat'));
        slice_dry = load(strcat(mat_files, dry_scan, '/raw/slice', num2str(i), '.mat'));

        %Convert the slices from STRUCT to MAT
        slice_wet = slice_wet.slice;
        slice_dry = slice_dry.slice;

        %Determine the x and y ranges within the slices that are within the sample.
        x_range = (centers(i, 1) - rad):(centers(i, 1) + rad);
        y_range = (centers(i, 2) - rad):(centers(i, 2) + rad);

        %Crop the slice down the x and y ranges specified above and then multiple it by 
        %the mask. This process sets all the data outside the sample equal to 0.
        slice_wet = (slice_wet(y_range, x_range) .* mask);
        slice_dry = (slice_dry(y_range, x_range) .* mask);

        %Create the file name for the modified WET and DRY slices and save them to their 
        %folders
        outfile_name = strcat(wet_path, 'slice', num2str(i), '.mat');
        save(outfile_name, 'slice_wet');

        outfile_name = strcat(dry_path, 'slice', num2str(i), '.mat');
        save(outfile_name, 'slice_dry');

        %Clear the WET and DRY slices
        clear slice*
    end
    
    %Iterate through each of the experimental scans
    for j = 1:n_scans
        %Skip the scan if it's found in the skip list
        if (sum(j == skip_list) == 1)
            continue
        end

        %Define and create the directory that will hold the masked slices
        path = strcat(mat_files, strjoin(scan_names(j, 1)), '/masked/');
        mkdir(path);

        %Iterate through the slices in the current experimental scan and apply the mask
        for i = 1:n_slices
            %Load the current slice
            load(strcat(mat_files, strjoin(scan_names(j, 1)), '/raw/slice', ...
                                                  num2str(i), '.mat'));

            %Determine the x and y ranges within the slices that are within the sample
            x_range = (centers(i, 1) - rad):(centers(i, 1) + rad);
            y_range = (centers(i, 2) - rad):(centers(i, 2) + rad);

            %Crop the slice down the x and y ranges specified above and then multiple it 
            %by the mask. This process sets all the data outside the sample equal to 0.
            slice = (slice(y_range, x_range) .* mask);

            %Create the file name for the modified slice and save it
            outfile_name = strcat(path, 'slice', num2str(i), '.mat');
            save(outfile_name, 'slice');

            %Clear the current slice
            clear slice*
        end
        
        %Display the completion progress of the current procedure
        fprintf('Percent Complete: %d\n', round(((j - 1) / n_scans) * 100))
    end
 
    %Display this procedure as complete
    fprintf('Percent Complete: %d\n', round(100))
    
    %Update this checkpoint status
    scan_mask = 'y';
end

fprintf('DONE\n')

%---------------------------------------------------------------------------------------------------------
%CHECKPOINT 4: DATA FILTERING
%This section applies a median filter to all the masked slices. This process reduces the
%instrument error and helps remove anomalous attenuation values, while retaining real
%heterogeneous structures

%Define the desired median filter size. For our system, a 10-pixel filter was determined 
%to be optimal. My thesis has additional information on how this was determined.
filt = 10;

fprintf('--- Filtering Wet/Dry/Experimental Scans --- \n')

while (scan_filt == 'n')
    %Creates new directories to hold the masked mat files of the WET and DRY scans
    wet_path = strcat(mat_files, wet_scan, '/filtered/');
    dry_path = strcat(mat_files, dry_scan, '/filtered/');
    mkdir(wet_path);
    mkdir(dry_path);
    
    %Iterate through the slices in the saturated and unsaturated scans and apply the 
    %median filter
    for i = 1:n_slices
        %Load a slice from both the WET and DRY scans
        load(strcat(mat_files, wet_scan, '/masked/slice', num2str(i), '.mat'));
        load(strcat(mat_files, dry_scan, '/masked/slice', num2str(i), '.mat'));

        %Apply the 2D median filter to each image to smooth out the image
        slice_wet = medfilt2(slice_wet, [filt filt]);
        slice_dry = medfilt2(slice_dry, [filt filt]);

        %Create the file name for the modified WET and DRY slices and save them to 
        %their folders
        outfile_name = strcat(wet_path, 'slice', num2str(i), '.mat');
        save(outfile_name, 'slice_wet');

        outfile_name = strcat(dry_path, 'slice', num2str(i), '.mat');
        save(outfile_name, 'slice_dry');

        %Clear the WET and DRY slices
        clear slice*
    end
        
    %Iterate through each of the experimental scans
    for j = 1:n_scans
        %Skip the scan if it's found in the skip list
        if (sum(j == skip_list) == 1)
            continue
        end

        %Define and create the directory that will hold the filtered slices
        path = strcat(mat_files, strjoin(scan_names(j, 1)), '/filtered/');
        mkdir(path);

        for i = 1:n_slices
            %Load the current slice
            load(strcat(mat_files, strjoin(scan_names(j, 1)), '/masked/slice', ...
                                              num2str(i), '.mat'));
            
            %Apply a 2D median filter to each image to smooth out the image
            slice = medfilt2(slice, [filt filt]);

            %Create the file name for the modified slice and save it
            outfile_name = strcat(path, 'slice', num2str(i), '.mat');
            save(outfile_name, 'slice');

            %Clear the current slice
            clear slice*
        end
        
        %Display the completion progress of the current procedure
        fprintf('Percent Complete: %d\n', round(((j - 1) / n_scans) * 100))
    end
 
    %Display this procedure as complete
    fprintf('Percent Complete: %d\n', round(100))
    
    %Update this checkpoint status
    scan_filt = 'y';
end

fprintf('DONE\n')

%---------------------------------------------------------------------------------------------------------
%CHECKPOINT 5: HEAT CORRECTION
%This section determine the appropriate correction value that needs to be applied to each
%slice in every scan in order to remove the effect of the CT tube heat on the attenuation

%SUBROUTINES:
%   imcircle.m

fprintf('--- Calculating/Applying CT Heat Correction Factor --- \n')

while (heat_correct == 'n')
    
    %STEP 1: Isolate the confining fluid attenuations
    %Adjust the internal and external diameters and centerpoint of the confining reservoir 
    %to remove data from all other materials
    
    %Define the slice to be used for this step
    test_slice = floor(n_slices / 2);
    
    %Take an initial guess at the internal (ID) and external (OD) diameters and the 
    %centerpoint
    OD = 200;
    ID = 140;
    conf_center = [258, 261];

    %Define the x, y range of the OD of the confining region in the scan from the current 
    %centerpoint
    OD_x_range = (conf_center(1, 1) - OD):(conf_center(1, 1) + OD);
    OD_y_range = (conf_center(1, 2) - OD):(conf_center(1, 2) + OD);

    %Define the x, y range of the ID of the confining region in the scan from the current 
    %centerpoint
    ID_x_range = (centers(test_slice, 1) - ID):(centers(test_slice, 1) + ID);
    ID_y_range = (centers(test_slice, 2) - ID):(centers(test_slice, 2) + ID);

    %Create a mask for the OD and ID
    mask_OD = int16(imcircle(2 * OD + 1));
    mask_ID = int16(imcircle(2 * ID + 1));

    %Define a blank confining mask of zeros to start with
    mask_conf = zeros(512, 512);

    %Add the OD mask over the correct range to change the value of the appropriate cells 
    %within that region to 1.
    mask_conf(OD_y_range, OD_x_range) = (mask_conf(OD_y_range, OD_x_range) + double(mask_OD));

    %Add the ID mask over the correct range to change the value of the appropriate cells 
    %within that region to 0. This will leave only the confining region with cells equal 
    %to 1.
    mask_conf(ID_y_range, ID_x_range) = (mask_conf(ID_y_range, ID_x_range) - double(mask_ID));

    %Load a slice from the middle of the WET scan
    load(strcat(mat_files, wet_scan, '/raw/slice', num2str(test_slice), '.mat'));

    %Apply the confining mask to the slice to isolate the confining fluid data
    slice = (double(slice) .* mask_conf);

    %Iterate through the slice and eliminate any value over 0 or under 200 to remove the 
    %thermistor and material at the edges of the confining fluid zone
    for j = 1:512
        for k = 1:512
            if (slice(j, k) > 0 | slice(j, k) < -200)
                slice(j, k) = 0;
            end
        end
    end

    %Plot the current slice to show the region that will be included in the heat correction 
    %calculation
    figure
    mesh(slice)
    view(2)
    
    %Ask the user if the selected region is OK
    reply = questdlg('Region Check', 'Region OK?', 'Yes', 'No', 'No');

    %If these values are not correct the user is asked to adjust them. A dialog box will 
    %open and the user will be able to enter new values for the ID, OD, and centerpoint.
    while strcmp(reply, 'No')
        
        %Display the dialog box that prompts the user for changes
        prompt = {'Adjust y-center?', 'Adjust x-center?', 'Adjust ID?', 'Adjust OD?'};
        dlg_title = 'Adjust Region?';
        num_lines = 1;
        def = {num2str(conf_centers(1)), num2str(conf_centers(2)), ID, OD};
        answer = inputdlg(prompt, dlg_title, num_lines, def);

        %Overwrite the old values with the new ones
        conf_centers(1) = str2num(char(answer(1)));
        conf_centers(2) = str2num(char(answer(2)));
        ID = str2num(char(answer(3)));
        OD = str2num(char(answer(4)));
        
        %Close the old figure before plotting the new one
        close all
        
        %Define the x, y range of the OD of the confining region in the scan from the 
        %current centerpoint
        OD_x_range = (conf_center(1, 1) - OD):(conf_center(1, 1) + OD);
        OD_y_range = (conf_center(1, 2) - OD):(conf_center(1, 2) + OD);

        %Define the x, y range of the ID of the confining region in the scan from the 
        %current centerpoint
        ID_x_range = (centers(test_slice, 1) - ID):(centers(test_slice, 1) + ID);
        ID_y_range = (centers(test_slice, 2) - ID):(centers(test_slice, 2) + ID);

        %Create a mask for the OD and ID
        mask_OD = int16(imcircle(2 * OD + 1));
        mask_ID = int16(imcircle(2 * ID + 1));

        %Define a blank confining mask of zeros to start with
        mask_conf = zeros(512, 512);

        %Add the OD mask over the correct range to change the value of the appropriate 
        %cells within that region to 1.
        mask_conf(OD_y_range, OD_x_range) = (mask_conf(OD_y_range, OD_x_range) + double(mask_OD));

        %Add the ID mask over the correct range to change the value of the appropriate 
        %cells within that region to 0. This will leave only the confining region with 
        %cells equal to 1.
        mask_conf(ID_y_range, ID_x_range) = (mask_conf(ID_y_range, ID_x_range) - double(mask_ID));

        %Load a slice from the middle of the WET scan
        load(strcat(mat_files, wet_scan, '/raw/slice', num2str(test_slice), '.mat'));

        %Apply the confining mask to the slice to isolate the confining fluid data
        slice = (double(slice) .* mask_conf);

        %Iterate through the slice and eliminate any value over 0 or under 200 to remove 
        %the thermistor and material at the edges of the confining fluid zone
        for j = 1:512
            for k = 1:512
                if (slice(j, k) > 0 | slice(j, k) < -200)
                    slice(j, k) = 0;
                end
            end
        end

        %Plot the current slice to show the region that will be included in the heat 
        %correction calculation
        figure
        mesh(slice)
        view(2)

        %Ask the user if the new circle parameters are OK, if not then it returns to the 
        %beginning of the while loop to ask the user to adjust the parameters again.
        reply = questdlg('Region Check', 'Region OK?', 'Yes', 'No', 'No');
    end
    
    %Close the figure
    close all
    
    %STEP 2: Determine heat correction
    %The correction value for each slice is determined using the change in the confining 
    %oil attenuation from the average attenuation
    
    %Initialize a matrix that will hold the average confining oil attenuations
    avg_conf_atten = zeros(n_slices, 1);
    
    %Iterate through the WET scan and determine the average attenuation of the confining 
    %fluid
    for i = 1:n_slices
        
        %Define the x, y range of the OD of the confining region in the scan from the 
        %current centerpoint
        OD_x_range = (conf_center(1, 1) - OD):(conf_center(1, 1) + OD);
        OD_y_range = (conf_center(1, 2) - OD):(conf_center(1, 2) + OD);

        %Define the x, y range of the ID of the confining region in the scan from the 
        %current centerpoint
        ID_x_range = (centers(i, 1) - ID):(centers(i, 1) + ID);
        ID_y_range = (centers(i, 2) - ID):(centers(i, 2) + ID);

        %Create a mask for the OD and ID
        mask_OD = int16(imcircle(2 * OD + 1));
        mask_ID = int16(imcircle(2 * ID + 1));

        %Define a blank confining mask of zeros to start with
        mask_conf = zeros(512, 512);

        %Add the OD mask over the correct range to change the value of the appropriate 
        %cells within that region to 1.
        mask_conf(OD_y_range, OD_x_range) = (mask_conf(OD_y_range, OD_x_range) + double(mask_OD));

        %Add the ID mask over the correct range to change the value of the appropriate 
        %cells within that region to 0. This will leave only the confining region with 
        %cells equal to 1.
        mask_conf(ID_y_range, ID_x_range) = (mask_conf(ID_y_range, ID_x_range) - double(mask_ID));
        
        %Load a slice from the WET scan
        load(strcat(mat_files, wet_scan, '/raw/slice', num2str(i), '.mat'));

        %Apply the confining mask to the scan to isolate the attenuations representative of 
        %the confining fluid
        slice = (double(slice) .* mask_conf);

        %Iterate through the slice and eliminate any value over 0 or under 200 to remove 
        %the thermistor and material at the edges of the confining fluid zone. Also, add 
        %up all the zeros in the slice to determine the number of cells that will not be 
        %included in the average.
        sum_zero = 0;
        for j = 1:512
            for k = 1:512
                if (slice(j, k) > 0 | slice(j, k) < -200)
                    slice(j, k) = 0;
                end
                
                if (slice(j, k) == 0)
                    sum_zero = sum_zero + 1;
                end
            end
        end
        
        %Add the entire confining mask together to determine the number of cells that will 
        %be included in the average
        mask_conf_sum = ((512 * 512) - sum_zero);
        
        %Determine the average attenuation of the current slice.
        avg_conf_atten(i, 1) = (sum(sum(slice)) / mask_conf_sum); 
    end
    
    %Calculate the overall average attenuation across all the slices
    tot_avg_conf_atten = (sum(avg_conf_atten) / length(avg_conf_atten));
    
    %Determine the heat correction as the difference between the average attenuation of 
    %each slice and the total average attenuation.
    heat_corr = (avg_conf_atten - tot_avg_conf_atten);
    
    %Save the heat correction to the results folder for later use.
    outfile = strcat(results, 'heat_corr.mat');
    save(outfile, 'heat_corr');
    
    %STEP 3: Apply heat correction to all slices
    
    %Create new directories to hold the corrected MAT files of the WET and DRY scans.
    wet_corr_path = strcat(mat_files, wet_scan, '/corrected/');
    dry_corr_path = strcat(mat_files, dry_scan, '/corrected/');
    mkdir(wet_corr_path);
    mkdir(dry_corr_path);
    
    %Iterate through the filtered WET and DRY scans and apply the heat correction to every 
    %slice
    for j = 1:n_slices
        %Load the current slice from the WET and DRY scans
        load(strcat(mat_files, wet_scan, '/filtered/slice', num2str(j), '.mat'));
        load(strcat(mat_files, dry_scan, '/filtered/slice', num2str(j), '.mat'));
        
        %Apply the heat correction to the entire slice and then reapply the mask to remove 
        %values outside the sample
        slice_wet_corr = ((double(slice_wet) + heat_corr(j, 1)) .* double(mask));
        slice_dry_corr = ((double(slice_dry) + heat_corr(j, 1)) .* double(mask));
        
        %Create the file name for the modified WET and DRY slices and save them to their 
        %folders
        outfile_name = strcat(wet_corr_path, 'slice', num2str(j), '.mat');
        save(outfile_name, 'slice_wet_corr');
        
        outfile_name = strcat(dry_corr_path, 'slice', num2str(j), '.mat');
        save(outfile_name, 'slice_dry_corr');

        %Clear the WET and DRY slices
        clear slice*
    end
    
    %Iterate through each of the experimental scans
    for j = 1:n_scans
        %Skip the scan if it's found in the skip list
        if (sum(j == skip_list) == 1)
            continue
        end

        %Define and create the directory that will hold the corrected slices
        path = strcat(mat_files, strjoin(scan_names(j, 1)), '/corrected/');
        mkdir(path);

        %Iterate through each slice and apply the heat correction
        for i = 1:n_slices            
            %Load the current slice
            load(strcat(mat_files, strjoin(scan_names(j, 1)), '/filtered/slice', ...
                num2str(i), '.mat'));

            %Apply the heat correction to the entire slice and then reapply the mask to 
            %remove values outside the sample
            slice_corr = ((double(slice) + heat_corr(i, 1)) .* double(mask));

            %Create the file name for the modified slice and save it
            outfile_name = strcat(path, 'slice', num2str(i), '.mat');
            save(outfile_name, 'slice_corr');

            %Clear the current slice
            clear slice*
        end
        
        %Display the completion progress of the current procedure
        fprintf('Percent Complete: %d\n', round(((j - 1) / n_scans) * 100))
    end
    
    %Display that this process is complete
    fprintf('Percent Complete: %d\n', round(100))
    
    %Update this checkpoint status
    heat_correct = 'y';
end

fprintf('DONE\n')

%---------------------------------------------------------------------------------------------------------
%CHECKPOINT 6: SAMPLE POROSITY CALCULATION
%This section takes the final attenuation data and calculates a porosity value for every
%voxel in the sample from the saturated and unsaturated scans

fprintf('--- Calculating the 3D Sample Porosity --- ')

while (porosity == 'n')
    %Create a directory to hold the porosity results
    phi_path = strcat(results, '/porosity/data/');
    mkdir(phi_path);
    
    %Initial arrays to hold calculation data
    atten_dry = zeros(1, n_slices);
    atten_wet = zeros(1, n_slices);
    phi_avg = zeros(1, n_slices);
           
    %STEP 1: Determine the denominator for the porosity calculation
    %The denominator used to calculate the porosity from the CT scans is calibrated such 
    %that the bulk average porosity of the sample is equal to the gravimetric porosity
    
    %Iterate through the slices in the saturated and unsaturated scans and determine the 
    %average CT attenuation
    for i = 1:n_slices
        %Load a slice from both the WET and DRY scans
        load(strcat(mat_files, wet_scan, '/corrected/slice', num2str(i), '.mat'));
        load(strcat(mat_files, dry_scan, '/corrected/slice', num2str(i), '.mat'));
        
        %Calculate average CT attenuation from each slice
        atten_dry(1, i) = (sum(sum(slice_dry_corr)) / mask_sum);
        atten_wet(1, i) = (sum(sum(slice_wet_corr)) / mask_sum);
        
        %Clear the WET and DRY slices
        clear slice*
    end
    
    %Calculate the overall average CT attenuations for the WET and DRY scans
    avg_atten_min = (sum(atten_dry) / n_slices);
    avg_atten_max = (sum(atten_wet) / n_slices);
    
    %Calculate the denominator that will result in a bulk average porosity that is 
    %equivalent to the gravimetric porosity
    denom_final = ((avg_atten_max - avg_atten_min) / phi);
    
    %STEP 2: Determine the 3D sample porosity
    
    %Iterate through the slices in the WET and DRY scans and calculate the porosity
    for i = 1:n_slices
        %Load a slice from both the WET and DRY scans
        load(strcat(mat_files, wet_scan, '/corrected/slice', num2str(i), '.mat'));
        load(strcat(mat_files, dry_scan, '/corrected/slice', num2str(i), '.mat'));

        %Calculate the porosity at every point in the current slice
        slice_phi = (double(slice_wet_corr - slice_dry_corr) / denom_final);
        
        %Record the average porosity in the slice
        phi_avg(1, i) = (sum(sum(slice_phi)) / mask_sum);
        
        %Create the file name for the modified WET and DRY slices and save them to their 
        %folders
        outfile_name = strcat(phi_path, 'slice', num2str(i), '.mat');
        save(outfile_name, 'slice_phi');

        %Clear the current slice
        clear slice*
    end
    
    %Calculate the overall average porosity for comparison to the gravimetric porosity
    avg_phi_final = (sum(phi_avg) / n_slices);
    
    %Update this checkpoint status
    porosity = 'y';
end

fprintf('DONE\n')

%---------------------------------------------------------------------------------------------------------
%CHECKPOINT 7: SAMPLE BULK DENSITY CALCULATION
%This section uses the porosity and CT attenuation data to calculate the bulk density and 
%change in bulk change for every voxel in the sample for every scan taken throughout the 
%experiment

fprintf('--- Calculating 3D Bulk Density --- \n')

while (density == 'n')
    
    %Define the path to the porosity and corrected WET and DRY attenuation data
    wet_corr_path = strcat(mat_files, wet_scan, '/corrected/');
    dry_corr_path = strcat(mat_files, dry_scan, '/corrected/');
    phi_path = strcat(results, 'porosity/data/');
    
    %Make directories to hold the density data the WET and DRY scans
    wet_den_path = strcat(mat_files, wet_scan, '/density/');
    dry_den_path = strcat(mat_files, dry_scan, '/density/');
    mkdir(wet_den_path);
    mkdir(dry_den_path);
    
    %Iterate through the slices in the porosity data and calculate the bulk density in the 
    %saturated and unsaturated scans
    for i = 1:n_slices        
        %Load the current porosity data
        load(strcat(phi_path, 'slice', num2str(i), '.mat'));
        
        %Calculate the density at every point in the sample
        slice_wet_den = ((((1 - slice_phi) * rhos) + (slice_phi * rhow)) .* double(mask));
        slice_dry_den = ((((1 - slice_phi) * rhos) + (slice_phi * rhog)) .* double(mask));

        %Create the file name for the WET and DRY bulk density data and save it
        outfile_name = strcat(wet_den_path, 'slice', num2str(i), '.mat');
        save(outfile_name, 'slice_wet_den');
        
        outfile_name = strcat(dry_den_path, 'slice', num2str(i), '.mat');
        save(outfile_name, 'slice_dry_den');

        %Clear the WET, DRY, and porosity data
        clear slice*
    end
    
    %Iterate through each scan
    for j = 1:n_scans
        %Skip the scan if it's found in the skip list
        if (sum(j == skip_list) == 1)
            continue
        end

        %Define and create the directory that will hold the bulk density data
        path = strcat(mat_files, strjoin(scan_names(j, 1)), '/density/');
        drho_path = strcat(mat_files, strjoin(scan_names(j, 1)), '/drho/');
        mkdir(path);
        mkdir(drho_path);

        %Iterate through each slice and calculate the bulk density data
        for i = 1:n_slices
            %Load the bulk density data for the current WET and DRY slice
            load(strcat(wet_den_path, 'slice', num2str(i), '.mat'));
            load(strcat(dry_den_path, 'slice', num2str(i), '.mat'));
            
            %Load the attenuation data for the current WET and DRY slice
            load(strcat(wet_corr_path, 'slice', num2str(i), '.mat'));
            load(strcat(dry_corr_path, 'slice', num2str(i), '.mat'));
            
            %Load the CT attenuation for the current slice
            load(strcat(mat_files, strjoin(scan_names(j, 1)), '/corrected/slice', ...
                                                                      num2str(i), '.mat'));

            %Calculate the bulk density using linear interpolation
            slice_den = (slice_wet_den + ((slice_dry_den - slice_wet_den) .* (slice_corr ...
                                         - slice_wet_corr) ./ (slice_dry_corr - ...
                                           slice_wet_corr)));
            slice_den(isnan(slice_den)) = 0;
            
            slice_drho = (slice_den - slice_wet_den);
            
            %Create the file name for the bulk density data and save it
            outfile_name = strcat(path, 'slice', num2str(i), '.mat');
            save(outfile_name, 'slice_den');
            
            outfile_name = strcat(drho_path, 'slice', num2str(i), '.mat');
            save(outfile_name, 'slice_drho');

            %Clear the WET, DRY, and experimental data
            clear slice*
        end
        
        %Display the completion progress of the current procedure
        fprintf('Percent Complete: %d\n', round(((j - 1) / n_scans) * 100))
    end
    
    %Display that this process is complete
    fprintf('Percent Complete: %d\n', round(100))
    
    %Update this checkpoint status
    density = 'y';
end

fprintf('DONE\n')

%---------------------------------------------------------------------------------------------------------
%CHECKPOINT 8: AFFECTED VOLUME CALCULATION
%This section determines the pore volume behind the front (the affected volume) using the
%change in bulk density and porosity in each voxel and the voxel volume.

fprintf('--- Calculating Affected Volume --- ')

while (v_aff == 'n')
    
    %Define the path to the porosity data
    phi_path = strcat(results, 'porosity/data/');
    
    %Define the path to where the affected volume will be stored
    Vaff_path = strcat(results, 'Vaff.mat');
            
    %Create the matricies that will hold the affected volume in mL
    V_aff = zeros(n_scans, 1);

    %Defined the pixel res and the voxel and slice volumes
    pixel_res = ((r * 2) / ((rad * 2) + 1)); %1D pixel resolution, [cm]
    voxel_vol = (slice_space * pixel_res * pixel_res); %Voxel volume, [mL]
    
    %Iterate through each scan
    for j = 2:n_scans
        %Skip the scan if it's found in the skip list
        if (sum(j == skip_list) == 1)
            continue
        end
        
        %Define the path to the change in bulk density data
        drho_path = strcat(mat_files, strjoin(scan_names(j, 1)), '/drho/');
        
        %Initialize the affected volume variable
        V_aff_sum = 0;
        
        %Iterate through each slice, calculate the affected volume, and add it to the total
        for i = 1:n_slices
            %Load the porosity and change in bulk density data for the current slice
            load(strcat(phi_path, 'slice', num2str(i), '.mat'));
            load(strcat(drho_path, 'slice', num2str(i), '.mat'));
            
            %Check which voxels have bulk a bulk density change exceeding the density error
            %and produce a matrix of ones and zeros marking the "affected" voxels.
            slice_Vaff = (slice_drho < -0.024);
            
            %Record the porosities of each affected voxel
            slice_Vaff = (slice_Vaff .* slice_phi);
            
            %Calculate the pore volume of each affected voxel
            slice_Vaff = (slice_Vaff * voxel_vol);
            
            %Sum all the pore volumes of the affected volumes in this slice and add it to
            %the running total
            V_aff_sum = (V_aff_sum + sum(sum(slice_Vaff)));
        end
        
        %Record the total affected volume for the scan
        V_aff(j, 1) = V_aff_sum;
        
        %Display the completion progress of the current procedure
        fprintf('Percent Complete: %d\n', round(((j - 1) / n_scans) * 100))
    end
    
    %Display that this process is complete
    fprintf('Percent Complete: %d\n', round(100))
       
    %Save the final front volume, location, and velocity results to the correct location
    save(Vaff_path, 'V_aff');
    
    %Update this checkpoint status
    v_aff = 'y';
end

fprintf('DONE\n')

%---------------------------------------------------------------------------------------------------------
%CHECKPOINT 9: ORTHOGONAL DATA RECORDING
%This section transforms the 2D bulk density data from axial slices to orthogonal slices by 
%recording the data along the centerline of each slice in each scan into a new matrix that
%shows the 2D distribution of the data along the centerline.

fprintf('--- Recording Orthogonal Slices --- \n')

while(ortho_data == 'n')
    
    %Define the paths to the bulk density and porosity data
    phi_path = strcat(results, 'porosity/data/');
    wet_den_path = strcat(mat_files, wet_scan, '/density/');
    dry_den_path = strcat(mat_files, dry_scan, '/density/');
    
    %Define the array that will hold the results
    ortho_phi_array = [];
    ortho_wet_array = [];
    ortho_dry_array = [];
    
    %Iterate through the slices and save the centerline data 
    for j = 1:n_slices
            %Load the current slice
            load(strcat(phi_path, 'slice', num2str(j), '.mat'));
            load(strcat(wet_den_path, 'slice', num2str(j), '.mat'));
            load(strcat(dry_den_path, 'slice', num2str(j), '.mat'));

            %Save the middle column of the array to a appropriate column in the 
            %orthogonal array
            ortho_phi_array(:, j) = slice_phi(:, rad);
            ortho_wet_array(:, j) = slice_wet_den(:, rad);
            ortho_dry_array(:, j) = slice_dry_den(:, rad);
    end
    
    %Save the orthogonal data to the appropriate directory
    outfile = strcat(results, 'ortho_porosity.mat');
    save(outfile, 'ortho_phi_array');
    
    outfile = strcat(results, 'ortho_wet_den.mat');
    save(outfile, 'ortho_wet_array');
    
    outfile = strcat(results, 'ortho_dry_den.mat');
    save(outfile, 'ortho_dry_array');
    
    %Create a directory to hold the orthogonal slice results
    den_ortho_path = strcat(results, 'ortho_data/density/');
    drho_ortho_path = strcat(results, 'ortho_data/drho/');
    mkdir(den_ortho_path);
    mkdir(drho_ortho_path);

    %Iterate through all the scans
    for i = 1:n_scans
        %Skip the scan if it's found in the skip list
        if (sum(i == skip_list) == 1)
            continue
        end

        %Create blank arrays to hold the density and saturation data
        den_ortho_array = [];
        drho_ortho_array = [];

        for j = 1:n_slices
            %Load the current slice
            load(strcat(mat_files, strjoin(scan_names(i, 1)), '/density/slice', ...
                num2str(j), '.mat'));
            load(strcat(mat_files, strjoin(scan_names(i, 1)), '/drho/slice', ...
                num2str(j), '.mat'));

            %Save the middle column of the array to a appropriate column in the orthogonal 
            %array
            den_ortho_array(:, j) = slice_den(:, rad);
            drho_ortho_array(:, j) = slice_drho(:, rad);
        end
        
        %Save the orthogonal data to the appropriate directory
        outfile = strcat(den_ortho_path, strjoin(scan_names(i, 1)), '.mat');
        save(outfile, 'den_ortho_array');
        
        outfile = strcat(drho_ortho_path, strjoin(scan_names(i, 1)), '.mat');
        save(outfile, 'drho_ortho_array');

        %Display the completion progress of the current procedure
        fprintf('Percent Complete: %d\n', round(((i - 1) / n_scans) * 100))
    end

    %Update this procedural checkpoint
    ortho_data = 'y';
end

fprintf('DONE\n')

%---------------------------------------------------------------------------------------------------------
%PLOTTING THE RESULTS
%This section creates contour plots of both the axial and orthogonal slices for the
%porosity, bulk density, and front location data for every scan

%SUBROUTINES:
%   cbarf.m

fprintf('--- Plotting Results --- \n')

%Define the number of rows and columns in the experimental data. Columns are the number of 
%slices in each sample and rows are the number of pixels across the sample.
rows = ((2* rad) + 1);
columns = n_slices;

%Load the orthogonal data
load(strcat(results, 'ortho_porosity.mat'));
load(strcat(results, 'ortho_wet_den.mat'));
load(strcat(results, 'ortho_dry_den.mat'));

%Define the minimum and maximum porosities in the sample
max_phi = max((ceil((max(max(ortho_phi_array([3:207], :))) * 100)) / 100));
min_phi = min((floor((min(min(ortho_phi_array([3:207], :))) * 100)) / 100));

%Define the minimum and maximum densities in the sample
max_tot_den = (ceil((max(max(ortho_wet_array([3:207], :))) * 100)) / 100);
min_tot_den = (floor((min(min(ortho_dry_array([3:207], :))) * 100)) / 100);

%Define the contours to be used in this plot. You may need to play with this value to get 
%a good image depending on the porosity range. Change the middle value below to change 
%the contour interval.
contours_phi = [min_phi:0.03:max_phi];
contours_den = [min_tot_den:0.1:max_tot_den];
contours_drho = [-0.024:-0.047:-0.4];

%Define the axis parameters for the orthogonal plots
X_ticks = [1 50 100 150 rows];
Y_ticks = fliplr([[columns:-5:5] 1]);
X_tick_labels = num2cell([0 1.2 2.4 3.6 5.1]);
Y_tick_labels = num2cell([L [9:-1.5:0]]);

%STEP 1: Create contour plots of all the axial slices from every scan
while (axial_plots == 'y')
    %Define the paths to the porosity and the WET and DRY density slice data
    phi_path = strcat(results, '/porosity/data/');
    wet_den_path = strcat(mat_files, wet_scan, '/density/');
    dry_den_path = strcat(mat_files, dry_scan, '/density/');

    %Make the directories to hold the porosity and the WET and DRY bulk density slice     %images
    wet_plot_path = strcat(mat_files, wet_scan, '/plots/');
    dry_plot_path = strcat(mat_files, dry_scan, '/plots/');
    phi_plot = strcat(results, 'porosity/plots/');
    mkdir(wet_plot_path);
    mkdir(dry_plot_path);
    mkdir(phi_plot);

    %Iterate through each slice and plot the porosity and the WET and DRY bulk density
    for j = 1:n_slices
        %Load the current slice
        load(strcat(wet_den_path, 'slice', num2str(j), '.mat'));

        %Open a new figure window with a white background
        figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

        %Plot the current slices on a contour plot
        contourf(slice_wet_den, contours_den);

        %Manipulate the plot properties to reduce post-processing later
        colormap(flipud(colormap));
        caxis([min(contours_den) max(contours_den)]);
        set(gca, 'box', 'off')
        axis square %Make the image square
        axis off
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', '');
        set(gca, 'YTickLabel', '');

        %Save the image to the correct place
        outfile = strcat(wet_plot_path, 'slice', num2str(j), '.jpg');
        img = getframe(gcf);
        imwrite(img.cdata, outfile)

        %Close the figure
        close all

        %Load the current slice
        load(strcat(dry_den_path, 'slice', num2str(j), '.mat'));

        %Open a new figure window with a white background
        figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

        %Plot the current slices on a contour plot
        contourf(slice_dry_den, contours_den);

        %Manipulate the plot properties to reduce post-processing later
        colormap(flipud(colormap));
        caxis([min(contours_den) max(contours_den)]);
        set(gca, 'box', 'off')
        axis square
        axis off
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', '');
        set(gca, 'YTickLabel', '');

        %Save the image to the correct place
        outfile = strcat(dry_plot_path, 'slice', num2str(j), '.jpg');
        img = getframe(gcf);
        imwrite(img.cdata, outfile)

        %Close the figure
        close all

        %Load the current dry porosity slice
        load(strcat(phi_path, 'slice', num2str(j), '.mat'));

        %Open a new figure window with a white background
        figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

        %Plot the current slices on a contour plot
        contourf(slice_phi, contours_phi);

        %Manipulate the plot properties to reduce post-processing later
        caxis([min(contours_phi) max(contours_phi)]);
        set(gca, 'box', 'off')
        axis square
        axis off
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', '');
        set(gca, 'YTickLabel', '');

        %Save the image to the correct place
        outfile = strcat(phi_plot, 'slice', num2str(j), '.jpg');
        img = getframe(gcf);
        imwrite(img.cdata, outfile)

        %Close the figure
        close all
    end

    for i = 1:n_scans
        %Skip the scan if it's found in the skip list
        if (sum(i == skip_list) == 1)
            continue
        end

        %Concatenate and create the directory that will hold the processed slices
        den_path = strcat(mat_files, strjoin(scan_names(i, 1)), '/plots/density/');
        drho_path = strcat(mat_files, strjoin(scan_names(i, 1)), '/plots/drho/');
        mkdir(den_path);
        mkdir(drho_path);

        %Iterate through the slices in the current scan and plot them on the sample figure
        for j = 1:n_slices
            
            %BULK DENSITY
            %Load the current slice
            load(strcat(mat_files, strjoin(scan_names(i, 1)), '/density/slice', num2str(j), '.mat'));

            %Open a new figure window with a white background
            figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

            %Plot the current slices on a contour plot
            contourf(slice_den, contours_den);

            %Manipulate the plot properties to reduce post-processing later
            colormap(flipud(colormap));
            caxis([min(contours_den) max(contours_den)]);
            set(gca, 'box', 'off')
            axis square
            axis off
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
            set(gca, 'XTickLabel', '');
            set(gca, 'YTickLabel', '');

            %Save the image to the correct place
            outfile = strcat(den_path, 'slice', num2str(j), '.jpg');
            img = getframe(gcf);
            imwrite(img.cdata, outfile)

            %Close the figure
            close all

            %BULK DENSITY CHANGE
            %Load the current slice
            load(strcat(mat_files, strjoin(scan_names(i, 1)), '/drho/slice', num2str(j), '.mat'));

            %Open a new figure window with a white background
            figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

            %Plot the current slices on a contour plot
            contourf(slice_drho, contours_drho);

            %Manipulate the plot properties to reduce post-processing later
            colormap(flipud(colormap));
            caxis([min(contours_drho) max(contours_drho)]);
            set(gca, 'box', 'off')
            axis square
            axis off
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
            set(gca, 'XTickLabel', '');
            set(gca, 'YTickLabel', '');

            %Save the image to the correct place
            outfile = strcat(drho_path, 'slice', num2str(j), '.jpg');
            img = getframe(gcf);
            imwrite(img.cdata, outfile)

            %Close the figure
            close all

        end

        %Display the completion progress of the current procedure
        fprintf('Percent Complete: %d\n', round(((i - 1) / n_scans) * 100))
    end

    %Display this process as complete
    fprintf('Percent Complete: %d\n', round(100))

    fprintf('DONE\n')

    %Update this checkpoint status
    axial_plots = 'n';
end

%STEP 2: Create contour plots of all the orthogonal slices from every scan
while (ortho_plots == 'y')
        
    %Open a new figure with a white background
    figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

    %Plot the orthogonal porosity data
    contourf(ortho_phi_array', contours_phi);

    %Manipulate the plot to reduce post-processing time
    xlim([1 rows])
    xlabel('Cross-Sectional Distance, y (cm)')
    ylabel('Distance from Inlet, x (cm)')
    set(gca, 'XTick', X_ticks);
    set(gca, 'YTick', Y_ticks);
    set(gca, 'XTickLabel', X_tick_labels);
    set(gca, 'YTickLabel', Y_tick_labels);
    caxis([min(contours_phi) max(contours_phi)]);
    set(gca, 'PlotBoxAspectRatio', [0.486 1 1]);
    set(gca, 'box', 'on', 'LineWidth', 2)

    %Save the plot as a JPEG into the correct folder
    outfile = strcat(results, 'ortho_porosity.jpg');
    img = getframe(gcf);
    imwrite(img.cdata, outfile)

    %Close the figure
    close all

    %Open a new figure with a white background
    figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

    %Plot the orthogonal wet bulk density data
    contourf(ortho_wet_array', contours_den);

    %Manipulate the plot to reduce post-processing time
    xlim([1 rows])
    xlabel('Cross-Sectional Distance, y (cm)')
    ylabel('Distance from Inlet, x (cm)')
    set(gca, 'XTick', X_ticks);
    set(gca, 'YTick', Y_ticks);
    set(gca, 'XTickLabel', X_tick_labels);
    set(gca, 'YTickLabel', Y_tick_labels);
    colormap(flipud(colormap));
    caxis([min(contours_den) max(contours_den)]);
    set(gca, 'PlotBoxAspectRatio', [0.486 1 1]);
    set(gca, 'box', 'on', 'LineWidth', 2)

    %Save the plot as a JPEG into the correct folder
    outfile = strcat(results, 'ortho_wet_den.jpg');
    img = getframe(gcf);
    imwrite(img.cdata, outfile)

    %Close the figure
    close all

    %Open a new figure with a white background
    figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

    %Plot the orthogonal dry bulk density data
    contourf(ortho_dry_array', contours_den);

    %Manipulate the plot to reduce post-processing time
    xlim([1 rows])
    xlabel('Cross-Sectional Distance, y (cm)')
    ylabel('Distance from Inlet, x (cm)')
    set(gca, 'XTick', X_ticks);
    set(gca, 'YTick', Y_ticks);
    set(gca, 'XTickLabel', X_tick_labels);
    set(gca, 'YTickLabel', Y_tick_labels);
    colormap(flipud(colormap));
    caxis([min(contours_den) max(contours_den)]);
    set(gca, 'PlotBoxAspectRatio', [0.486 1 1]);
    set(gca, 'box', 'on', 'LineWidth', 2)

    %Save the plot as a JPEG into the correct folder
    outfile = strcat(results, 'ortho_dry_den.jpg');
    img = getframe(gcf);
    imwrite(img.cdata, outfile)

    %Close the figure
    close all

    %Define the path to the orthogonal data
    den_ortho_path = strcat(results, 'ortho_data/density/');
    drho_ortho_path = strcat(results, ortho_data/drho');

    %Make the directory to hold the results
    den_ortho_exp_path = strcat(results, 'ortho_plots/density/');
    drho_ortho_exp_path = strcat(results, 'ortho_plots/drho/');
    mkdir(den_ortho_exp_path);
    mkdir(drho_ortho_exp_path);

    %Create contour plots of the experimental scans
    for i = 1:n_scans
        %Skip the scan if it's found in the skip list
        if (sum(i == skip_list) == 1)
            continue
        end

        %BULK DENSITY
        %Load the current orthogonal data
        load(strcat(den_ortho_path, strjoin(scan_names(i, 1)), '.mat'));

        %Open a new figure window with a white background
        figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

        %Plot the current ortho data
        contourf(den_ortho_array', contours_den);

        %Manipulate the plot to reduce post-processing time
        xlim([1 rows])
        title_str = strcat('t = ', num2str(round(scan_times(i))), 'hrs');
        title(title_str);
        xlabel('Cross-Sectional Distance, y (cm)')
        ylabel('Distance from Inlet, x (cm)')
        set(gca, 'XTick', X_ticks);
        set(gca, 'YTick', Y_ticks);
        set(gca, 'XTickLabel', X_tick_labels);
        set(gca, 'YTickLabel', Y_tick_labels);
        set(gca, 'PlotBoxAspectRatio', [0.486 1 1]);
        colormap(flipud(colormap));
        caxis([min(contours_den) max(contours_den)]);

        %Save the image out to the correct location
        outfile = strcat(den_ortho_exp_path, strjoin(scan_names(i, 1)), '.jpg');
        img = getframe(gcf);
        imwrite(img.cdata, outfile)

        %Close the figure
        close all

        %BULK DENSITY CHANGE
        %Load the current orthogonal data
        load(strcat(drho_ortho_path, strjoin(scan_names(i, 1)), '.mat'));

        %Open a new figure window with a white background
        figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

        %Plot the current ortho data
        contourf(drho_ortho_array, contours_drho);

        %Manipulate the plot to reduce post-processing time
        xlim([1 rows])
        title_str = strcat('t = ', num2str(round(scan_times(i))), 'hrs');
        title(title_str);
        xlabel('Cross-Sectional Distance, y (cm)')
        ylabel('Distance from Inlet, x (cm)')
        set(gca, 'XTick', X_ticks);
        set(gca, 'YTick', Y_ticks);
        set(gca, 'XTickLabel', X_tick_labels);
        set(gca, 'YTickLabel', Y_tick_labels);
        set(gca, 'PlotBoxAspectRatio', [0.486 1 1]);
        colormap(flipud(colormap));
        caxis([min(contours_drho) max(contours_drho)]);

        %Save the image out to the correct location
        outfile = strcat(drho_ortho_exp_path, strjoin(scan_names(i, 1)), '.jpg');
        img = getframe(gcf);
        imwrite(img.cdata, outfile)

        %Close the figure
        close all
    end

    fprintf('DONE\n')

    %Update this procedural checkpoint
    ortho_plots = 'n';
end

%STEP 3: Create contoured colorbars for each data type plotted above
%This process is done separately from the plotting in order to preserve the correct aspect %ratio in the axial and orthogonal plots
while (colorbars == 'y')
   
    %Open a new figure with a white background
    figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

    %Plot the orthogonal porosity data
    contourf(ortho_phi_array', contours_phi);

    %Manipulate the plot to reduce post-processing time
    xlim([1 rows])
    xlabel('Cross-Sectional Distance, y (cm)')
    ylabel('Distance from Inlet, x (cm)')
    set(gca, 'XTick', X_ticks);
    set(gca, 'YTick', Y_ticks);
    set(gca, 'XTickLabel', X_tick_labels);
    set(gca, 'YTickLabel', Y_tick_labels);
    set(gca, 'PlotBoxAspectRatio', [0.486 1 1]);
    set(gca, 'box', 'on', 'LineWidth', 2)
    caxis([min(contours_phi) max(contours_phi)])
    cbarf([(min(contours_phi) - 0.1) contours_phi (max(contours_phi) + 0.1)], contours_phi);

    %Save the plot as a JPEG into the correct folder
    outfile = strcat(results, 'porosity_colorbar.jpg');
    saveas(gcf, outfile, 'jpeg');

    %Close the figure
    close all

    %Open a new figure with a white background
    figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

    %Plot the orthogonal dry bulk density data
    contourf(ortho_dry_array', contours_den);

    %Manipulate the plot to reduce post-processing time
    xlim([1 rows])
    xlabel('Cross-Sectional Distance, y (cm)')
    ylabel('Distance from Inlet, x (cm)')
    set(gca, 'XTick', X_ticks);
    set(gca, 'YTick', Y_ticks);
    set(gca, 'XTickLabel', X_tick_labels);
    set(gca, 'YTickLabel', Y_tick_labels);
    set(gca, 'PlotBoxAspectRatio', [0.486 1 1]);
    set(gca, 'box', 'on', 'LineWidth', 2)
    colormap(flipud(colormap));
    caxis([min(contours_den) max(contours_den)])
    cbarf([(min(contours_den) - 0.1) contours_den (max(contours_den) + 0.1)], contours_den);

    %Save the plot as a JPEG into the correct folder
    outfile = strcat(results, 'density_colorbar.jpg');
    saveas(gcf, outfile, 'jpeg');

    %Close the figure
    close all

    %Open a new figure with a white background
    figure('Color', 'W', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])

    %Plot the orthogonal dry bulk density data
    contourf(drho_ortho_array', contours_drho);

    %Manipulate the plot to reduce post-processing time
    xlim([1 rows])
    xlabel('Cross-Sectional Distance, y (cm)')
    ylabel('Distance from Inlet, x (cm)')
    set(gca, 'XTick', X_ticks);
    set(gca, 'YTick', Y_ticks);
    set(gca, 'XTickLabel', X_tick_labels);
    set(gca, 'YTickLabel', Y_tick_labels);
    set(gca, 'PlotBoxAspectRatio', [0.486 1 1]);
    set(gca, 'box', 'on', 'LineWidth', 2)
    colormap(flipud(colormap));
    caxis([min(contours_drho) max(contours_drho)])
    cbarf([(min(contours_drho) - 0.1) contours_drho (max(contours_drho) + 0.1), contours_drho]);

    %Save the plot as a JPEG into the correct folder
    outfile = strcat(results, 'drho_colorbar.jpg');
    saveas(gcf, outfile, 'jpeg');

    %Close the figure
    close all
    
    %Update this checkpoint status
    colorbars = 'n';
end
