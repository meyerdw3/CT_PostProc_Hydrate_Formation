%This function converts the PILE files, output by the computed tomography
%scanner in the Petroleum and Geosystem Engineering Department at the
%University of Texas at Austin, from binary format to a set of 2D matrices
%of raw attenuation values.

%Inputs:
%   exp_name = experiment/sample name
%   file_name = location of the PILE file
%   scan_name = complete name of the PILE file
%   n_slices = number of CT slices per scan

%Outputs:
%   slice.mat = 2D matrix of attenuation data for a particular slice. 
%               Separate files will be created for each slice in the scan.
%   slice.tif = Grayscale image of attenuation data for a particular slice. 
%               Separate files will be created for each slice in the scan.

function [rawdata_3D] = Pile_to_tif_mat(exp_name, file_name, scan_name, n_slices)

%Load file information
s = dir(file_name);

%Extract file size
filesize = s.bytes;

%Open the file
fid = fopen (file_name);

%Read in data from the file in the correct format
rawdata_integer = fread (fid,512*512*n_slices,'int16','ieee-le');

%Recast the data as a double
rawdata_double = cast(rawdata_integer,'int16');
rawdata_3D2 =rawdata_double;

%Reshape the data into a 3D matrix
rawdata_3D = reshape(rawdata_3D2,512,512,n_slices);

%Create directories to hold the output files
tif_dir = strcat('./', exp_name, '/tif_files/', scan_name);
mat_dir = strcat('./', exp_name, '/mat_files/', scan_name, '/raw');
mkdir(tif_dir);
mkdir(mat_dir);

%Iterate through the reshaped matrix and save the .mat and .tif files
temp = 'slice';
for i=1:n_slices;
    %Show athe .tif file in grayscale and save it
    imshow(rawdata_3D(:,:,i),[500 1500]);
    colormap('gray');
    name = strcat(tif_dir,'/',temp,int2str(i),'.tif');
    X=getframe(gca);
    imwrite(X.cdata, name);

    %Save individual .mat files
    slice = rawdata_3D(:,:,i);
    name = strcat(mat_dir,'/',temp,int2str(i),'.mat');  
    save(name,'slice');
end
 