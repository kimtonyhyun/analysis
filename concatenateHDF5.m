function [hdf5Name] = concatenateHDF5(tifDir)
% Concatenates all tif files in the specified directory into an
% hdf5 file. A file structure similar to the one below should be used. The
% stored files are all space binned by the factor set in 'downsmpFactor'.
%
% Example file structure: "G:\cohort9-herrin224\mouse5\day10_allo-south"
% Example use: "concatenateHDF5(tifDir)" where 'tifDir' points to the
%   folder containing the tif files
%
% 2015-01-06 Jessica Maxey

cd(tifDir);
list = ls('*.tif');
firstName = strtrim(list(1,:));
[rows,cols] = size(imread(firstName));
numFiles = size(list,1);
total_frames = count_frames_in_tif(tifDir);

totalFrames = 1;
[drive,buffer] = strtok(tifDir,'\');
[cohort,buffer] = strtok(buffer,'\');
[mouse,buffer] = strtok(buffer,'\');
[day,buffer] = strtok(buffer,'\');

hdf5Name = [cohort,'_',mouse,'_',day,'.h5'];
dataset = ['/',mouse];

downsmpFactor = 0.5;
chunkSize = [rows*downsmpFactor cols*downsmpFactor 1];
h5create(hdf5Name,dataset,[Inf Inf Inf],'ChunkSize',chunkSize);
for i=1:size(list,1)
    tifName = strtrim(list(i,:));
    tifInfo = imfinfo(tifName);
    tifFile = Tiff(tifName,'r');
    numFrames = length(tifInfo);
    for j=1:numFrames
        tifFile.setDirectory(j);
        imageStack(:,:,j) = imresize(tifFile.read(),downsmpFactor,'bilinear');
    end
    [dRows,dCols] = size(imageStack(:,:,1));

    h5write(hdf5Name,dataset,imageStack,[1,1,totalFrames],[dRows,dCols,numFrames]);
   
    fprintf('File %d of %d stored\n',i,size(list,1));
    totalFrames = totalFrames+numFrames;

    clear imageStack tifName tifInfo numFrames
end

h5disp(hdf5Name);