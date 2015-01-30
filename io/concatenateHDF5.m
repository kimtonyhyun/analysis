function concatenateHDF5(tifDir,outputDir,hdf5Name,downsmpFactor,plusmazeName,trim)

% Concatenates all tif files in the specified directory into an
% hdf5 file. The stored files are all downsampled by the factor set in 
% 'downsmpFactor'. The number of frames to be dropped from the beginning 
% and end of each trial is set in 'trim'. Bad trials are removed (using the
% plusmaze text file) and if frames were dropped during a good trial, the 
% previous frame is used to fill the gap. A new text file is produced
% indicating the start and stop frames for each trial (as well as the gate
% up/down frames, the direcitons of the trial and time to complete each
% trial).
%
% Example arguments:
% tifDir = '/Volumes/COHORT9/cohort9-herrin224/mouse5/day20_ego-right';
% outputDir = '/Users/jmaxey/Documents/MATLAB/PreFrontal';
% hdf5Name = 'test.hdf5';
% downsmpFactor = 0.5;
% plusmazeName = 'mouse5_day20_ego-right.txt';
% trim = [10,5];
%
% 2015-01-27 Jessica Maxey

outputFileName = fullfile(outputDir,[plusmazeName(1:end-4),'_new.txt']);
outputFile = fopen(outputFileName,'wt+');

[frame_indices,location_info,time] = getTrialIndicies(fullfile(tifDir,plusmazeName));
startFrames = frame_indices(:,1);

list = dir(fullfile(tifDir,'*.tif'));
firstName = fullfile(tifDir,list(1).name);
[rows,cols] = size(imread(firstName));
numFiles = length(list);

totalFrames = 1;

downsmpRows = floor(rows*downsmpFactor);
downsmpCols = floor(cols*downsmpFactor);
chunkSize = [downsmpRows downsmpCols 1];
dataset = '/Data/Images';
h5create(fullfile(outputDir,hdf5Name),dataset,[Inf Inf Inf],'ChunkSize',chunkSize,'Datatype','uint16');
badTrial = 0;
skipCount = 0;
testFrames = 1;
trialCount = 0;

for i=1:length(list)
    if(i ~= 1)
        %%% Find if the start frame of the next recording matches one of
        %%% the start frame indicies
        idx = find(testFrames == startFrames);
        if(isempty(idx))
            %%% Bad trial
            badTrial = 1;
            tif_info = imfinfo(fullfile(tifDir,list(i).name));
            numBadFrames = length(tif_info);
            testFrames = testFrames + numBadFrames;
        else
            %%% Good trial
            badTrial = 0;
        end
    end

    if(badTrial == 0)
        %%% Add the good trial to the hdf5 file
        trialCount = trialCount+1;
        tifName = fullfile(tifDir,list(i).name);
        tifInfo = imfinfo(tifName);
        tifFile = Tiff(tifName,'r');
        oriFrames = length(tifInfo);
        
        xmlName = [tifName(1:end-3),'xml'];
        xmlData = parse_miniscope_xml(xmlName);
        
        %%% Determine if frames were dropped by the miniscope during the
        %%% trial
        if(str2num(xmlData.dropped_count) ~= 0)
            droppedFrames = str2num(xmlData.dropped);
        else
            numFrames = oriFrames;
            droppedFrames = 0;
        end
        
        %%% Read in the images, trim the frames at the beginning and end 
        %%% of the trial (set by the 'trim' argument), and downsample by the 
        %%% 'downsmpFactor' argument
        imageStack = uint16(zeros(floor(rows*downsmpFactor),floor(cols*downsmpFactor),oriFrames-trim(1)-trim(2)));
        if(downsmpFactor == 1)
            for j=1:oriFrames-trim(1)-trim(2)
                tifFile.setDirectory(j+trim(1));
                imageStack(:,:,j) = uint16(tifFile.read());
            end
        else
            for j=1:oriFrames-trim(1)-trim(2)
                tifFile.setDirectory(j+trim(1));
                imageStack(:,:,j) = uint16(imresize(tifFile.read(),downsmpFactor,'bilinear'));
            end
        end
        
        %%% Replace any dropped frames with the frame immediately
        %%% preceeding it
        if(droppedFrames ~=0)
            for j=1:length(droppedFrames)
                disp('DroppedFrame');
                droppedFrame = droppedFrames(j);
                if((droppedFrame > trim(1)) && (droppedFrame < oriFrames-trim(2)))
                    droppedFrame = droppedFrame+trim(1);
                    frontStack = imageStack(:,:,1:droppedFrame-1);
                    frontStack = cat(3,frontStack,imageStack(:,:,droppedFrame-1));
                    backStack = imageStack(:,:,droppedFrame:end);
                    newImageStack = cat(3,frontStack,backStack);
                    clear imageStack
                    imageStack = newImageStack;
                end
            end
        end
        
        numFrames = size(imageStack,3);
        [dRows,dCols] = size(imageStack(:,:,1));
        
        h5write(fullfile(outputDir,hdf5Name),dataset,imageStack,[1,1,totalFrames],[dRows,dCols,numFrames]);
        
        %%%%%%%%%%%%%% NEED TO FIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Currently does not take into account added frames if some were
        %%% originally dropped - only tested on cases without dropped
        %%% frames
        gateUpFrame = totalFrames + (frame_indices(trialCount,2) - frame_indices(trialCount,1) - trim(1));
        gateDownFrame = totalFrames + (frame_indices(trialCount,3) - frame_indices(trialCount,1) - trim(1));
        fprintf(outputFile,'%s %s %s %f %d %d %d %d\n',location_info{trialCount,1},...
            location_info{trialCount,2},location_info{trialCount,3},time(trialCount),totalFrames,gateUpFrame,gateDownFrame,totalFrames+numFrames-1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fprintf('File %d of %d stored\n',i,size(list,1));
        %%% Total frame count stored in hdf5 file
        totalFrames = totalFrames+numFrames;
        
        %%% Total frame count corresponding to behavior text file
        testFrames = testFrames+oriFrames;
        
        clear imageStack tifName tifInfo tifFile droppedFrames newImageStack
    end
end

%%% Remove the extra frame count needed to index the hdf5 file
totalFrames = totalFrames-1;

h5create(fullfile(outputDir,hdf5Name),'/Parameters/Number of Frames',1);
h5write(fullfile(outputDir,hdf5Name),'/Parameters/Number of Frames',totalFrames);
h5create(fullfile(outputDir,hdf5Name),'/Parameters/Number of Rows',1,'Datatype','uint16');
h5write(fullfile(outputDir,hdf5Name),'/Parameters/Number of Rows',dRows);
h5create(fullfile(outputDir,hdf5Name),'/Parameters/Number of Columns',1,'Datatype','uint16');
h5write(fullfile(outputDir,hdf5Name),'/Parameters/Number of Columns',dCols);
h5create(fullfile(outputDir,hdf5Name),'/Parameters/Trim Values',[1 2],'Datatype','uint16');
h5write(fullfile(outputDir,hdf5Name),'/Parameters/Trim Values',trim);
h5create(fullfile(outputDir,hdf5Name),'/Parameters/Downsample Factor',1,'Datatype','uint16');
h5write(fullfile(outputDir,hdf5Name),'/Parameters/Downsample Factor',downsmpFactor);

h5disp(fullfile(outputDir,hdf5Name));