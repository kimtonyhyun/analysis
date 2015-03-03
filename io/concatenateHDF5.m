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
% To Do: 
% Write a more efficient method to replace dropped frames
%
% 2015-02-03 Jessica Maxey

[frame_indices,location_info,time] = parse_plusmaze(fullfile(tifDir,plusmazeName));
startFrames = frame_indices(:,1);

list = dir(fullfile(tifDir,'*.tif'));
firstName = fullfile(tifDir,list(1).name);
[rows,cols] = size(imread(firstName));

totalFrames = 1;

locLUT = {'north';'south';'east';'west'};

downsmpRows = floor(rows*downsmpFactor);
downsmpCols = floor(cols*downsmpFactor);
chunkSize = [downsmpRows downsmpCols 1];
h5create(fullfile(outputDir,hdf5Name),'/Data/Images',[Inf Inf Inf],'ChunkSize',chunkSize,'Datatype','uint16');
h5create(fullfile(outputDir,hdf5Name),'/TrialInfo/Frames',[Inf,4],'ChunkSize',[1,4]);
h5create(fullfile(outputDir,hdf5Name),'/TrialInfo/Locations',[Inf,3],'ChunkSize',[1,3],'Datatype','uint8');
h5create(fullfile(outputDir,hdf5Name),'/TrialInfo/Time',[Inf,1],'ChunkSize',[1,1]);
badTrial = 0;
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

    if(~badTrial)
        %%% Add the good trial to the hdf5 file
        trialCount = trialCount+1;
        tifName = fullfile(tifDir,list(i).name);
        tifInfo = imfinfo(tifName);
        tifFile = Tiff(tifName,'r');
        oriFrames = length(tifInfo);
        
        [path,name,ext] = fileparts(tifName);
        xmlName = [fullfile(tifDir,name),'.xml'];
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
        imageStack = zeros(downsmpRows,downsmpCols,oriFrames-trim(1)-trim(2),'uint16');
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
        cropDroppedFrames = 0;
        if(droppedFrames ~= 0)
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
                else
                    cropDroppedFrames = cropDroppedFrames+1;
                end
            end
        end
        
        numFrames = size(imageStack,3);
        [dRows,dCols] = size(imageStack(:,:,1));
        
        h5write(fullfile(outputDir,hdf5Name),'/Data/Images',imageStack,[1,1,totalFrames],[dRows,dCols,numFrames]);
        
        %%% Find frames that correspond to the gate going up and down for
        %%% the current trial
        gateUpFrame = totalFrames + (frame_indices(trialCount,2) - frame_indices(trialCount,1) - trim(1));
        gateDownFrame = totalFrames + (frame_indices(trialCount,3) - frame_indices(trialCount,1) - trim(1));
        
        frameData = [totalFrames,gateUpFrame,gateDownFrame,totalFrames+numFrames-1];
        
        %%% Convert location information to numeric representation ('north' = 1, 'south' = 2, 'east' = 3, 'west' = 4)
        startLoc = uint8(find(ismember(locLUT,location_info{trialCount,1})));
        correctEnd = uint8(find(ismember(locLUT,location_info{trialCount,2})));
        actualEnd = uint8(find(ismember(locLUT,location_info{trialCount,3})));
        locData = [startLoc, correctEnd, actualEnd];
        
        %%% Write the trial info to the hdf5 file
        h5write(fullfile(outputDir,hdf5Name),'/TrialInfo/Frames',frameData,[trialCount,1],[1,4]);
        h5write(fullfile(outputDir,hdf5Name),'/TrialInfo/Locations',locData,[trialCount,1],[1,3]);
        h5write(fullfile(outputDir,hdf5Name),'/TrialInfo/Time',time(trialCount),[trialCount,1],[1,1]);
        
        fprintf('File %d of %d stored\n',i,size(list,1));
        
        %%% Total frame count stored in hdf5 file
        totalFrames = totalFrames+numFrames;
        
        %%% Total frame count corresponding to behavior text file
        testFrames = testFrames+oriFrames+cropDroppedFrames;
        
        clear imageStack tifName tifInfo tifFile droppedFrames newImageStack
    end
end

[path,name,ext] = fileparts(list(1).name);
xmlName = [fullfile(tifDir,name),'.xml'];
xmlData = parse_miniscope_xml(xmlName);
frameRate = str2num(xmlData.fps);

%%% Remove the extra frame count needed to index the hdf5 file
totalFrames = totalFrames-1;

h5create(fullfile(outputDir,hdf5Name),'/Params/NumFrames',1);
h5write(fullfile(outputDir,hdf5Name),'/Params/NumFrames',totalFrames);
h5create(fullfile(outputDir,hdf5Name),'/Params/NumRows',1,'Datatype','uint16');
h5write(fullfile(outputDir,hdf5Name),'/Params/NumRows',uint16(dRows));
h5create(fullfile(outputDir,hdf5Name),'/Params/NumCols',1,'Datatype','uint16');
h5write(fullfile(outputDir,hdf5Name),'/Params/NumCols',uint16(dCols));
h5create(fullfile(outputDir,hdf5Name),'/Params/TrimVals',[1 2],'Datatype','uint16');
h5write(fullfile(outputDir,hdf5Name),'/Params/TrimVals',uint16(trim));
h5create(fullfile(outputDir,hdf5Name),'/Params/DownsmpFactor',1,'Datatype','double');
h5write(fullfile(outputDir,hdf5Name),'/Params/DownsmpFactor',downsmpFactor);
h5create(fullfile(outputDir,hdf5Name),'/Params/FrameRate',1,'Datatype','double');
h5write(fullfile(outputDir,hdf5Name),'/Params/FrameRate',frameRate);

h5disp(fullfile(outputDir,hdf5Name));