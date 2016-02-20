function concatenateHDF5(tifDir,outputDir,hdf5Name,plusmazeName,trim)
% Concatenates all tif files in the specified directory into a single
% hdf5 file. The number of frames to be dropped from the beginning 
% and end of each trial is set in 'trim'. Bad trials are removed (using the
% plusmaze text file) and if frames were dropped during a good trial, the 
% previous frame is used to fill the gap.
%
% Example arguments:
% tifDir = '/Volumes/COHORT9/cohort9-herrin224/mouse5/day20_ego-right';
% outputDir = '/Users/jmaxey/Documents/MATLAB/PreFrontal';
% hdf5Name = 'test.hdf5';
% plusmazeName = 'mouse5_day20_ego-right.txt';
% trim = [10,5];
%
% To Do: 
% Write a more efficient method to replace dropped frames
%

[frame_indices, ~, ~] = parse_plusmaze(fullfile(tifDir,plusmazeName));

list = dir(fullfile(tifDir,'*.tif'));
num_files = length(list);

% Determine the image size and recording FPS
%------------------------------------------------------------
firstName = fullfile(tifDir,list(1).name);
[rows,cols] = size(imread(firstName));

xmlName = convert_extension(firstName, 'xml');
xmlData = parse_miniscope_xml(xmlName);
frameRate = str2double(xmlData.fps);

% Initialize output HDF5 datset
movie_dataset = '/Data/Images';
chunkSize = [rows cols 1];
h5create(fullfile(outputDir,hdf5Name),movie_dataset,[Inf Inf Inf],'ChunkSize',chunkSize,'Datatype','uint16');

% Indexing variables
trialCount = 0;
totalFrames = 0;

startFrames = frame_indices(:,1);
testFrames = startFrames(1,1);

for i=1:num_files
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
    else
        if(testFrames ~= 1)
            %%% Bad trial
            badTrial = 1;
            tif_info = imfinfo(fullfile(tifDir,list(i).name));
            numBadFrames = length(tif_info);
            testFrames = 1 + numBadFrames;
        else
            %%% Good trial
            badTrial = 0;
        end
    end

    if (badTrial)
        fprintf('  %d: File "%s" skipped\n', i, list(i).name);
    else
        % Add the good trial to the hdf5 file
        trialCount = trialCount+1;
        
        % Read in the images
        tifName = fullfile(tifDir,list(i).name);
        imageStack = load_movie_from_tif(tifName);
        oriFrames = size(imageStack, 3);
        
        xmlName = convert_extension(tifName, 'xml');
        xmlData = parse_miniscope_xml(xmlName);
        
        %%% Determine if frames were dropped by the miniscope during the
        %%% trial
        if(str2num(xmlData.dropped_count) ~= 0)
            droppedFrames = str2num(xmlData.dropped);
        else
            droppedFrames = 0;
        end
        
        %%% Replace any dropped frames with the frame immediately
        %%% preceeding it
        numDroppedFrames = 0;
        if(droppedFrames ~= 0)
            for j=1:length(droppedFrames)
                fprintf('DroppedFrame: %i\n',droppedFrames(j));
                droppedFrame = droppedFrames(j);
                frontStack = imageStack(:,:,1:droppedFrame-1);
                frontStack = cat(3,frontStack,imageStack(:,:,droppedFrame-1));
                backStack = imageStack(:,:,droppedFrame:end);
                newImageStack = cat(3,frontStack,backStack);
                clear imageStack
                imageStack = newImageStack;
                numDroppedFrames = numDroppedFrames+1;
            end
        end
        
        % Trim frames from the beginning and end
        finalImageStack = imageStack(:,:,1+trim(1):end-trim(2));
        numFrames = size(finalImageStack,3);
        [dRows,dCols] = size(finalImageStack(:,:,1));
        
        h5write(fullfile(outputDir,hdf5Name),movie_dataset,finalImageStack,[1,1,1+totalFrames],[dRows,dCols,numFrames]);      
        fprintf('  %d: File "%s" stored\n', i, list(i).name);
        
        %%% Total frame count stored in hdf5 file
        totalFrames = totalFrames+numFrames;
        
        %%% Total frame count corresponding to behavior text file
        testFrames = testFrames+oriFrames+numDroppedFrames;
        
        clear imageStack tifName tifInfo tifFile droppedFrames newImageStack
    end
end

h5create(fullfile(outputDir,hdf5Name),'/Params/TrimVals',[1 2],'Datatype','double');
h5write(fullfile(outputDir,hdf5Name),'/Params/TrimVals',trim);
h5create(fullfile(outputDir,hdf5Name),'/Params/FrameRate',1,'Datatype','double');
h5write(fullfile(outputDir,hdf5Name),'/Params/FrameRate',frameRate);

h5disp(fullfile(outputDir,hdf5Name));