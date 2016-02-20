function concatenateHDF5(tifDir, outputDir, hdf5Name, plusmazeName, trim, varargin)
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

ignore_xml = 0;
for k = 1:length(varargin)
    if ischar(varargin{k})
        vararg = lower(varargin{k});
        switch vararg
            case {'noxml', 'ignorexml'}
                fprintf('Ignoring XMLs!\n');
                ignore_xml = 1;
        end
    end
end

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
h5create(fullfile(outputDir,hdf5Name), movie_dataset, [Inf Inf Inf],...
         'ChunkSize', [rows cols 1],...
         'Datatype', 'uint16');

% Indexing variables
trialCount = 0;  % Number of trials saved to HDF5
totalFrames = 0; % Number of frames saved to HDF5

startFrames = frame_indices(:,1);
testFrame = 1;

for i=1:num_files
    % Load recording
    tifName = fullfile(tifDir,list(i).name);
    imageStack = load_movie_from_tif(tifName);
    
    % Optionally, check XML for dropped frames
    if ~ignore_xml
        xmlName = convert_extension(tifName, 'xml');
        xmlData = parse_miniscope_xml(xmlName);

        numDroppedFrames = str2double(xmlData.dropped_count);
        if(numDroppedFrames ~= 0)
            droppedFrames = str2num(xmlData.dropped);
        else
            droppedFrames = [];
        end

        % Fill in dropped frames by its previous frame. Remarks:
        %   1) Code will break if the first frame of file has been dropped.
        %   2) Assumes that 'droppedFrames' is ascending.
        for j=1:length(droppedFrames)
            droppedFrame = droppedFrames(j);
            fprintf('DroppedFrame: %i\n',droppedFrame);

            frontStack = imageStack(:,:,1:droppedFrame-1);
            backStack = imageStack(:,:,droppedFrame:end);
            prev_frame = imageStack(:,:,droppedFrame-1);

            imageStack = cat(3, frontStack, prev_frame, backStack);
        end
        clear frontStack backStack prev_frame;
    end
    
    % Save frames to HDF5, if part of a good trial
    %------------------------------------------------------------
    if (testFrame == startFrames(1+trialCount)) % GOOD TRIAL        
        trialCount = trialCount+1;

        % Trim frames from the beginning and end
        frames_to_save = imageStack(:,:,1+trim(1):end-trim(2));
        numFrames = size(frames_to_save,3);
        
        h5write(fullfile(outputDir,hdf5Name), movie_dataset,...
                frames_to_save,...
                [1,1,1+totalFrames],...
                size(frames_to_save));      
        fprintf('%d: File "%s" stored\n', i, list(i).name);
        
        % Total frame count stored in hdf5 file
        totalFrames = totalFrames+numFrames;
    else % BAD TRIAL
        fprintf('%d: File "%s" skipped\n', i, list(i).name);
    end

    % Increment
    testFrame = testFrame + size(imageStack,3);
end

h5create(fullfile(outputDir,hdf5Name),'/Params/TrimVals',[1 2],'Datatype','double');
h5write(fullfile(outputDir,hdf5Name),'/Params/TrimVals',trim);
h5create(fullfile(outputDir,hdf5Name),'/Params/FrameRate',1,'Datatype','double');
h5write(fullfile(outputDir,hdf5Name),'/Params/FrameRate',frameRate);

h5disp(fullfile(outputDir,hdf5Name));

end % concatenateHDF5