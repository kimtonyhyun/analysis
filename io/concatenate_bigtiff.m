function concatenate_bigtiff(plusmaze_source, path_to_xml, trim)
% Concatenate the miniscope recordings into a single BigTiff file.
%   Compatible with the plus maze output as of Cohort 9 (where we
%   explicitly turn off the scope during the inter-trial intervals).
%
% Inputs:
%   plusmaze_source: Text file output from the plus maze
%   path_to_xml: Path to folder that contains Miniscope XML/TIF pairs
%   trim: Number of frames to drop from the beginning and end of each 
%       trial (i.e. to account for dark frames at beginning)
%
% Example usage:
%   concatenate_bigtiff('mouse7_day21_ego-left.txt', pwd, [10 5]);
%
% Tony Hyun Kim (2015 Jan 17)

% Get the trial frames (beginning and end) according to plusmaze output
% and trim the beginning and end
frame_indices = get_trial_frame_indices(plusmaze_source);
trial_frame_indices = [frame_indices(:,1)+trim(1) frame_indices(:,4)-trim(2)];

% Read all XML files at the specified path
xml_files = dir(fullfile(path_to_xml,'recording_*.xml'));
num_files = length(xml_files);

% Open one TIF file, so that we can determine width and height
xml_file = fullfile(path_to_xml, xml_files(1).name);
tif_file = convert_extension(xml_file, 'tif');
info = imfinfo(tif_file);

% TIFF parameters taken from Mosaic output
tagstruct.ImageLength = info(1).Height;
tagstruct.ImageWidth = info(1).Width;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip = 5;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Orientation = Tiff.Orientation.TopLeft;

% Setup output in BigTiff format
bigtiff_output = convert_extension(plusmaze_source, 'tif');
t = Tiff(bigtiff_output, 'w8');
t.setTag(tagstruct);

% Begin write
frame_idx = 1;
write_idx = 1;
for file_idx = 1:num_files
    % TODO: Handle dropped frames as indicated in XML
    xml_file = fullfile(path_to_xml, xml_files(file_idx).name);
    
    tif_file = convert_extension(xml_file, 'tif');
    info = imfinfo(tif_file);
    num_frames_in_file = length(info);
    fprintf('%s: File %s (%d of %d)\n', datestr(now),...
        xml_file, file_idx, num_files);
    
    for k = 1:num_frames_in_file
        % FIXME: checking individual frames whether it belongs in the
        %   trial frames via 'is_trial_frame' is inefficient...
        if (is_trial_frame(trial_frame_indices, frame_idx))
            frame = imread(tif_file, k);
            
            % Save into BigTiff
            if (write_idx ~= 1)
                t.writeDirectory();
                t.setTag(tagstruct);
            end
            t.write(frame);
            write_idx = write_idx + 1;
        end
        frame_idx = frame_idx + 1;
    end
end
t.close();