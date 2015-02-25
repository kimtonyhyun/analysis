function norm_by_disk_filter(inputDir,inputName,inputDataset,outputDir,outputName,outputDataset,varargin)
%Normalize every frame in the movie by a disk filtered version of itself
%
%   inputDir, inputName, inputDataset = Path, name, and dataset name of the
%   input hdf5 movie
%   outputDir, outputName, outputDataset = Path, name, and dataset name of
%   the output hdf5 movie
%   disk_radius may be passed as an input argument. Default for disk_radius
%   is 15.
%
% 2015 01 31 Tony Hyun Kim (Latest Revision: Jessica Maxey, 15-Feb-23)

if isempty(varargin)
    disk_radius = 15;
elseif length(varargin) == 1
    disk_radius = varargin{1};
else
    error('Only 1 variable input argument is allowed');
end

h5att = h5info(fullfile(inputDir,inputName),inputDataset);
imageStackSize = h5att.Dataspace.Size;
rows = imageStackSize(1);
cols = imageStackSize(2);
num_frames = imageStackSize(3);

h5create(fullfile(outputDir,outputName),outputDataset,[rows,cols,num_frames],'ChunkSize',[rows,cols,1],'Datatype','single');

% Apply spatial normalization
hDisk = fspecial('disk', disk_radius);

for i=1:num_frames
    if (mod(i,1000)==0)
        fprintf('  Frames %d of %d done\n', i, num_frames);
    end
    m = single(h5read(fullfile(inputDir,inputName),inputDataset,[1,1,i],[rows,cols,1]));
    m_f = imfilter(m, hDisk, 'replicate');
    M_norm = m./m_f;
    h5write(fullfile(outputDir,outputName),outputDataset,M_norm,[1,1,i],[rows,cols,1]);
end
