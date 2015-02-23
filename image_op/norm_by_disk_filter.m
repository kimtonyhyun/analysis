function norm_by_disk_filter(hdf5Dir,hdf5Name,dataset,varargin)
%Normalize every frame in the movie by a disk filtered version of itself
%
%   movie: movie matrix, [h x w x num_frames]
%   disk_radius may be passed as an input argument. Default for disk_radius
%   is 15.
% 2015 01 31 Tony Hyun Kim (Latest Revision: Jessica Maxey, 15-Feb-23)
%

if isempty(varargin)
    disk_radius = 15;
elseif length(varargin) == 1
    disk_radius = varargin{1};
else
    error('Only 1 variable input argument is allowed');
end

h5att = h5info(fullfile(hdf5Dir,hdf5Name),'/Data/Images');
imageStackSize = h5att.Dataspace.Size;
rows = imageStackSize(1);
cols = imageStackSize(2);
num_frames = imageStackSize(3);

h5create(fullfile(hdf5Dir,hdf5Name),'/Data/Norm',[rows,cols,num_frames],'ChunkSize',[rows,cols,1],'Datatype','single');

% Apply spatial normalization
hDisk = fspecial('disk', disk_radius);

for i=1:num_frames
    if (mod(i,1000)==0)
        fprintf('  Frames %d of %d done\n', i, num_frames);
    end
    m = single(h5read(fullfile(hdf5Dir,hdf5Name),dataset,[1,1,i],[rows,cols,1]));
    m_f = imfilter(m, hDisk, 'replicate');
    M_norm = m./m_f;
    h5write(fullfile(hdf5Dir,hdf5Name),'/Data/Norm',M_norm,[1,1,i],[rows,cols,1]);
end
