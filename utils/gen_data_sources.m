function sources = gen_data_sources()

datafiles = dir('_data');

for i = 1:length(datafiles)
	
	fname = datafiles(i).name;
	[~,name,ext] = fileparts(fname);
	
	if strcmp(ext,'.txt')
		sources.maze = ['_data/',fname];
	end
	if strcmp(ext,'.mp4')
		%sources.behavior = ['_data/',fname];
	end
	if strcmp(ext,'.xy')
		sources.tracking = ['_data/',fname];
	end
	if strcmp(ext,'.hdf5')
		sources.miniscope = ['_data/',fname];
	end

end

end
