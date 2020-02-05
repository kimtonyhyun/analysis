function import_extract(extract_struct)

filters = extract_struct.spatial_weights;
traces = extract_struct.temporal_weights;
info.num_pairs = size(filters, 3);
info.type = 'extract';

info.extract.info = extract_struct.info;
info.extract.config = extract_struct.config;

save_rec(info, filters, traces);