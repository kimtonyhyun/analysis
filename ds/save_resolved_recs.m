function rec_savename = save_resolved_recs(res_list, md)

info.type = 'resolved';
info.num_pairs = md.num_cells;

% We assume that image sizes and trace lengths are all the same for each
% "day" of the provided MultiDay instance. TODO: assert this assumption.
ind = md.valid_days(1); % Dummy index
[height, width] = size(md.day(ind).cells(1).im);
num_frames = md.day(ind).full_num_frames;

filters = zeros(height, width, md.num_cells);
traces = zeros(num_frames, md.num_cells);

idx = 0;
for k = 1:md.num_cells
    res = res_list(k,:);
    if (res(1) > 0)
        idx = idx + 1;
        ds = md.day(res(1));
        filters(:,:,idx) = ds.cells(res(2)).im;
        traces(:,idx) = ds.get_trace(res(2));
    end
end
info.num_pairs = idx;
filters = filters(:,:,1:idx);
traces = traces(:,1:idx);

rec_savename = save_rec(info, filters, traces);
fprintf('Saved resolved rec to "%s"!\n', rec_savename);