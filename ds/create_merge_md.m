function md = create_merge_md(ds_array)
% Helper function to create MultiDay instance to represent multiple
% extraction runs on the same dataset, to be merged into one set. Used, for
% example, with iterative CNMF or creating union filter sets across days.
%
% Input 'ds_array' may be an array of DaySummary's (e.g. [ds1 ds2 ds3]) or
% a cell array (e.g. {ds1, ds2, ds3}).
%
% Usage:
%   md = create_merge_md([ds1 ds2 ds3]);
%   res_list = resolve_merged_recs(md, M); % 'M' is the fluorescence movie
%   save_resolved_recs(res_list, md);
%
if iscell(ds_array)
    get_kth_ds = @(k) ds_array{k};
else
    get_kth_ds = @(k) ds_array(k);
end
num_ds = length(ds_array);
ds_list = cell(num_ds, 2);
for k = 1:num_ds
    ds = get_kth_ds(k);
    if all(cellfun(@isempty, ds.get_class))
        fprintf('create_merge_md: DS%d is completely unlabeled. Setting all sources to be cells!\n', k);
        ds.set_labels;
    end
    
    ds_list{k,1} = k;
    ds_list{k,2} = ds;
end

match_list = cell(num_ds*(num_ds-1)/2, 4);
idx = 1;
for i = 1:(num_ds-1)
    for j = (i+1):num_ds
        [m_itoj, m_jtoi] = run_alignment(get_kth_ds(i), get_kth_ds(j), 'notrans', 'noprompt');
        close all;
        match_list{idx,1} = i;
        match_list{idx,2} = j;
        match_list{idx,3} = m_itoj;
        match_list{idx,4} = m_jtoi;
        idx = idx + 1;
    end
end
md = MultiDay(ds_list, match_list, 'keepall');
