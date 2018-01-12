function md = create_merge_md(ds_array)
% See also 'resolve_recs'.
%
% Usage:
%   md = create_merge_md([ds1 ds2 ds3]);

num_ds = length(ds_array);
ds_list = cell(num_ds, 2);
for k = 1:num_ds
    ds_list{k,1} = k;
    ds_list{k,2} = ds_array(k);
end

num_match = num_ds*(num_ds-1)/2;
match_list = cell(num_match, 4);
idx = 1;
for k = 1:(num_ds-1)
    for l = (k+1):num_ds
        [m_k2l, m_l2k] = run_alignment(ds_array(k), ds_array(l), 'notrans', 'noprompt');
        close all;
        match_list{idx,1} = k;
        match_list{idx,2} = l;
        match_list{idx,3} = m_k2l;
        match_list{idx,4} = m_l2k;
        idx = idx + 1;
    end
end
md = MultiDay(ds_list, match_list, 'keepall');
