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

match_list = cell(num_ds-1, 4);
for k = 1:(num_ds-1)
    [m_itoj, m_jtoi] = run_alignment(ds_array(k), ds_array(k+1), 'notrans', 'noprompt');
    close all;
    match_list{k,1} = k;
    match_list{k,2} = k+1;
    match_list{k,3} = m_itoj;
    match_list{k,4} = m_jtoi;
end
md = MultiDay(ds_list, match_list, 'keepall');
