function [X,cell_idx,trial_idx] = multiday_tensor(sessions,trial_type)
% MULTIDAY_TENSOR, construct a 3d tensor (neurons x time x trial) from
% several behavioral sessions.
%
% Example usage:
%   m1d12 = Day
%   ds_list = {12, m1d12; 13, m1d13; 14, m1d14};


if nargin == 1
    trial_type = 'all';
end

n_sessions = size(ds_list,1);
