function val = is_trial_frame(frame_indices, test_frame)

trial_idx = find(test_frame >= frame_indices(:,1), 1, 'last');
val = test_frame <= frame_indices(trial_idx,2);