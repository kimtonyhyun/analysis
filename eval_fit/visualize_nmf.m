function visualize_nmf(T,k,trial_meta)

Tk = T(:,:,k)'; % [Trials x Time]
[num_trials, num_samples] = size(Tk);

figure;
subplot(3,4,[1 2]);
hold on;
for i = 1:num_trials
    plot(Tk(i,:), '.');
end
xlim([1 num_samples]);
title(sprintf('Factor k=%d', k));

trace_avg = mean(Tk,1);
plot(trace_avg, 'r', 'LineWidth', 3);
grid on;
xlabel('Within-trial sample');
ylabel('Factor amplitude');

subplot(3,4,[5 6 9 10]);
imagesc(Tk);
xlim([1 num_samples]);
ylabel([1 num_trials]);
xlabel('Within-trial sample');
ylabel('Trial index');

% Separate trials by start location
filtered_trials = (trial_meta.correct == 1);
subplot(3,4,3);
imagesc(Tk(filtered_trials,:));
title('Correct trials');

filtered_trials = ~filtered_trials;
subplot(3,4,4);
imagesc(Tk(filtered_trials,:));
title('Incorrect trials');

% Separate trials by start location
filtered_trials = cellfun(@(x) strcmp(x, 'east'), trial_meta.start);
subplot(3,4,7);
imagesc(Tk(filtered_trials,:));
title('East start');

filtered_trials = cellfun(@(x) strcmp(x, 'west'), trial_meta.start);
subplot(3,4,8);
imagesc(Tk(filtered_trials,:));
title('West start');

% Separate trials by end location
filtered_trials = cellfun(@(x) strcmp(x, 'north'), trial_meta.end);
subplot(3,4,11);
imagesc(Tk(filtered_trials,:));
title('North end');

filtered_trials = cellfun(@(x) strcmp(x, 'south'), trial_meta.end);
subplot(3,4,12);
imagesc(Tk(filtered_trials,:));
title('South end');

end % visualize_nmf

