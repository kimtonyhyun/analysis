function visualize_nmf(T, k)

Tk = T(:,:,k);
[num_samples, num_trials] = size(Tk);

figure;
subplot(121);
hold on;
for i = 1:num_trials
    plot(Tk(:,i), '.');
end
xlim([1 num_samples]);

trace_avg = mean(Tk,2);
plot(trace_avg, 'r', 'LineWidth', 3);
grid on;
xlabel('Within-trial sample');
ylabel('Factor amplitude');

subplot(122);
imagesc(Tk');
xlim([1 num_samples]);
ylabel([1 num_trials]);
xlabel('Within-trial sample');
ylabel('Trial index');

suptitle(sprintf('Factor k=%d', k));