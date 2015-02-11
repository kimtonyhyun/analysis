function view_ic_pair(time, trace, filter)
% Plot IC trace and IC filter (as a 3D graph)
%
% 2015 02 08 Tony Hyun Kim

subplot(3,1,[2 3]);
surf(double(filter)); % 'surf' doesn't take single
colormap jet;
shading flat;
[h, w] = size(filter);
xlim([1 w]);
ylim([1 h]);
z_min = min(filter(:));
z_max = max(filter(:));
z_range = z_max - z_min;
zlim([z_min z_max] + 0.1*z_range*[-1 1]);
daspect([1 1 z_range/max(h,w)]);
xlabel('x [px]');
ylabel('y [px]');
set(gca, 'YDir', 'Reverse');
zlabel('IC filter [a.u.]');

subplot(3,1,1);
plot(time, trace);
t_min = min(trace);
t_max = max(trace);
t_range = t_max - t_min;
xlim([time(1) time(end)]);
ylim([t_min t_max] + 0.1*t_range*[-1 1]);
xlabel('Time [s]');
ylabel('IC trace [a.u.]');