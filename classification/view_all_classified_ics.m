function view_all_classified_ics(class_file,ica_mat,label_option)
% View ICs on top of cell map, where color of outline indicates
%   cell (green) or not (cell)
%   
%   class_file = classification text file
%   ica_mat = ica mat file
%   label_option = 'y' if you want IC number labels; 'n' if not
%
% 2015 02 19 Fori Wang

figure;

% load ICs
load(ica_mat);

% Pull out classification data
class = load_classification(class_file);

% Draw background image (sum of all IC filters)
h = imagesc(sum(ica_filters,3));
colormap gray;
axis image;
xlabel('x [px]');
ylabel('y [px]');
title_pwd = strrep(pwd, '_','\_');
title_class_file = strrep(class_file, '_','\_');
title(['Classified ICs for ',title_pwd, ' from ',title_class_file])

hold on;
colors = ['g'; 'r'];

% Draw all filter outlines (green = cell; red = not cell)
for ic_idx = 1:length(class)
    % Generate the outline of the IC filter
    ic_filter_threshold = 0.3;
    ic_filter = ica_filters(:,:,ic_idx);
    B = threshold_ic_filter(ic_filter, ic_filter_threshold);
    boundaries = bwboundaries(B, 'noholes');

    % Setup colors
    if strcmp(class(ic_idx),'cell')
        line_color = colors(1);
    else
        line_color = colors(2);
    end
    
    % Draw filter outline, color depending on classification
    for i = 1:length(boundaries)
        boundary = boundaries{i};
        plot(boundary(:,2), boundary(:,1), line_color, 'LineWidth', 1);
    end
    
    % label outline with IC #
    if strcmp(label_option,'y')
        text(max(boundary(:,2)),max(boundary(:,1)),int2str(ic_idx),...
                'Color',line_color);
    end
end