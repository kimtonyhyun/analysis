function view_all_classified_ics(ic_dir, label_option)
% View ICs on top of cell map, where color of outline indicates
%   cell (green) or not (cell)
%   
%   class_file = classification text file
%   ica_mat = ica mat file
%   label_option = 'y' if you want IC number labels; 'n' if not
%
% 2015 02 19 Fori Wang

% Load data
ica_file = get_most_recent_file(ic_dir, 'ica_*.mat');
class_file = get_most_recent_file(ic_dir, 'class_*.txt');

load(ica_file); % 'ica_info', 'ica_traces', 'ica_filters'
class = load_classification(class_file);

% Draw background image (sum of all IC filters)
figure;
h = imagesc(sum(ica_filters,3));
colormap gray;
axis image;
xlabel('x [px]');
ylabel('y [px]');
title_pwd = strrep(pwd, '_','\_');
title_class_file = strrep(class_file, '_','\_');
title(['Classified ICs for ',title_pwd, ' from ',title_class_file])

hold on;

% Draw all filter outlines (green = cell; red = not cell)
ic_filter_threshold = 0.3;
for ic_idx = 1:length(class)
    ic_filter = ica_filters(:,:,ic_idx);
    boundaries = compute_ic_boundary(ic_filter, ic_filter_threshold);
    boundary = boundaries{1}; % Longest boundary
    
    % Setup colors
    if strcmp(class(ic_idx),'not a cell')
        line_color = 'r';
    else
        line_color = 'g';
    end
    
    % Draw filter outline, color depending on classification
    plot(boundary(:,1), boundary(:,2), line_color, 'LineWidth', 1);
    
    % label outline with IC #
    if strcmp(label_option,'y')
        text(max(boundary(:,1)),max(boundary(:,2)),int2str(ic_idx),...
                'Color',line_color);
    end
end