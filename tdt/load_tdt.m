function tdt = load_tdt(source)

[~, ~, ext] = fileparts(source);
if isempty(ext)
    % 'source' is path to folder containing tdt classification
    filename = get_most_recent_file(source, 'tdt_*.mat');
    
    if isempty(filename)
        cprintf('blue', 'No tdt classification exists at "%s"\n', source);
        tdt.pos = [];
        tdt.neg = [];
    else
        tdt = load(filename);
    end
else
    % Assume full specification of the MAT file
    tdt = load(source);
end