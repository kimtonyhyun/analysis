function outfile = save_tdt(ds, outfile)

if ~exist('outfile', 'var')
    outfile = sprintf('tdt_%s.mat', datestr(now, 'yymmdd-HHMMSS'));
end

[pos, neg] = collect_tdt_labels(ds); %#ok<ASGLU>
save(outfile, 'pos', 'neg');

fprintf('%s: tdTomato labels saved to "%s"\n', datestr(now), outfile);