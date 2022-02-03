function load_tdt(ds, source)

tdt = load(source);

ds.reset_labels;

for k = tdt.pos
    ds.cells(k).label = 'positive';
end

for k = tdt.neg
    ds.cells(k).label = 'negative';
end