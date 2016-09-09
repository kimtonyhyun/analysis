% replace with your path, or comment this out and load data yourself
if ismac
    p_m1d12 = '/Users/alex/Dropbox/strategy_switch/data/c11m1d12';
    p_m1d13 = '/Users/alex/Dropbox/strategy_switch/data/c11m1d13';
    p_m1d14 = '/Users/alex/Dropbox/strategy_switch/data/c11m1d14';
    figdir = '/Users/alex/Dropbox/strategy_switch/figs/';
    load /Users/alex/Dropbox/strategy_switch/data/match_m1d12-14.mat
elseif isunix
    p_m1d12 = '/home/alex/Dropbox/strategy_switch/data/c11m1d12';
    p_m1d13 = '/home/alex/Dropbox/strategy_switch/data/c11m1d13';
    p_m1d14 = '/home/alex/Dropbox/strategy_switch/data/c11m1d14';
    figdir = '/home/alex/Dropbox/strategy_switch/figs/';
    load /home/alex/Dropbox/strategy_switch/data/match_m1d12-14.mat
else
    error('Platform not supported')
end

src_m1d12 = gen_data_sources(p_m1d12);
m1d12 = DaySummary(src_m1d12, [p_m1d12,'/cm01']);
          
md = MultiDay({12,m1d12},{});

save('alex_12','md','-v7.3')
clearvars -except md
