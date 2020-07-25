function resolve_depths(ds1, ds2)
% Uses the 'browse_corrlist' mechanism underneath

corrlist = compute_corrlist(ds1, ds2);
browse_corrlist(corrlist, ds1, ds2,...
                'app_name', 'Resolve depths',...
                'zsc', 'overlay');