function rsData = memResample(ftdata,rs_fs,NSplit)
% Memory Efficient Resample
datr = {};
Nblc = ceil(numel(ftdata.label)./NSplit);
p = 1;
for i=1:Nblc
    list = p:p+NSplit-1;
    if list(end)>numel(ftdata.label)
        list = p:numel(ftdata.label);
    end
    
    cfgp         = [];
    cfgp.channel = [p:(p+(NSplit-1))];
    datp         = ft_preprocessing(cfgp,ftdata);
    
    cfgr            = [];
    cfgr.resamplefs = rs_fs;
    datr{i}         = ft_resampledata(cfgr, datp);
    clear datp
    p = p+NSplit;
end

cfg = [];
rsData = ft_appenddata(cfg, datr{:}); % this expands all cells into input variables
clear datr;