function ftdata_cat = appendFTData(ftdata_rep,repN,block)
flag = 1; % signals first rep to initialize
cfg.keepsampleinfo  = 'no';
cfg.appenddim = 'rpt';
if nargin == 3
    hdr = ftdata_rep{repN(1),block}.hdr;
    for rep = repN
        if flag % do once
            ftdata_cat = ftdata_rep{rep,block};
            flag = 0; % flag down
        else
            ftdata_cat = ft_appenddata(cfg,ftdata_cat,ftdata_rep{rep,block});
        end
    end
    
elseif nargin<3
    hdr = ftdata_rep{repN(1)}.hdr;
    
    for  rep = 1:repN
        if flag % do once
            ftdata_cat = ftdata_rep{rep};
            flag = 0; % flag down
        else
            if rep ==36
                a = 1;
            end
            ftdata_cat = ft_appenddata(cfg,ftdata_cat, ftdata_rep{rep});
        end
    end
end



ftdata_cat.hdr = hdr;



%%% SCRIPT GRAVE
% if repN == 1
%     ftdata_cat = ftdata_rep{1,block};
% elseif repN == 2
%     ftdata_cat = ft_appenddata([],ftdata_rep{1,block},ftdata_rep{2,block});
% elseif repN == 3
%     ftdata_cat = ft_appenddata([],ftdata_rep{1,block},ftdata_rep{2,block},ftdata_rep{3,block});
% elseif repN == 4
%     ftdata_cat = ft_appenddata([],ftdata_rep{1,block},ftdata_rep{2,block},ftdata_rep{3,block},ftdata_rep{4,block});
% end