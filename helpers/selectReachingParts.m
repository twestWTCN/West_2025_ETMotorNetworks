function ftdata_cat = selectReachingParts(ftdata_cat,part)
% This selects the right epochs from trans data to ensure you get the right
% side.
% % switch part
% %     case 1 % Rest
% %         cfg = [];
% %         cfg.latency = [1.5 3]; % Take pre posture
% %         ftdata_cat = ft_selectdata(cfg,ftdata_cat);
% %     case 2 % Posture
% %         cfg = [];
% %         cfg.latency = [0.75 3.75]; % Take post posture
% %         ftdata_cat = ft_selectdata(cfg,ftdata_cat);
% %     case 3 % Prep
% %         cfg = [];
% %         cfg.latency = [0 2]; % Take post cue
% %         ftdata_cat = ft_selectdata(cfg,ftdata_cat);
% %     case 4 % Reach
% %         cfg = [];
% %         cfg.latency = [-0.5 1]; % Take post cue
% %         ftdata_cat = ft_selectdata(cfg,ftdata_cat);
% %     case 5 % Hold
% %         cfg = [];
% %         cfg.latency = [0.3 2]; % Take post cue
% %         ftdata_cat = ft_selectdata(cfg,ftdata_cat);
% % end

switch part
    case 1 % Rest
        cfg = [];
        cfg.latency = [-2 0]; % Take pre posture
        ftdata_cat = ft_selectdata(cfg,ftdata_cat);
    case 2 % Posture
        cfg = [];
        cfg.latency = [-0.5 1.5]; % Take post posture
        ftdata_cat = ft_selectdata(cfg,ftdata_cat);
    case 3 % Prep
        cfg = [];
        cfg.latency = [0 2]; % Take post cue
        ftdata_cat = ft_selectdata(cfg,ftdata_cat);
    case 4 % Reach
        cfg = [];
        cfg.latency = [-0.5 1.5]; % Take post exec init
        ftdata_cat = ft_selectdata(cfg,ftdata_cat);
    case 5 % Hold
        cfg = [];
        cfg.latency = [0.5 1.5]; % Take post exec term
        ftdata_cat = ft_selectdata(cfg,ftdata_cat);
end
