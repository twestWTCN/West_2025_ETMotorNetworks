function [nantr,list] = ft_nancheck(ftdata)
list = []; nantr = [];
for tr = 1:numel(ftdata.trial)
    trX = ftdata.trial{tr};
    if ~any(isnan(trX(:)))
        list = [list tr];
    else
        nantr = [nantr tr];
    end
end

if any(nantr)
    warning('Your data structure contains NaNs!!!')
    a = [];
end
