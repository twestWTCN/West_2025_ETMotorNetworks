function AVData = readGermanAVData(filename)
% This interprets the German AV data - tricky due to comma delimited +
% comma decimal. Interpreted case by case.
load(filename,'Data');
Data = table2array(Data);

for i = 2:size(Data,1)-1
    XStr = Data(i,1);
    comLoc = strfind(XStr,',');
    comN = numel(comLoc);
    if comN == 13
        startp = [1 comLoc([1 2 4 6 8 10 12])+1]';
        endp = [comLoc([1 2 4 6 8 10 12])-1 strlength(XStr)]';
    elseif comN == 14
        startp = [1 comLoc([1 3 5 7 9 11 13])+1]';
        endp = [comLoc([1 3 5 7 9 11 13])-1 strlength(XStr)]';
    elseif comN == 12
        startp = [1 comLoc([1 2 3 5 7 9 11])+1]';
        endp = [comLoc([1 2 3 5 7 9 11])-1 strlength(XStr)]';
    else
        a = 1;
        error('Case not defined!')
    end
    
    for col = 1:numel(endp)
        tab(col,i-1) = str2double( regexprep(extractBetween(XStr,startp(col),endp(col)), ',', '.') );
    end
    
    % weird conditions
    if tab(3,i-1) > 1e13
        % odd case when time vector has no decimal
        startp = [1 comLoc([1 3 4 6 8 10 12])+1]';
        endp = [comLoc([1 3 4 6 8 10 12])-1 strlength(XStr)]';
        for col = 1:numel(endp)
            tab(col,i) = str2double( regexprep(extractBetween(XStr,startp(col),endp(col)), ',', '.') );
        end
    end
    
    if rem(i,1000) == 0
       disp((i/ (size(Data,1)-1))*100)
    end
end

AVData = array2table(tab','VariableNames',{'Condition','Coder','Timer','xScreenPosition','yScreenPosition','trackingX','trackingY','trackingZ'});
