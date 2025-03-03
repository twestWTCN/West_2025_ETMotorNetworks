function ftdata_flipped = flipSensors(ftdata)
% This function flips sensor labels to switch between right and left
% dominant
% Acticap layout convention: labels containing an odd/even number are
% left/right
ftdata_flipped = ftdata;
% Flip the main channel labels
if isfield(ftdata,'hdr')
    eegSel = find(strcmp(ftdata.hdr.chantype,'eeg') |  strcmp(ftdata.hdr.chantype,'EEG') |  strcmp(ftdata.hdr.chantype,'vc') );
    labelsFlip = eegLabelFlip(ftdata.label(eegSel));
    ftdata_flipped.label(eegSel) = labelsFlip;
    display('Fliping data labels!')
end
% Flip also elec definition
if isfield(ftdata,'elec')
    ftdata_flipped.elec.label =  eegLabelFlip(ftdata_flipped.elec.label);
    display('Fliping elec labels!')
end

if isfield(ftdata,'powspctrm')
    ftdata_flipped.label =  eegLabelFlip(ftdata.label);
    display('Fliping powspctrm labels!')
end


function labelsFlip = eegLabelFlip(labels)
if any(startsWith(labels,'Cerebellum')) % then VC
    for ch = 1:numel(labels)
        chtmp = labels{ch};
        leftchan = regexp(chtmp(2:end) ,'L'); % start from 2nd to avoid capital
        rightchan = regexp(chtmp(2:end),'R');
        
        if any(leftchan)
            chtmp(leftchan+1) = 'R';
        end
        
        if any(rightchan)
            chtmp(rightchan+1) = 'L';
        end
        labelsFlip{ch,1} = chtmp;
    end
else % normal 1020
    
    
    oddInd =[]; evenInd = [];
    for ch = 1:numel(labels)
        chtmp = labels{ch};
        
        oddnum = regexp(chtmp ,'[1 3 5 7 9]'); % left channel
        evennum = regexp(chtmp ,'[2 4 6 8 10]'); % right channel
        
        % account for "10" comprising two digits
        if numel(evennum) ==2
            evenum = evennum(1); oddnum = [];
        end
        
        % Now flip
        if any(oddnum)
            oddInd = [oddInd ch];
            if numel(chtmp)>oddnum
                chtmp = [chtmp(1:oddnum(1)-1)  num2str(str2double(chtmp(oddnum))+1) chtmp(oddnum(1)+1:end)];
            else
                chtmp = [chtmp(1:oddnum(1)-1)  num2str(str2double(chtmp(oddnum))+1)];
            end
        elseif any(evennum)
            evenInd = [evenInd ch];
            if numel(chtmp)>evennum
                if str2double(chtmp(evennum))==10
                    chtmp = [chtmp(1:evennum(1)-1) '09' chtmp(evennum(1)+2:end)];
                else
                    chtmp = [chtmp(1:evennum(1)-1)  num2str(str2double(chtmp(evennum))-1) chtmp(evennum(1)+1:end)];
                end
            else
                chtmp = [chtmp(1:evennum(1)-1)  num2str(str2double(chtmp(evennum))-1)];
            end
        else
            chtmp = chtmp; % dont flip as midline
        end
        labelsFlip{ch,1} = chtmp;
    end
end
