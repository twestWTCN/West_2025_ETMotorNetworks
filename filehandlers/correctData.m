function correctData(R,cordat)
ip = 0; varPartFlag = 0;
for sub = cordat.sub
    for block = cordat.block % Block of the Experiment
        switch R.import.blockpart{block}
            case 'Task'
                for cond = cordat.cond % Condion type
                    for part = cordat.part % part of task
                        ip = ip + 1;
                        fileappend = eval(cordat.fileapend);
                        ftdata_tmp = loadExpData(R,sub{1},'Task',[],fileappend,eval(cordat.subfold));
                        if strcmp(cordat.vargin{end},'part') || varPartFlag
                            varPartFlag = 1;
                            cordat.vargin{end} = part;
                        end
                        ftdata_tmp = cordat.corrfx(ftdata_tmp,cordat.vargin);
%                         fileappend = [fileappend '_CM'];
                        saveExpData(R,sub{1},'Task',[],fileappend,ftdata_tmp,eval(cordat.subfold));
                    end
                end
            otherwise
                ip = ip + 1;
                fileappend = eval(cordat.fileapend);
                ftdata_tmp = loadExpData(R,sub{1},R.import.blockpart{block},[],fileappend,cordat.subfold);
                
                ftdata_tmp = cordat.corrfx(ftdata_tmp,cordat.vargin);
%                 fileappend = [fileappend '_CM'];
                saveExpData(R,sub{1},R.import.blockpart{block},[],fileappend,ftdata_tmp,cordat.subfold);
        end
    end
end
