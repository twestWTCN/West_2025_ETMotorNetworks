function [] = correctMovementAxes(R)
% CORRECTMOVEMENTAXES ensures that movement axes are correctly labelled.
% Corrects the axes for consistency across subjects and conditions.
% X = L/R; Y = Down/Up; Z = Near/Far.

subI = 0;
for sub = R.corrMoveAxes.subsel
    subI = subI + 1;

    % Get normalization factor for everything
    retdat.sub = sub;
    retdat.block = 3; % task
    retdat.cond = 1:4; % condition
    retdat.part = 2:5; % movement
    retdat.fileapend = "[R.epoch.names{part} '_trans_pp_arV3']";
    retdat.subfold = "['Condition' num2str(cond) '_epoched']";
    retdat.ftflag = 1; % nonFT concat
    ftdata_cat = retrieveData(R, retdat);

    motNames = {'MotivePos', 'motivePos'};
    for mt = 1:2
        % Compute factor
        if mt == 1
            motind = [find(strcmp(ftdata_cat.label, 'PosX')),...
                      find(strcmp(ftdata_cat.label, 'PosY')),...
                      find(strcmp(ftdata_cat.label, 'PosZ'))];
        elseif mt == 2
            motind = [find(strcmp(ftdata_cat.label, 'xFinger1Pos')),...
                      find(strcmp(ftdata_cat.label, 'yFinger1Pos')),...
                      find(strcmp(ftdata_cat.label, 'zFinger1Pos'))];
        end

        if ~isempty(motind)
            cr = {'', 'r', 'b', 'g', 'k'};
            for prt = 2:5
                % Determine the Y axis from Posture - predominant positive Y movement
                retdat.part = prt; % movement
                ftdata_cat = retrieveData(R, retdat);

                if mt == 1
                    motind = [find(strcmp(ftdata_cat.label, 'PosX')),...
                              find(strcmp(ftdata_cat.label, 'PosY')),...
                              find(strcmp(ftdata_cat.label, 'PosZ'))];
                elseif mt == 2
                    motind = [find(strcmp(ftdata_cat.label, 'xFinger1Pos')),...
                              find(strcmp(ftdata_cat.label, 'yFinger1Pos')),...
                              find(strcmp(ftdata_cat.label, 'zFinger1Pos'))];
                end

                accind = find(strncmp(ftdata_cat.label, 'Acc', 3));

                % Normalize the signals
                cfg = [];
                cfg.channel = ftdata_cat.label(motind);
                motdata = ft_selectdata(cfg, ftdata_cat);
                motdataPrt = selectReachingParts(motdata, prt);

                cfg = [];
                cfg.channel = ftdata_cat.label(accind);
                accdata = ft_selectdata(cfg, ftdata_cat);
                accdataPrt = selectReachingParts(accdata, prt);

                XTMot = []; XTAcc = []; Xmot = []; Xacc = [];
                for tr = 1:numel(motdataPrt.trial)
                    XTMot(:, :, tr) = motdataPrt.trial{tr};
                    XTAcc(:, :, tr) = cumsum(accdataPrt.trial{tr});
                    Xmot(:, tr) = range(motdataPrt.trial{tr}, 2)';
                    Xacc(:, tr) = range(accdataPrt.trial{tr}, 2)';
                end

                figure(10 + mt);
                X = squeeze(mean(XTMot(1, :, :), 3));
                Y = squeeze(mean(XTMot(2, :, :), 3));
                Z = squeeze(mean(XTMot(3, :, :), 3));
                plot3(X, Y, Z, cr{prt}); hold on;
                scatter3(X(1), Y(1), Z(1), 100, cr{prt}, 'o');
                scatter3(X(end), Y(end), Z(end), 100, cr{prt}, 'x');
                xlabel(motdata.label{1}); ylabel(motdata.label{2}); zlabel(motdata.label{3});
                hold on;
                title('MOT');
                axis equal;
                if prt == 2 || prt == 5
                    text(X(1), Y(1), Z(1), sub{1}(1:4));
                end
                
                figure(2);
                X = squeeze(mean(XTAcc(1, :, :), 3));
                Y = squeeze(mean(XTAcc(2, :, :), 3));
                Z = squeeze(mean(XTAcc(3, :, :), 3));
                plot3(X, Y, Z, cr{prt}); hold on;
                scatter3(X(1), Y(1), Z(1), 100, cr{prt}, 'o');
                scatter3(X(end), Y(end), Z(end), 100, cr{prt}, 'x');
                xlabel('X'); ylabel('Y'); zlabel('Z');
                hold on;
                title('ACC');
                axis equal;

                if prt == 4
                    % Get Motive Ranges
                    [~, maxReachMotAx] = max(mean(Xmot, 2));
                    maxReachMotdelta = squeeze(mean(XTMot(maxReachMotAx, :, :), 3));
                    maxReachMotdelta = diff(maxReachMotdelta([1, end]));

                    % Get ACC Ranges
                    [~, maxReachAccAx] = max(mean(Xacc, 2));
                    maxReachAccdelta = squeeze(mean(XTAcc(maxReachAccAx, :, :), 3));
                    maxReachAccdelta = diff(maxReachAccdelta([1, end]));
                elseif prt == 5
                    trtar = ftdata_cat.trialinfo(:, 2);
                    XAx = setdiff(1:3, [maxReachMotAx, maxPostMotAx]);

                    % Left and Right targets
                    trLeft = find(trtar == 45 | trtar == 360 | trtar == 315);
                    trRight = find(trtar == 135 | trtar == 180 | trtar == 225);

                    % Get Positive X direction
                    maxRightHoldMotDelta = diff([squeeze(mean(mean(XTMot(XAx, :, trLeft))))', squeeze(mean(mean(XTMot(XAx, :, trRight))))']);
                elseif prt == 2
                    % Get Posture Ranges
                    [~, maxPostMotAx] = max(mean(Xmot, 2));
                    maxPostMotdelta = squeeze(mean(XTMot(maxPostMotAx, :, :), 3));
                    maxPostMotdelta = diff(maxPostMotdelta([1, end]));

                    % Get ACC Ranges
                    [~, maxPostAccAx] = max(mean(Xacc, 2));
                    maxPostAccdelta = squeeze(mean(XTAcc(maxPostAccAx, :, :), 3));
                    maxPostAccdelta = diff(maxPostAccdelta([1, end]));
                end
            end

            % Redefine the axes
            ZMot = maxReachMotAx;
            YMot = maxPostMotAx;
            XMot = setdiff(1:3, [ZMot, YMot]);

            ZFlip = double(maxReachMotdelta < 0);
            YFlip = double(maxPostMotdelta < 0);
            XFlip = double(maxRightHoldMotDelta < 0);

            if any(strcmp(motdata.label(XMot)', {'xFinger1Pos', 'PosX'})) && ...
               any(strcmp(motdata.label(YMot)', {'yFinger1Pos', 'PosY'})) && ...
               any(strcmp(motdata.label(ZMot)', {'zFinger1Pos', 'PosZ'}))
                disp('Subjects motive axes are properly aligned');
            else
                % Correct Steady State Data
                corDat = [];
                corDat.sub = sub;
                corDat.block = 1:2; % task
                corDat.fileapend = "'pp_arV3'";
                corDat.subfold = 'epoched';
                corDat.corrfx = @relabelMovementAxes;
                corDat.vargin = {[], [XMot, YMot, ZMot], [XFlip, YFlip, ZFlip]};
                if all(ftdata_cat.label{motind(1)}(1:3) == 'Pos')
                    corDat.vargin{1} = {'PosX', 'PosY', 'PosZ'};
                end
                if all(ftdata_cat.label{motind(1)}(end-2:end) == 'Pos')
                    corDat.vargin{1} = {'xFinger1Pos', 'yFinger1Pos', 'zFinger1Pos'};
                end
                correctData(R, corDat);

                % Correct Blocked Task Data
                corDat = [];
                corDat.sub = sub;
                corDat.block = 3; % task
                corDat.cond = 1:4; % condition
                corDat.part = 1:5; % movement
                corDat.fileapend = "[R.epoch.names{part} '_block_pp_arV3']";
                corDat.subfold = "['Condition' num2str(cond) '_epoched']";
                corDat.corrfx = @relabelMovementAxes;
                corDat.vargin = {[], [XMot, YMot, ZMot], [XFlip, YFlip, ZFlip]};
                if all(ftdata_cat.label{motind(1)}(1:3) == 'Pos')
                    corDat.vargin{1} = {'PosX', 'PosY', 'PosZ'};
                end
                if all(ftdata_cat.label{motind(1)}(end-2:end) == 'Pos')
                    corDat.vargin{1} = {'xFinger1Pos', 'yFinger1Pos', 'zFinger1Pos'};
                end
                correctData(R, corDat);

                % Correct Time Locked Task Data
                corDat = [];
                corDat.sub = sub;
                corDat.block = 3; % task
                corDat.cond = 1:4; % condition
                corDat.part = 2:5; % movement
                corDat.fileapend = "[R.epoch.names{part} '_trans_pp_arV3']";
                corDat.subfold = "['Condition' num2str(cond) '_epoched']";
                corDat.corrfx = @relabelMovementAxes;
                corDat.vargin = {[], [XMot, YMot, ZMot], [XFlip, YFlip, ZFlip]};
                if all(ftdata_cat.label{motind(1)}(1:3) == 'Pos')
                    corDat.vargin{1} = {'PosX', 'PosY', 'PosZ'};
                end
                if all(ftdata_cat.label{motind(1)}(end-2:end) == 'Pos')
                    corDat.vargin{1} = {'xFinger1Pos', 'yFinger1Pos', 'zFinger1Pos'};
                end
                correctData(R, corDat);
            end
        end
    end
end % Sub loop

end