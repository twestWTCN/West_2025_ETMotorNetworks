function [dataX_shift, dataY_shift, sampShift] = realignData(X, FX, Y, FY, fsamp, dataX, dataY, plotop, tarInd)
    % REALIGNDATAV2: Shifts vectors F(X) and F(Y) by introducing a NaN padding to maximize their cross-correlation.
    % Parameters:
    %   X, Y     - Time vectors for signals FX and FY respectively.
    %   FX, FY   - Signals to be aligned.
    %   fsamp    - Sampling frequency.
    %   dataX, dataY - Data arrays to be shifted accordingly.
    %   plotop   - (Optional) Set to 1 to plot results, 0 otherwise (default: 0).
    %   tarInd   - (Optional) Indices of channels to plot for verification.
    % Returns:
    %   dataX_shift, dataY_shift - Shifted versions of input data.
    %   sampShift - Number of samples that data has been shifted by.

    if nargin < 8
        plotop = 0;
    end

    %% Step 1: Initialize Plotting and Compute Cross-Correlation
    flag = 0; % Flag if variables need to be switched

    if plotop == 1
        figure;
        subplot(2, 2, 1);
        plot(X, FX); hold on;
        plot(Y, FY);
        title('Original Signals');
    end

    % Compute Cross-Correlation
    N = numel(FX) - 1;
    [xc, clag] = crosscorr(FX, FY, N);
    [a, loc] = max(abs(xc));

    if plotop == 1
        subplot(2, 2, 2);
        plot(clag, xc);
        title('Cross-Correlation');
    end

    %% Step 2: Switch Variables if Necessary
    if clag(loc) < 0
        % Switch X and Y if lag is negative (X should lead)
        [X, Y] = switchVariables(X, Y);
        [FX, FY] = switchVariables(FX, FY);

        % Recompute Cross-Correlation
        N = min([numel(FX) - 1, numel(FY) - 1]);
        [xc, clag] = crosscorr(FX, FY, N);
        [a, loc] = max(abs(xc));

        flag = 1; % Set flag indicating switch
        varN = 1; % Indicates X is being shifted
    else
        varN = 2; % Indicates Y is being shifted (default behavior)
    end

    %% Step 3: Handle Already-Aligned Signals
    if clag(loc) == 0
        warning('Signals already appear to be aligned!');
        sampShift = clag(loc);
        dataX_shift = [];
        dataY_shift = [];
        return;
    end

    %% Step 4: Shift Y and Define New Time Vectors
    Y_Shift = Y - (clag(loc) / fsamp); % Adjust time axis
    FY_Shift = FY(clag(loc):end);
    Y_Shift = Y_Shift(1:numel(FY_Shift));

    % Define Shift Amount
    if varN == 1
        sampShift = -clag(loc);
    elseif varN == 2
        sampShift = clag(loc);
    end

    % Truncate FX
    FX_Shift = FX(1:numel(Y_Shift));
    X_Shift = Y_Shift;

    %% Step 5: Plot Shifted Signals
    if plotop == 1
        subplot(2, 2, 3);
        plot(X_Shift, FX_Shift); hold on;
        plot(Y_Shift, FY_Shift);
        title('Aligned Signals');
    end

    % Recompute Cross-Correlation to Verify Alignment
    [xc, clag2] = crosscorr(FY_Shift, FX_Shift, numel(FY_Shift) - 1);
    [a, loc2] = max(abs(xc));

    if plotop == 1
        subplot(2, 2, 4);
        plot(clag2 ./ fsamp, xc);
        title('Cross-Correlation After Alignment');
    end

    if abs(clag2(loc2) / fsamp) > 1 / fsamp
        warning('Realignment failed! Significant residual lag detected.');
    end

    %% Step 6: Shift Data Arrays (dataX, dataY)
    if ~isempty(dataX) && ~isempty(dataY)
        if varN == 1
            % Shift dataX and adjust dataY accordingly
            dataX_shift = dataX(:, clag(loc):end);
            dataY_shift = dataY(:, 1:size(FY_Shift, 2));
        elseif varN == 2
            % Shift dataY and adjust dataX accordingly
            dataY_shift = dataY(:, clag(loc) + 1:end);
            dataX_shift = dataX(:, 1:size(dataY_shift, 2));
            X_Shift = Y_Shift;
        end

        %% Step 7: Plot Verification (if plotop is enabled)
        if plotop
            figure;
            subplot(2, 1, 1);
            plot(X_Shift, dataX_shift(tarInd(1), :));
            yyaxis right;
            plot(Y_Shift, dataY_shift(tarInd(2), :));
            title('Shifted Data Alignment');

            % Recompute Cross-Correlation for Verification
            [xc, clag2] = crosscorr(dataX_shift(tarInd(1), :), dataY_shift(tarInd(2), :), size(dataX_shift, 2) - 1);
            subplot(2, 1, 2);
            plot(clag2 ./ fsamp, xc);
            title('Cross-Correlation of Shifted Data');
        end

        % Print Realignment Status
        if abs(clag2(loc2) / fsamp) > 0.05
            warning('Realignment failed!');
        else
            fprintf('Realignment success: error of %.3f ms is within accepted limit\n', abs(clag2(loc2) / fsamp) * 1000);
        end
    else
        dataX_shift = [];
        dataY_shift = [];
    end
end

%% Helper Function to Switch Variables
function [var1, var2] = switchVariables(var1, var2)
    temp = var1;
    var1 = var2;
    var2 = temp;
end
