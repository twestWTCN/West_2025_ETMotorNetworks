function ftdata = addHistoryField(ftdata, entry, prsphist)
    % ADDHISTORYFIELD: Adds a history entry to the FieldTrip data structure.
    % Ensures that a history field exists and appends the provided entry to it.
    % Parameters:
    %   ftdata - FieldTrip data structure
    %   entry - Entry to add to the history field
    %   prsphist - (Optional) Pre-specified history to set

    %% Step 1: Initialize History Field if Missing
    if ~isfield(ftdata, 'hdr')
        ftdata.hdr = [];
    end

    if ~isfield(ftdata.hdr, 'history')
        ftdata.hdr.history = {};
    end

    %% Step 2: Set History to Pre-specified if Provided
    if nargin > 2
        ftdata.hdr.history = prsphist;
    end

    %% Step 3: Remove Existing Entry if Duplicate
    histlist = strncmp(ftdata.hdr.history, entry, 4);
    if any(histlist)
        ftdata.hdr.history(histlist) = [];
    end

    %% Step 4: Add New Entry with Current Date
    ftdata.hdr.history = [ftdata.hdr.history, [entry ' ' date]];
end
