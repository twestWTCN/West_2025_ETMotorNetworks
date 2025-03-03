function local_peak = findLocalSourcePeak(source, roi_center, radius,par)
    % source: FieldTrip source structure
    % roi_center: [x, y, z] coordinates of the center of the ROI
    % radius: radius of the ROI in the same units as the source.pos (usually cm)
 
    % Extract the source positions and the coherence values
    pos = source.pos;  % Nx3 matrix of [x, y, z] positions
    coherence = source.(par);  % replace 'avg.coh' with the correct field for your data
 
    % Initialize variables
    local_peak = struct('pos', [], 'value', -inf);
 
    % Iterate over each source point
    for i = 1:size(pos, 1)
        % Calculate the distance from the ROI center
        distance = norm(pos(i, :) - roi_center);
 
        % Check if the point is within the ROI
        if distance <= radius
            % Check if this point has a higher coherence value than the current max
            if coherence(i) > local_peak.value
                % Update the local peak information
                local_peak.pos = pos(i, :);
                local_peak.value = coherence(i);
            end
        end
    end
 
    if isempty(local_peak.pos)
        disp('No local peak found within the specified ROI.');
    else
        fprintf('Local peak found at [%f, %f, %f] with a coherence value of %f.\n', ...
                local_peak.pos(1), local_peak.pos(2), local_peak.pos(3), local_peak.value);
    end
end