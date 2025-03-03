function [mean_velocity, max_velocity, pathlength, movement_duration, smoothness_index, hold_stability, magnitude_velocityMove, moveIndex, holdIndex] = computeKinematics(XYZ, fsamp, reachThreshold, holdThreshold)
% Compute kinematic features for movement analysis
% Inputs:
%   XYZ - 3D position data (3 x time)
%   fsamp - sampling frequency
%   reachThreshold - threshold for movement onset detection
%   holdThreshold - threshold for hold detection
% Outputs:
%   mean_velocity - average velocity during movement
%   max_velocity - maximum velocity during movement
%   pathlength - total path length of the movement
%   movement_duration - duration of movement in ms
%   smoothness_index - integrated squared jerk for smoothness assessment
%   hold_stability - variability during hold phase
%   magnitude_velocityMove - velocity magnitude during movement
%   moveIndex - indices representing the movement period
%   holdIndex - indices representing the hold period (if any)

% Smooth the position data to reduce noise
for i = 1:3
    XYZ(i, :) = smooth(XYZ(i, :), 256);
end

% Compute velocity as the change in position over time
velocity = diff(XYZ, 1, 2) * fsamp;
for i = 1:3
    velocity(i, :) = smooth(velocity(i, :), 64);
end

% Calculate the magnitude of velocity
magnitude_velocity = sqrt(sum(velocity.^2, 1));
magnitude_velocity_z = magnitude_velocity / std(magnitude_velocity);

% Detect movement onset
indexes = find(magnitude_velocity_z > reachThreshold);
indexes = SplitVec(indexes, 'consecutive');
index_length = cellfun(@length, indexes);
[~, indmax] = max(index_length);
onset_index = indexes{indmax}(1);

% Detect movement offset
offset_index = onset_index + find(magnitude_velocity_z(onset_index:end) < reachThreshold, 1, 'first');
if isempty(offset_index) || offset_index > size(magnitude_velocity, 2)
    offset_index = size(magnitude_velocity, 2) - 2;
end
moveIndex = onset_index:offset_index;

% Extract movement trajectory and velocity during movement
movement_trajectory = XYZ(:, onset_index:offset_index);
magnitude_velocityMove = magnitude_velocity(onset_index:offset_index);

% Calculate path length
pathlength = sum(sqrt(sum(diff(movement_trajectory, 1, 2).^2)));

% Calculate maximum and mean velocity
max_velocity = max(magnitude_velocityMove);
mean_velocity = mean(magnitude_velocityMove);

% Calculate movement duration in milliseconds
movement_duration = numel(onset_index:offset_index) / fsamp * 1000;

% Calculate smoothness using integrated squared jerk
jerk = diff(magnitude_velocityMove, 1); % Jerk is the derivative of acceleration
smoothness_index = sum(jerk.^2) / length(magnitude_velocityMove); % Integrated squared jerk

% Detect hold phase after movement
indexes = offset_index + find(magnitude_velocity_z(offset_index:end) < holdThreshold);
if ~isempty(indexes)
    indexes = SplitVec(indexes, 'consecutive');
    bridgeInds = bridgeCellSplits(indexes, fsamp / 2);
    index_length = cellfun(@length, bridgeInds);
    [~, indmax] = max(index_length);
    holdIndex = bridgeInds{indmax};
    magnitude_velocityHold = magnitude_velocity(holdIndex - 1); % Index minus 1 as it is a derivative
    hold_stability = std(magnitude_velocityHold) / sqrt(numel(magnitude_velocityHold));
else
    % No hold phase detected
    holdIndex = nan;
    magnitude_velocityHold = nan;
    hold_stability = nan;
end

% Warn if the movement or hold index is very short
if numel(holdIndex) < 256
    warning('Hold Index is very short!');
end

if numel(moveIndex) < 256
    warning('Move Index is very short!');
end

end
