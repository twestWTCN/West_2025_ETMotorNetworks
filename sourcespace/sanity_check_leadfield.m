function sanity_check_leadfield(elec_realigned, ...
    leadfield, headmodel, mesh_eeg)

%   This function displays leadfields along with the headmodel and
%   the electrodes to figure out problems with the estimation; adapted from:
%   https://www.fieldtriptoolbox.org/workshop/oslo2019/forward_modeling/#head-models-component-1

%   ## Version 1.0

%   Copyright (C) November 2021
%   D. Pedrosa
%   University Hospital of Gießen and Marburg
%
%   This software may be used, copied, or redistributed as long as it is
%   not sold and this copyright notice is reproduced on each copy made.
%   This routine is provided as is without any express or implied
%   warranties whatsoever.


figure('units', 'normalized', 'outerposition', [0 0 0.5 0.5])
source_index = 1000; %% random source
sensory_dipole_current = 100e-9; % Am (realistic)

n_sensors = length(elec_realigned.label);

inside_sources = find(leadfield.inside);
inside_index = inside_sources(source_index);
lead = leadfield.leadfield{inside_index};
xs = zeros(1, n_sensors);
ys = zeros(1, n_sensors);
zs = zeros(1, n_sensors);
voltages = zeros(1, n_sensors);
titles = {'Lead field (x)' 'Lead field (y)' 'Lead field (z)'};

% get the xyz and norm

for sensor_index = 1:n_sensors
    this_x = lead(sensor_index, 1);
    this_y = lead(sensor_index, 2);
    this_z = lead(sensor_index, 3);
    this_norm = norm(lead(sensor_index, :));
    xs(sensor_index) = this_x * sensory_dipole_current;
    ys(sensor_index) = this_y * sensory_dipole_current;
    zs(sensor_index) = this_z * sensory_dipole_current;
    voltages(sensor_index) = this_norm * sensory_dipole_current;
end

% plot xyz
axes = {xs ys zs};

for axis_index = 1:3
    this_axis = axes{axis_index};
    subplot(1, 3, axis_index)
    hold on
    ft_plot_topo3d(elec_realigned.chanpos, this_axis, 'facealpha', 0.8)
%     if strcmp(headmodel.type, 'dipoli')
%         caxis([-10e-6, 10e-6])
%     elseif strcmp(headmodel.type, 'openmeeg')
        caxis([0, max(voltages)])
%     end
    c = colorbar('location', 'southoutside');
    c.Label.String = 'Lead field (V)';
    axis tight
    ft_plot_mesh(mesh_eeg, 'facealpha', 0.10);
    ft_plot_sens(elec_realigned, 'elecsize', 20);
    title(titles{axis_index})
    plot3(leadfield.pos(inside_index, 1), ...
        leadfield.pos(inside_index, 2), ...
        leadfield.pos(inside_index, 3), 'bo', ...
        'markersize', 20, 'markerfacecolor', 'r')
end

% plot norm
figure('units', 'normalized', 'outerposition', [0 0 0.5 0.85])
hold on
ft_plot_topo3d(elec_realigned.chanpos, voltages, 'facealpha', 0.8)
% if strcmp(headmodel.type, 'dipoli')
%     caxis([0, 10e-6])
% elseif strcmp(headmodel.type, 'openmeeg')
    caxis([0, max(voltages)])
% end
c = colorbar('location', 'eastoutside');
c.Label.String = 'Lead field (V)';
axis tight
ft_plot_mesh(mesh_eeg, 'facealpha', 0.10);
ft_plot_sens(elec_realigned, 'elecsize', 20);
title('Leadfield magnitude')
plot3(leadfield.pos(inside_index, 1), ...
    leadfield.pos(inside_index, 2), ...
    leadfield.pos(inside_index, 3), 'bo', ...
    'markersize', 20, 'markerfacecolor', 'r')
view(-90, 0)
end

