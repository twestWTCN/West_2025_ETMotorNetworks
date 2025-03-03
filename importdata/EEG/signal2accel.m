function Acc = signal2accel(X)

% Define the constants
amplifier_resolution = 1; % microvolts per unit
conversion_factor = 300; % mV/g
g_to_m_s2 = 9.81; % m/s^2 per g

% Step 1: Convert amplifier units to microvolts
signal_microvolts = X.* amplifier_resolution;

% Step 2: Convert microvolts to millivolts
signal_millivolts = signal_microvolts./1000;

% Step 3: Convert millivolts to acceleration in g
acceleration_g = signal_millivolts./conversion_factor;

% Step 4: Convert acceleration from g to m/s^2
acceleration_m_s2 = acceleration_g.*g_to_m_s2;
Acc = acceleration_m_s2;
