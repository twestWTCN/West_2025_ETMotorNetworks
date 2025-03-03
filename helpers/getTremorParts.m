function tremrng = getTremorParts(part)  
% times of at which tremor is computed
if part == 1
    tremrng = [-2 0]; % prior to posture (rest tremor?)
elseif part == 2
    tremrng = [1 3]; % after posture (postural tremor)
elseif part == 3
    tremrng = [0 2]; % after cue
elseif part == 4
    tremrng = [0 2]; % kinetic tremor
elseif part == 5
    tremrng = [0.5 2.5]; % tremor at hold
end