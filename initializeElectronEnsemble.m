function ensemble = initializeElectronEnsemble(params, scattering)
% INITIALIZEELECTRONENSEMBLE - Set up initial electron states
%
% Inputs:
%   params - Structure from initializePhysicalConstants
%   scattering - Structure from calculateScatteringRates
%
% Outputs:
%   ensemble - Structure containing initial electron states
%     .energy_electrons - Initial energies (thermal distribution)
%     .initial_valley - Valley assignments (1, 2, or 3)
%     .kx, .ky, .kz - Initial wave vectors
%     .time_flight - Initial time to next scattering event
%
% Dependencies: None

n_electrons = params.number_electrons;

% Initialize energies from thermal distribution
ensemble.energy_electrons = -1 * (1.5 * params.V_t) * log(rand(1, n_electrons));

% Random valley assignment (equiprobable)
random_array = 3.0 * rand(1, n_electrons);
ensemble.initial_valley = zeros(1, n_electrons);
ensemble.initial_valley(random_array <= 1) = 1;
ensemble.initial_valley(random_array <= 2 & random_array > 1) = 2;
ensemble.initial_valley(random_array <= 3 & random_array > 2) = 3;

% Wave vector calculation
k_magnitude = params.smh * sqrt(ensemble.energy_electrons .* ...
    (1 + params.alpha_nonpara * ensemble.energy_electrons));

phi = 2 * pi * rand(1, n_electrons);
ct = 1 - 2 * rand(1, n_electrons);  % cos(theta)
st = sqrt(1 - ct.^2);               % sin(theta)

ensemble.kx = k_magnitude .* st .* cos(phi);
ensemble.ky = k_magnitude .* st .* sin(phi);
ensemble.kz = k_magnitude .* ct;

% Time of flight calculation
random_array = rand(1, n_electrons);
while any(random_array <= 1e-5)
    random_array(random_array <= 1e-5) = rand(1, sum(random_array <= 1e-5));
end

tau_max = 1 / scattering.max_cumu_sacttering;
ensemble.time_flight = -1 * log(random_array) * tau_max;

end