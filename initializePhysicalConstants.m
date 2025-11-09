function params = initializePhysicalConstants()
% INITIALIZEPHYSICALCONSTANTS - Set up all physical and material parameters
%
% Outputs:
%   params - Structure containing all simulation parameters
%     .energy_levels - Number of energy discretization points
%     .density_crystal - Crystal density (kg/m^3)
%     .doping - Doping concentration (m^-3)
%     .sound_vel - Sound velocity (m/s)
%     .h, .h_cross - Planck constants
%     .m_not - Electron rest mass
%     .mass_longitudinal, .mass_transverse - Effective mass ratios
%     .effective_mass_DOS, .effective_mass_CB - Computed effective masses
%     .tm1, tm2, tm3 - Mass transformation factors
%     .hm1, hm2, hm3 - Velocity prefactors
%     .q, .kb - Elementary charge, Boltzmann constant
%     .sigma - Deformation potential
%     .temperature - Lattice temperature (K)
%     .energy_max - Maximum energy for discretization
%     .e - Energy grid vector
%     .de - Energy step
%     .alpha_nonpara - Non-parabolicity parameter
%     .V_t - Thermal voltage
%     .smh - k-magnitude prefactor
%     .w0_g, w0_f - Phonon energies for g/f valleys
%     .delta_fi_zero_g, delta_fi_zero_f - Valley energy offsets
%     .number_final_valley_g_absorb, number_final_valley_f_absorb
%     .coupling_const_g_absorb, coupling_const_f_absorb
%     .cnst_absorb_g, cnst_absorb_f - Intervalley coupling constants
%     .rnq_g_absorb, rnq_f_absorb - Phonon occupation numbers
%
% Dependencies: None

% Simulation parameters
params.energy_levels = 5000;

% Material properties (Silicon)
params.density_crystal = 2329;  % kg/m^3
params.doping = 1.0e21;         % m^-3
params.sound_vel = 9040;        % m/s

% Physical constants
params.h = 6.625e-34;           % Planck constant (J*s)
params.h_cross = params.h / (2 * pi);
params.m_not = 9.11e-31;        % Electron rest mass (kg)

% Effective mass ratios
params.mass_longitudinal = 0.96;
params.mass_transverse = 0.196;

% Effective mass calculations
params.effective_mass_DOS = (params.mass_longitudinal * ...
    params.mass_transverse * params.mass_transverse)^(1/3) * params.m_not;
params.effective_mass_CB = 3.0 / ((1.0 / params.mass_longitudinal) + ...
    (2.0 / params.mass_transverse)) * params.m_not;

% Mass transformation factors
params.tm1 = sqrt(params.effective_mass_CB / params.mass_longitudinal / params.m_not);
params.tm2 = sqrt(params.effective_mass_CB / params.mass_transverse / params.m_not);
params.tm3 = sqrt(params.effective_mass_CB / params.mass_transverse / params.m_not);

% Velocity prefactors
params.hm1 = (params.h_cross * params.tm1) / params.effective_mass_CB;
params.hm2 = (params.h_cross * params.tm2) / params.effective_mass_CB;
params.hm3 = (params.h_cross * params.tm3) / params.effective_mass_CB;

% Fundamental constants
params.q = 1.602e-19;           % Elementary charge (C)
params.kb = 1.38064e-23;        % Boltzmann constant (J/K)
params.sigma = 9 * params.q;    % Deformation potential (J)
params.temperature = 300;       % Temperature (K)
params.V_t = params.kb * params.temperature;  % Thermal voltage (J)

% Energy grid
params.energy_max = 5 * params.q;
params.e = linspace(0, params.energy_max, params.energy_levels);
params.de = params.e(2) - params.e(1);

% Non-parabolicity parameter
params.alpha_nonpara = 1.5 / params.q;

% k-magnitude prefactor
params.smh = sqrt(2.0 * params.effective_mass_CB) / params.h_cross;

% Intervalley scattering parameters - g-type
params.number_final_valley_g_absorb = 1;
params.coupling_const_g_absorb = 5.23e10 * params.q;
params.w0_g = 0.063 * params.q;
params.delta_fi_zero_g = 0.0;
params.cnst_absorb_g = (params.number_final_valley_g_absorb * ...
    params.coupling_const_g_absorb^2 * params.effective_mass_DOS^1.5) / ...
    (sqrt(2) * pi * params.density_crystal * params.w0_g * params.h_cross^2);
params.rnq_g_absorb = 1.0 / (exp(params.w0_g / params.V_t) - 1.0);

% Intervalley scattering parameters - f-type
params.number_final_valley_f_absorb = 4;
params.coupling_const_f_absorb = 5.23e10 * params.q;
params.w0_f = 0.059 * params.q;
params.delta_fi_zero_f = 0.0 * params.q;
params.cnst_absorb_f = (params.number_final_valley_f_absorb * ...
    params.coupling_const_f_absorb^2 * params.effective_mass_DOS^1.5) / ...
    (sqrt(2) * pi * params.density_crystal * params.w0_f * params.h_cross^2);
params.rnq_f_absorb = 1.0 / (exp(params.w0_f / params.V_t) - 1.0);

% Monte Carlo simulation parameters
params.number_electrons = 100;
params.time_simulation = 10e-12;  % seconds
params.time_step = 1e-14;         % seconds
params.time_observation = 0:params.time_step:params.time_simulation;
params.time_index = int32(params.time_simulation / params.time_step);
params.fieldy_start = 1e6;        % V/m
params.fieldy_end = 8e6;          % V/m
params.number_fieldy = 10;
params.fieldy = linspace(params.fieldy_start, params.fieldy_end, params.number_fieldy);

% Derived parameters
params.hhm = params.h_cross^2 / (params.effective_mass_CB * 2);

end