function scattering = calculateScatteringRates(params)
% CALCULATESCATTERINGRATES - Compute all scattering rate mechanisms
%
% Inputs:
%   params - Structure from initializePhysicalConstants
%
% Outputs:
%   scattering - Structure containing scattering rate arrays
%     .acoustic_scatter - Acoustic phonon scattering rate
%     .rate_intervally_g_absorb_nonpara - g-valley absorption
%     .rate_intervally_g_emission_nonpara - g-valley emission
%     .rate_intervally_f_absorb_nonpara - f-valley absorption
%     .rate_intervally_f_emission_nonpara - f-valley emission
%     .total_scattering - Sum of all scattering rates
%     .max_cumu_sacttering - Maximum scattering rate (for normalization)
%
% Dependencies: None

%% Acoustic Scattering
Cnst_nume = params.effective_mass_DOS^1.5 * params.sigma^2 * params.V_t * sqrt(2);
Cnst_denom = pi * params.h_cross^4 * params.density_crystal * params.sound_vel^2;
Cnst = Cnst_nume / Cnst_denom;

final_e_nonpara = params.e .* (1 + params.alpha_nonpara * params.e);
scattering.acoustic_scatter = (Cnst * sqrt(final_e_nonpara)) .* ...
    (1 + 2 * params.alpha_nonpara * params.e);

%% Intervalley Scattering - g-type
% Absorption
energy_final_state_inter_absorb_g = params.e + params.w0_g - params.delta_fi_zero_g;
intervally_g_absorb_nonpara = energy_final_state_inter_absorb_g .* ...
    (1 + params.alpha_nonpara * energy_final_state_inter_absorb_g);

scattering.rate_intervally_g_absorb_nonpara = zeros(size(params.e));
temperory_absorb = zeros(size(params.e));

valid_idx = energy_final_state_inter_absorb_g > 0;
temperory_absorb(valid_idx) = params.cnst_absorb_g * params.rnq_g_absorb;
scattering.rate_intervally_g_absorb_nonpara(valid_idx) = ...
    sqrt(intervally_g_absorb_nonpara(valid_idx)) .* ...
    (1 + 2 * params.alpha_nonpara * energy_final_state_inter_absorb_g(valid_idx)) .* ...
    temperory_absorb(valid_idx);

% Emission
const_emission = (1 + params.rnq_g_absorb) * params.cnst_absorb_g;
energy_final_state_inter_emission_g = params.e - params.w0_g - params.delta_fi_zero_g;
intervally_g_emission_nonpara = energy_final_state_inter_emission_g .* ...
    (1 + params.alpha_nonpara * energy_final_state_inter_emission_g);

scattering.rate_intervally_g_emission_nonpara = zeros(size(params.e));
temperory_emission_g = zeros(size(params.e));

valid_idx = energy_final_state_inter_emission_g > 0;
temperory_emission_g(valid_idx) = ...
    sqrt(intervally_g_emission_nonpara(valid_idx)) .* ...
    (1 + 2 * params.alpha_nonpara * energy_final_state_inter_emission_g(valid_idx));
scattering.rate_intervally_g_emission_nonpara(valid_idx) = ...
    const_emission * temperory_emission_g(valid_idx);

%% Intervalley Scattering - f-type
% Absorption
energy_final_state_inter_absorb_f = params.e + params.w0_f - params.delta_fi_zero_f;
intervally_f_absorb_nonpara = energy_final_state_inter_absorb_f .* ...
    (1 + params.alpha_nonpara * energy_final_state_inter_absorb_f);

scattering.rate_intervally_f_absorb_nonpara = zeros(size(params.e));
temperory_absorb = zeros(size(params.e));

valid_idx = energy_final_state_inter_absorb_f > 0;
temperory_absorb(valid_idx) = params.cnst_absorb_f * params.rnq_f_absorb;
scattering.rate_intervally_f_absorb_nonpara(valid_idx) = ...
    sqrt(intervally_f_absorb_nonpara(valid_idx)) .* ...
    (1 + 2 * params.alpha_nonpara * energy_final_state_inter_absorb_f(valid_idx)) .* ...
    temperory_absorb(valid_idx);

% Emission
const_emission = (1 + params.rnq_f_absorb) * params.cnst_absorb_f;
energy_final_state_inter_emission_f = params.e - params.w0_f - params.delta_fi_zero_f;
intervally_f_emission_nonpara = energy_final_state_inter_emission_f .* ...
    (1 + params.alpha_nonpara * energy_final_state_inter_emission_f);

scattering.rate_intervally_f_emission_nonpara = zeros(size(params.e));
temperory_emission_f = zeros(size(params.e));

valid_idx = energy_final_state_inter_emission_f > 0;
temperory_emission_f(valid_idx) = ...
    sqrt(intervally_f_emission_nonpara(valid_idx)) .* ...
    (1 + 2 * params.alpha_nonpara * energy_final_state_inter_emission_f(valid_idx));
scattering.rate_intervally_f_emission_nonpara(valid_idx) = ...
    const_emission * temperory_emission_f(valid_idx);

%% Total Scattering
scattering.total_scattering = scattering.acoustic_scatter + ...
    scattering.rate_intervally_f_emission_nonpara + ...
    scattering.rate_intervally_f_absorb_nonpara + ...
    scattering.rate_intervally_g_emission_nonpara + ...
    scattering.rate_intervally_g_absorb_nonpara;

scattering.max_cumu_sacttering = max(scattering.total_scattering);

end