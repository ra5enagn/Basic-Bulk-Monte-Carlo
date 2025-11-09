function scatter_table = buildScatteringTable(scattering, params)
% BUILDSCATTERINGTABLE - Create cumulative normalized scattering table
%
% Inputs:
%   scattering - Structure from calculateScatteringRates
%   params - Structure from initializePhysicalConstants
%
% Outputs:
%   scatter_table - 6 x energy_levels matrix of cumulative probabilities
%     Row 1: Acoustic scattering
%     Row 2: + f-valley emission
%     Row 3: + f-valley absorption
%     Row 4: + g-valley emission
%     Row 5: + g-valley absorption
%     Row 6: Total (=1)
%
% Dependencies: None

max_cumu = scattering.max_cumu_sacttering;

% Cumulative normalized scattering rates
normalize_acoustic = scattering.acoustic_scatter / max_cumu;

normalize_inter_emm_f = (scattering.acoustic_scatter + ...
    scattering.rate_intervally_f_emission_nonpara) / max_cumu;

normalize_inter_abs_f = (scattering.acoustic_scatter + ...
    scattering.rate_intervally_f_emission_nonpara + ...
    scattering.rate_intervally_f_absorb_nonpara) / max_cumu;

normalize_inter_emm_g = (scattering.acoustic_scatter + ...
    scattering.rate_intervally_f_emission_nonpara + ...
    scattering.rate_intervally_f_absorb_nonpara + ...
    scattering.rate_intervally_g_emission_nonpara) / max_cumu;

normalize_inter_abs_g = (scattering.acoustic_scatter + ...
    scattering.rate_intervally_f_emission_nonpara + ...
    scattering.rate_intervally_f_absorb_nonpara + ...
    scattering.rate_intervally_g_emission_nonpara + ...
    scattering.rate_intervally_g_absorb_nonpara) / max_cumu;

% Build table
scatter_table = ones(6, params.energy_levels);
scatter_table(1, :) = normalize_acoustic;
scatter_table(2, :) = normalize_inter_emm_f;
scatter_table(3, :) = normalize_inter_abs_f;
scatter_table(4, :) = normalize_inter_emm_g;
scatter_table(5, :) = normalize_inter_abs_g;

end