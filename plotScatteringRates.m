function plotScatteringRates(scattering, params)
% PLOTSCATTERINGRATES - Visualize scattering rates and cumulative table
%
% Inputs:
%   scattering - Structure from calculateScatteringRates
%   params - Structure from initializePhysicalConstants
%
% Outputs: None (generates figures)
%
% Dependencies: None

%% Figure 1: Individual scattering rates
figure(1)
hold on
grid on
semilogy(params.e / params.q, scattering.acoustic_scatter, ...
    'DisplayName', 'acoustic', 'LineWidth', 3)
semilogy(params.e / params.q, scattering.rate_intervally_g_absorb_nonpara, ...
    'DisplayName', 'g absorption', 'LineStyle', '--', 'LineWidth', 3)
semilogy(params.e / params.q, scattering.rate_intervally_g_emission_nonpara, ...
    'DisplayName', 'g emission', 'LineStyle', '-.', 'LineWidth', 3)
semilogy(params.e / params.q, scattering.rate_intervally_f_absorb_nonpara, ...
    'DisplayName', 'f absorption', 'LineStyle', '--', 'LineWidth', 3)
semilogy(params.e / params.q, scattering.rate_intervally_f_emission_nonpara, ...
    'DisplayName', 'f emission', 'LineStyle', '-.', 'LineWidth', 3)
hold off
xlabel('Energy (eV)')
ylabel('Scattering rate (s^{-1})')
title('Scattering Rates')
set(gca, 'yscale', 'log')
legend show

%% Figure 2: Normalized cumulative scattering rates
max_cumu = scattering.max_cumu_sacttering;

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

figure(2)
hold on
grid on
plot(params.e / params.q, normalize_acoustic, ...
    'DisplayName', 'acoustic', 'LineWidth', 3)
plot(params.e / params.q, normalize_inter_abs_g, ...
    'DisplayName', 'g absorption', 'LineStyle', '--', 'LineWidth', 3)
plot(params.e / params.q, normalize_inter_emm_g, ...
    'DisplayName', 'g emission', 'LineStyle', '-.', 'LineWidth', 3)
plot(params.e / params.q, normalize_inter_abs_f, ...
    'DisplayName', 'f absorption', 'LineStyle', '--', 'LineWidth', 3)
plot(params.e / params.q, normalize_inter_emm_f, ...
    'DisplayName', 'f emission', 'LineStyle', '-.', 'LineWidth', 3)
xlabel('Energy (eV)')
ylabel('Normalized Cumulative Scattering Rate')
title('Normalized Scattering Rates')
hold off
set(gca, 'yscale', 'linear')
legend show

end