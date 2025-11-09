% SMOKE_TEST - Quick validation that refactored code runs end-to-end
%
% This script runs a minimal version of the simulation to verify
% that all modules are properly connected and produce results.

clear;
close all;
clc;

fprintf('=== Smoke Test: Refactored Monte Carlo Simulation ===\n\n');

%% Modify parameters for quick test
fprintf('Step 1: Initialize constants...\n');
params = initializePhysicalConstants();

% Override for faster test
params.number_electrons = 10;       % Reduced from 100
params.time_simulation = 1e-12;     % Reduced from 10e-12
params.number_fieldy = 3;           % Reduced from 10
params.fieldy = linspace(params.fieldy_start, params.fieldy_end, params.number_fieldy);
params.time_observation = 0:params.time_step:params.time_simulation;

fprintf('  Using %d electrons, %d fields, %.2e s simulation time\n', ...
    params.number_electrons, params.number_fieldy, params.time_simulation);

%% Run pipeline
fprintf('Step 2: Calculate scattering rates...\n');
scattering = calculateScatteringRates(params);

fprintf('Step 3: Build scattering table...\n');
scatter_table = buildScatteringTable(scattering, params);

fprintf('Step 4: Test ensemble initialization...\n');
ensemble = initializeElectronEnsemble(params, scattering);
fprintf('  Initial energy range: %.3f - %.3f eV\n', ...
    min(ensemble.energy_electrons)/params.q, max(ensemble.energy_electrons)/params.q);

fprintf('\nStep 5: Run Monte Carlo (may take a moment)...\n');
try
    sim_results = runMonteCarloSimulation(params, scatter_table, scattering);
    fprintf('  ✓ Simulation completed successfully\n');
catch ME
    if contains(ME.message, 'drift.m')
        fprintf('  ⚠ Simulation failed: drift.m not found (expected)\n');
        fprintf('    drift.m is an external dependency that must be provided\n');
        return;
    else
        rethrow(ME);
    end
end

fprintf('Step 6: Average results...\n');
avg_results = averageResults(sim_results, params);

fprintf('Step 7: Generate plots...\n');
plotScatteringRates(scattering, params);
plotSimulationResults(avg_results, params, sim_results);

fprintf('\n=== Smoke Test Complete ===\n');
fprintf('All modules executed successfully.\n');
fprintf('Final drift velocity at max field: %.2e m/s\n', ...
    -avg_results.vy_avg_time(end));