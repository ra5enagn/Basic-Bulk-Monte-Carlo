% MAIN - Monte Carlo electron transport simulation orchestrator
%
% This script coordinates the full workflow:
%   1. Initialize physical constants and material parameters
%   2. Calculate scattering rates (acoustic, intervalley)
%   3. Build cumulative scattering table
%   4. Run Monte Carlo simulation over multiple electric fields
%   5. Average results and visualize
%
% Dependencies:
%   - initializePhysicalConstants.m
%   - calculateScatteringRates.m
%   - buildScatteringTable.m
%   - initializeElectronEnsemble.m
%   - runMonteCarloSimulation.m
%   - averageResults.m
%   - plotScatteringRates.m
%   - plotSimulationResults.m
%   - drift.m (external, must be in path)

clear;
close all;
clc;

%% Step 1: Initialize Physical Constants
fprintf('Initializing physical constants...\n');
params = initializePhysicalConstants();

%% Step 2: Calculate Scattering Rates
fprintf('Calculating scattering rates...\n');
scattering = calculateScatteringRates(params);

%% Step 3: Build Scattering Table
fprintf('Building cumulative scattering table...\n');
scatter_table = buildScatteringTable(scattering, params);

%% Step 4: Plot Scattering Rates
fprintf('Plotting scattering rates...\n');
plotScatteringRates(scattering, params);

%% Step 5: Run Monte Carlo Simulation
fprintf('Running Monte Carlo simulation...\n');
sim_results = runMonteCarloSimulation(params, scatter_table, scattering);

%% Step 6: Average Results
fprintf('Averaging results over time and field...\n');
avg_results = averageResults(sim_results, params);

%% Step 7: Plot Simulation Results
fprintf('Plotting simulation results...\n');
plotSimulationResults(avg_results, params, sim_results);

fprintf('Simulation complete.\n');