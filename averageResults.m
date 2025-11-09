function avg_results = averageResults(sim_results, params)
% AVERAGERESULTS - Compute time-averaged quantities for each field
%
% Inputs:
%   sim_results - Structure from runMonteCarloSimulation
%   params - Structure from initializePhysicalConstants
%
% Outputs:
%   avg_results - Structure containing time-averaged results
%     .e_avg_time - Average energy for each field (scalar per field)
%     .vx_avg_time, .vy_avg_time, .vz_avg_time - Average velocities
%
% Dependencies: None

n_fields = params.number_fieldy;
n_time = size(sim_results.e_avg, 2);

% Calculate time averages over second half of simulation (steady state)
start_idx = int32(n_time / 2);
end_idx = n_time - 1;

avg_results.e_avg_time = zeros(n_fields, 1);
avg_results.vx_avg_time = zeros(n_fields, 1);
avg_results.vy_avg_time = zeros(n_fields, 1);
avg_results.vz_avg_time = zeros(n_fields, 1);

for fy_i = 1:n_fields
    avg_results.e_avg_time(fy_i) = mean(sim_results.e_avg(fy_i, start_idx:end_idx));
    avg_results.vx_avg_time(fy_i) = mean(sim_results.vx_avg(fy_i, start_idx:end_idx));
    avg_results.vy_avg_time(fy_i) = mean(sim_results.vy_avg(fy_i, start_idx:end_idx));
    avg_results.vz_avg_time(fy_i) = mean(sim_results.vz_avg(fy_i, start_idx:end_idx));
end

end