function plotSimulationResults(avg_results, params, sim_results)
% PLOTSIMULATIONRESULTS - Visualize velocity-field curves and time evolution
%
% Inputs:
%   avg_results - Structure from averageResults
%   params - Structure from initializePhysicalConstants
%   sim_results - Structure from runMonteCarloSimulation
%
% Outputs: None (generates figures)
%
% Dependencies: None

%% Figure 3: Velocity-field characteristics
figure(3)

subplot(2, 2, 1)
grid on
plot(params.fieldy, avg_results.e_avg_time / params.q)
title('Average Energy vs Electric Field')
xlabel('Electric Field (V/m)')
ylabel('Energy (eV)')
grid off

subplot(2, 2, 2)
grid on
plot(params.fieldy, avg_results.vx_avg_time)
title('v_x vs Electric Field')
xlabel('Electric Field (V/m)')
ylabel('Velocity (m/s)')
grid off

subplot(2, 2, 3)
grid on
plot(params.fieldy, -1 * avg_results.vy_avg_time)
title('v_y vs Electric Field (drift direction)')
xlabel('Electric Field (V/m)')
ylabel('Drift Velocity (m/s)')
grid off

subplot(2, 2, 4)
grid on
plot(params.fieldy, avg_results.vz_avg_time)
title('v_z vs Electric Field')
xlabel('Electric Field (V/m)')
ylabel('Velocity (m/s)')
grid off

sgtitle('Velocity-Field Characteristics')

%% Figure 4: Time evolution at maximum field
figure(4)

fy_i = params.number_fieldy;  % Last field index
time_ps = params.time_observation / 1e-12;

subplot(2, 2, 1)
grid on
plot(time_ps, sim_results.e_avg(fy_i, :) / params.q)
title('Energy Evolution at Maximum Field')
xlabel('Time (ps)')
ylabel('Energy (eV)')
grid off

subplot(2, 2, 2)
grid on
plot(time_ps, sim_results.vx_avg(fy_i, :))
title('v_x Evolution at Maximum Field')
xlabel('Time (ps)')
ylabel('Velocity (m/s)')
grid off

subplot(2, 2, 3)
grid on
plot(time_ps, -1 * sim_results.vy_avg(fy_i, :))
title('v_y Evolution at Maximum Field')
xlabel('Time (ps)')
ylabel('Drift Velocity (m/s)')
grid off

subplot(2, 2, 4)
grid on
plot(time_ps, sim_results.vz_avg(fy_i, :))
title('v_z Evolution at Maximum Field')
xlabel('Time (ps)')
ylabel('Velocity (m/s)')
grid off

sgtitle(sprintf('Time Evolution at Field = %.2e V/m', params.fieldy(end)))

end