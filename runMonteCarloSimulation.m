function sim_results = runMonteCarloSimulation(params, scatter_table, scattering)
% RUNMONTECARLOSIMULATION - Execute Monte Carlo transport simulation
%
% Inputs:
%   params - Structure from initializePhysicalConstants
%   scatter_table - Cumulative scattering table from buildScatteringTable
%   scattering - Structure from calculateScatteringRates
%
% Outputs:
%   sim_results - Structure containing simulation results
%     .e_avg - Average energy vs time for each field (fieldy x time)
%     .kx_avg, .ky_avg, .kz_avg - Average k-vectors
%     .vx_avg, .vy_avg, .vz_avg - Average velocities
%     .e_history_final - Energy history for last field (time x electrons)
%     .vx_history_final, etc. - Velocity histories for last field
%
% Dependencies:
%   - drift.m (external function, must be in path)

n_electrons = params.number_electrons;
n_fields = params.number_fieldy;
n_time = length(params.time_observation);

% Preallocate result arrays
sim_results.e_avg = zeros(n_fields, n_time);
sim_results.kx_avg = zeros(n_fields, n_time);
sim_results.ky_avg = zeros(n_fields, n_time);
sim_results.kz_avg = zeros(n_fields, n_time);
sim_results.vx_avg = zeros(n_fields, n_time);
sim_results.vy_avg = zeros(n_fields, n_time);
sim_results.vz_avg = zeros(n_fields, n_time);

delete(gcp('nocreate')); % Close existing pool
parpool(4); % Start a new pool with 4 workers
feature('numcores')

% Loop over electric fields
for fy_i = 1:n_fields
    fprintf('  Field iteration number: %d\n', fy_i);
    fy = params.fieldy(fy_i);
    fprintf('  Field Value (y-axis): %.2e V/m\n', fy);
    
    % Initialize electron ensemble
    ensemble = initializeElectronEnsemble(params, scattering);
    
    % Extract initial states
    energy_electrons = ensemble.energy_electrons;
    initial_valley = ensemble.initial_valley;
    kx = ensemble.kx;
    ky = ensemble.ky;
    kz = ensemble.kz;
    dte = ensemble.time_flight;
    
    % Preallocate history matrices
    e_history = zeros(n_time, n_electrons);
    kx_history = zeros(n_time, n_electrons);
    ky_history = zeros(n_time, n_electrons);
    kz_history = zeros(n_time, n_electrons);
    vx_history = zeros(n_time, n_electrons);
    vy_history = zeros(n_time, n_electrons);
    vz_history = zeros(n_time, n_electrons);
    
    vx = zeros(1, n_electrons);
    vy = zeros(1, n_electrons);
    vz = zeros(1, n_electrons);
    
    % Electric field components
    fx = 0.0;
    fz = 0.0;
    
    % Time evolution loop
    i = 1;
    while i <= n_time
        parfor j = 1:n_electrons
            % Determine drift time
            if dte(j) >= params.time_step
                dt2 = params.time_step;
            else
                dt2 = dte(j);
            end
            
            % Drift step
            [nkx, nky, nkz, energy] = drift(dt2, initial_valley(j), ...
                kx(j), ky(j), kz(j), fx, fy, fz, ...
                params.tm1, params.tm2, params.tm3, params.hhm, ...
                energy_electrons(j), params.alpha_nonpara);
            
            kx(j) = nkx;
            ky(j) = nky;
            kz(j) = nkz;
            energy_electrons(j) = energy;
            
            % Scattering loop
            while dte(j) < params.time_step
                dte2 = dte(j);
                
                % Determine scattering type
                loc = floor(energy_electrons(j) / params.de) + 1;
                if loc == 0
                    loc = 1;
                elseif loc > params.energy_levels
                    loc = params.energy_levels;
                end
                
                % Select scattering mechanism
                random_sc_number = rand();
                m1 = 6;
                for m = 1:size(scatter_table, 1)
                    if scatter_table(m, loc) > random_sc_number
                        m1 = m;
                        break;
                    end
                end
                
                % Apply scattering
                if m1 <= 5
                [new_energy, new_kx, new_ky, new_kz, new_valley] = ...
                    applyScattering(m1, energy_electrons(j), params, initial_valley(j));
                
                % Update state
                energy_electrons(j) = new_energy;
                kx(j) = new_kx;
                ky(j) = new_ky;
                kz(j) = new_kz;
                initial_valley(j) = new_valley;
                end
                
                % Generate new time to next scattering
                rr = rand();
                while rr <= 1e-6
                    rr = rand();
                end
                dt3 = -1 * log(rr) / scattering.max_cumu_sacttering;
                dtp = params.time_step - dte2;
                
                if dt3 <= dtp
                    dt2 = dt3;
                else
                    dt2 = dtp;
                end
                
                % Drift again
                [nkx, nky, nkz, energy] = drift(dt2, initial_valley(j), ...
                    kx(j), ky(j), kz(j), fx, fy, fz, ...
                    params.tm1, params.tm2, params.tm3, params.hhm, ...
                    energy_electrons(j), params.alpha_nonpara);
                
                kx(j) = nkx;
                ky(j) = nky;
                kz(j) = nkz;
                energy_electrons(j) = energy;
                dte2 = dte2 + dt3;
                dte(j) = dte2;
            end
            
            dte(j) = dte(j) - params.time_step;
            
            % Calculate velocities
            [vx(j), vy(j), vz(j)] = calculateVelocity(initial_valley(j), ...
                kx(j), ky(j), kz(j), energy_electrons(j), params);
        end
        
        % Store history
        e_history(i, :) = energy_electrons;
        kx_history(i, :) = kx;
        ky_history(i, :) = ky;
        kz_history(i, :) = kz;
        vx_history(i, :) = vx;
        vy_history(i, :) = vy;
        vz_history(i, :) = vz;
        
        i = i + 1;
    end
    
    % Average over ensemble
    sim_results.e_avg(fy_i, :) = mean(e_history, 2);
    sim_results.kx_avg(fy_i, :) = mean(kx_history, 2);
    sim_results.ky_avg(fy_i, :) = mean(ky_history, 2);
    sim_results.kz_avg(fy_i, :) = mean(kz_history, 2);
    sim_results.vx_avg(fy_i, :) = mean(vx_history, 2);
    sim_results.vy_avg(fy_i, :) = mean(vy_history, 2);
    sim_results.vz_avg(fy_i, :) = mean(vz_history, 2);
end

% Store final field history for detailed plotting
sim_results.e_history_final = e_history;
sim_results.vx_history_final = vx_history;
sim_results.vy_history_final = vy_history;
sim_results.vz_history_final = vz_history;

end

%% Local helper function: applyScattering
function [new_energy, new_kx, new_ky, new_kz, new_valley] = ...
    applyScattering(scatter_type, energy, params, valley)
% APPLYSCATTERING - Apply scattering mechanism to electron state
% EXACTLY matches monolithic code logic
%
% Inputs:
%   scatter_type - Integer 1-5 indicating scattering mechanism
%   energy - Current electron energy
%   params - Parameter structure
%   valley - Current valley (1, 2, or 3)
%
% Outputs:
%   new_energy - Energy after scattering
%   new_kx, new_ky, new_kz - k-vector components after scattering
%   new_valley - Valley after scattering (may change for intervalley)

new_valley = valley;
new_energy = energy;

switch scatter_type
    case 1  % Acoustic scattering
        new_energy = energy;
        rknew = params.smh * sqrt(new_energy * (1 + params.alpha_nonpara * new_energy));
        phi = 2 * pi * rand();
        ct = 1 - 2 * rand();
        st = sqrt(1 - ct^2);
        new_kx = rknew * st * cos(phi);
        new_ky = rknew * st * sin(phi);
        new_kz = rknew * ct;
        
    case 2  % f-valley emission
        new_energy = energy - params.w0_f - params.delta_fi_zero_f;
        rknew = params.smh * sqrt(new_energy * (1 + params.alpha_nonpara * new_energy));
        phi = 2 * pi * rand();
        ct = 1 - 2 * rand();
        st = sqrt(1 - ct^2);
        new_kx = rknew * st * cos(phi);
        new_ky = rknew * st * sin(phi);
        new_kz = rknew * ct;
        % Valley transition AFTER k-vector generation (critical for rand sequence!)
        rr = rand();
        if valley == 1
            if rr <= 0.5
                new_valley = 2;
            else
                new_valley = 3;
            end
        elseif valley == 2
            if rr <= 0.5
                new_valley = 1;
            else
                new_valley = 3;
            end
        elseif valley == 3
            if rr <= 0.5
                new_valley = 1;
            else
                new_valley = 2;
            end
        end
        
    case 3  % f-valley absorption
        new_energy = energy + params.w0_f - params.delta_fi_zero_f;
        rknew = params.smh * sqrt(new_energy * (1 + params.alpha_nonpara * new_energy));
        phi = 2 * pi * rand();
        ct = 1 - 2 * rand();
        st = sqrt(1 - ct^2);
        new_kx = rknew * st * cos(phi);
        new_ky = rknew * st * sin(phi);
        new_kz = rknew * ct;
        % Valley transition AFTER k-vector generation (critical for rand sequence!)
        rr = rand();
        if valley == 1
            if rr <= 0.5
                new_valley = 2;
            else
                new_valley = 3;
            end
        elseif valley == 2
            if rr <= 0.5
                new_valley = 1;
            else
                new_valley = 3;
            end
        elseif valley == 3
            if rr <= 0.5
                new_valley = 1;
            else
                new_valley = 2;
            end
        end
        
    case 4  % g-valley emission
        new_energy = energy - params.w0_g - params.delta_fi_zero_g;
        rknew = params.smh * sqrt(new_energy * (1 + params.alpha_nonpara * new_energy));
        phi = 2 * pi * rand();
        ct = 1 - 2 * rand();
        st = sqrt(1 - ct^2);
        new_kx = rknew * st * cos(phi);
        new_ky = rknew * st * sin(phi);
        new_kz = rknew * ct;
        % No valley change for g-type scattering
        
    case 5  % g-valley absorption
        new_energy = energy + params.w0_g - params.delta_fi_zero_g;
        rknew = params.smh * sqrt(new_energy * (1 + params.alpha_nonpara * new_energy));
        phi = 2 * pi * rand();
        ct = 1 - 2 * rand();
        st = sqrt(1 - ct^2);
        new_kx = rknew * st * cos(phi);
        new_ky = rknew * st * sin(phi);
        new_kz = rknew * ct;
        % No valley change for g-type scattering
        
    otherwise
        % Monolithic code: does nothing except a=1
        % Energy stays the same, k-vectors are NOT recalculated
        % Return current energy and DO NOT generate new k-vectors
        new_kx = NaN;  % Signal that k-vectors should not be updated
        new_ky = NaN;
        new_kz = NaN;
end

end

%% Local helper function: randomValleyTransition
function new_valley = randomValleyTransition(current_valley)
% RANDOMVALLEYTRANSITION - Randomly select new valley (equiprobable)
%
% For f-type scattering, electron transitions to one of the other two valleys

rr = rand();

if current_valley == 1
    if rr <= 0.5
        new_valley = 2;
    else
        new_valley = 3;
    end
elseif current_valley == 2
    if rr <= 0.5
        new_valley = 1;
    else
        new_valley = 3;
    end
else  % current_valley == 3
    if rr <= 0.5
        new_valley = 1;
    else
        new_valley = 2;
    end
end

end

%% Local helper function: calculateVelocity
function [vx, vy, vz] = calculateVelocity(valley, kx, ky, kz, energy, params)
% CALCULATEVELOCITY - Compute velocity components from k-vector and valley
%
% Inputs:
%   valley - Valley index (1, 2, or 3)
%   kx, ky, kz - Wave vector components
%   energy - Electron energy
%   params - Parameter structure
%
% Outputs:
%   vx, vy, vz - Velocity components

denom = 1 + 2 * params.alpha_nonpara * energy;

if valley == 1
    vx = params.hm1 * kx / denom;
    vy = params.hm2 * ky / denom;
    vz = params.hm3 * kz / denom;
elseif valley == 2
    vx = params.hm2 * kx / denom;
    vy = params.hm1 * ky / denom;
    vz = params.hm3 * kz / denom;
else  % valley == 3
    vx = params.hm3 * kx / denom;
    vy = params.hm2 * ky / denom;
    vz = params.hm1 * kz / denom;
end

end