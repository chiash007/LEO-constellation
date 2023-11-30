% topology generator
clear all;
clc;

%% I - GENERATE SIMULATION SCENARIO: Date, Total Simulation Time in Seconds, Snapshot %%
% Simulation = Simulate('2018-06-01 12:30:00', 12200, 100);
% Simulate(Date, TotalSimulationTime, SimulationInterval);
% TotalSimulationTime: 1 lap = 6100 seconds; 2 laps = 12200 seconds; 5 laps = 30000 seconds;
% SimulationInterval: snapshot interval (scenarios), 100 seconds;
tic

Simulation = Simulate('2020-07-20 13:00:00', 90000, 100);

Simulation.Lb = 0.3;  % 0.0045 Minimum of 100 J in the Route;
Simulation.Lc = 0;  % 0.001
Simulation.AdaptationType = 0; % 2 - without considering the eclipse

Simulation.Weights = [0.35 0.35 0.30];
Simulation.ThresholdWeights = [0.15 0.70 0.15];
Simulation.constantA = 0.8;

Simulation.TotalSources = 1000;       % Total sources to simulate - number of terminals for simulation
Simulation.CBRs = [1];         % CBR Rates (Data Rates) to Simulate: in Mbps - Megabits per second
Simulation.Metrics = {'PROPOSAL'};   % Metrics
Simulation.PowerCP = 117000;       % Battery power capacity => 117 KJ
Simulation.PowerTX = 7;            % Power consumption in transmission  => watt = J/s
Simulation.PowerRX = 3;            % Power consumption in reception  => watt = J/s
Simulation.PowerON = 4;            % Power consumption in normal operation => watt = J/s
Simulation.PowerCG = 20;           % Energy harvesting power: 20 watt = 20 J/s
Simulation.LinkCapacity = 10;      % Link capacity (Mbps - Megabits per second)
%TotalSimulations = 30;               % Total simulations for confidence interval
TotalSimulations = 1;

%% II - SAVE SIMULATION DATA - SCENARIO
save(strcat('Results/', num2str(Simulation.TotalSources), '_sources', '/Simulation.mat'), 'Simulation', '-v7.3');

%% III - PERFORM ROUTING
for i = 1:TotalSimulations
    Results = Simulation.Routing(i); % Routing
    save(strcat('Results/', num2str(Simulation.TotalSources), '_sources', '/Results_', num2str(i), '.mat'), 'Results', '-v7.3');
end

%% IV - CONSOLIDATE DATA TO GENERATE GRAPHS
ConsolidatedData = Consolidate(strcat('Results/', num2str(Simulation.TotalSources), '_sources/'), TotalSimulations);

time = toc;
save(strcat('Results/', num2str(Simulation.TotalSources), '_sources', '/ConsolidatedData', '.mat'), 'ConsolidatedData', 'time', '-v7.3');
