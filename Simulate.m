classdef Simulate

    properties
        Scenario;            % Simulation scenario data
        TotalTerminals;      % Number of terminals for simulation
        CBRs;                % Bit rates for simulation
        Metrics;             % Simulation metrics
        PowerCP;             % Power consumption in nominal operation
        PowerTX;             % Power consumption in transmission
        PowerRX;             % Power consumption in reception
        PowerON;             % Power consumption in normal operation
        PowerHarvesting;     % Power harvesting capability
        ISLCapacity;         % Link capacity
        Results;             % Results variable
        Lb;
        Lc;
        AdaptationType;
        Weights;
        ThresholdWeights;
        constantA;
    end
    
    methods
        
        % Constructor
        function obj = Simulate(Data, TotalSimulationTime, SimulationInterval)
            % Generate simulation scenario
            obj.Scenario = Scenario(Data, TotalSimulationTime, SimulationInterval);
        end

        function obj = set.Scenario(obj, val)
            obj.Scenario = val;
        end
        
        function obj = set.TotalTerminals(obj, val)
            obj.TotalTerminals = val;
        end     
        
        function obj = set.CBRs(obj, val)
            obj.CBRs = val;
        end  
        
        function obj = set.Metrics(obj, val)
            obj.Metrics = val;
        end  
        
        function obj = set.PowerCP(obj, val)
            obj.PowerCP = val;
        end   
        
        function obj = set.PowerON(obj, val)
            obj.PowerON = val;
        end 
        
        function obj = set.PowerTX(obj, val)
            obj.PowerTX = val;
        end        
        
        function obj = set.PowerRX(obj, val)
            obj.PowerRX = val;
        end
        
        function obj = set.PowerHarvesting(obj, val)
            obj.PowerHarvesting = val;
        end  
        
        function obj = set.ISLCapacity(obj, val)
            obj.ISLCapacity = val;
        end        
        
        function obj = set.Results(obj, val)
            obj.Results = val;
        end
        
        function obj = set.Lc(obj, val)
            obj.Lc = val;
        end
        
        function obj = set.Lb(obj, val)
            obj.Lb = val;
        end
        
        function Results = Routing(obj, loopId, Iterations)
                  
            if nargin == 2
                Iterations = size(obj.Scenario.Propagation, 2);
            end   
            Results = {};
            
            % ---------------- Simulation over time -------------------------------------
            LoopId = 1;
            
            while LoopId <= Iterations 
                fprintf('Simulation %d (Routing %d/%d)\n', loopId, LoopId, Iterations);       
                
                % ---------------------- Source/Destination Generation -----------------   
                Pairs = {};
                CkSrc = randsrc(obj.TotalTerminals, 1, 1:max(obj.Scenario.Zones.Continents(:)));
                
                for sourceId = 1: obj.TotalTerminals
                    % Source continent
                    Pairs(sourceId).SourceContinent = CkSrc(sourceId); 
                    IdCkSrc = find(obj.Scenario.Zones.Continents...
                        .*(obj.Scenario.Zones.HostDensity > 0) == CkSrc(sourceId));
                    ZkSrc = randsrc(1, 1, IdCkSrc');  
                    Pairs(sourceId).SourceZone = ZkSrc; 
                    
                    % Source satellites
                    SatSrc = obj.Scenario.Propagation(LoopId).Zones.Satellites(ZkSrc)';
                    Pairs(sourceId).SourceSatellite = SatSrc;
                    
                    % Destination continent
                    CkDst = randsrc(1, 1, [1:max(obj.Scenario.Zones.Continents(:));...
                        obj.Scenario.Propagation(LoopId).Zones.Flow(CkSrc(sourceId), :)]);                  
                    Pairs(sourceId).DestinationContinent = CkDst;      
                    
                    % Destination zone/satellite selection - avoid selecting the same
                    while(1)
                        IdCkDst = find(obj.Scenario.Zones.Continents...
                            .*(obj.Scenario.Zones.HostDensity > 0) == CkDst);
                        ZkDst = randsrc(1, 1, IdCkDst'); 
                        Pairs(sourceId).DestinationZone = ZkSrc;
                        SatDst = obj.Scenario.Propagation(LoopId).Zones.Satellites(ZkDst)';
                        Pairs(sourceId).DestinationSatellite = SatDst;
                        
                        % Avoid selecting the same satellite for source and destination
                        if SatSrc ~= SatDst
                            break;
                        end       
                    end
                end
                
                Results(LoopId).Sources = Pairs; 
                
                % ----------------------- Active Links ------------------------------
                MatDelay = ((obj.Scenario.Propagation(LoopId).Links...
                ./299792458).*1000000);   % propagation time matrix in milliseconds
                
                % Adjacency matrix of the network
                MatNetwork = (MatDelay > 0); % Binary matrix with 1 for MatDelay > 0; 0 otherwise
          
                % --------------------- Bit Rates --------------------------------
                for cbrId = 1: size(obj.CBRs, 2)
                    cbr = obj.CBRs(cbrId);
                    
                    % --------------------- Metrics ----------------------------------
                    for metricId = 1:size(obj.Metrics, 1)
                        Metric = char(obj.Metrics(metricId));   
                        
                        % Energy matrix
                        EnergyConsumption = zeros(obj.Scenario.TotalSatellites, 1);
                        LifeCycle = zeros(obj.Scenario.TotalSatellites, 1);
                        GD = zeros(obj.Scenario.TotalSatellites, 1);
                        DOD = zeros(obj.Scenario.TotalSatellites, 1);
                        DOD_ant = zeros(obj.Scenario.TotalSatellites, 1);
                        
                        if LoopId == 1 % Initial conditions  
                            InitialEnergy = ones(obj.Scenario.TotalSatellites, 1) * obj.PowerCP;
                            DOD = zeros(obj.Scenario.TotalSatellites, 1);
                        else % Non-initial conditions
                            % Initial energy for the time interval is equal to 
                            % the residual energy from the end of the previous interval
                            InitialEnergy = Results(LoopId-1).(Metric).FinalEnergy(:, cbrId);
                            DOD_ant =  Results(LoopId-1).(Metric).DOD(:, cbrId); % DOD from the previous iteration
                        end 
                        
                        Results(LoopId).(Metric).InitialEnergy(:, cbrId) = InitialEnergy;
                        Results(LoopId).(Metric).EnergyConsumption(:, cbrId) = EnergyConsumption;
                        
                        % Energy Capacity matrix already discounts the nominal operation energy of satellites
                        EnergyON = ones(obj.Scenario.TotalSatellites, 1)...
                            .*obj.PowerON.*obj.Scenario.SimulationInterval;
                        
                        % Current residual energy of the satellites
                        Energy =  InitialEnergy  -  EnergyON;
                        
                        % Update Energy Consumption
                        EnergyConsumption = EnergyON; 
                        
                        % --- Capacity, Delay, Energy - calculate - only once - common parameters
                        % for every couple of satellites
                        % Contiguous zone map
                        MapContig = zeros(size(MatNetwork));
                        for idx = 1:obj.Scenario.TotalSatellites
                            MapContig(idx, obj.Scenario.Propagation(LoopId).Zones.Contiguous(idx, :)) = 1;
                        end
                        
                        % --------------------- Formulation/Link Costs ------------------------------
                        [Costs, Coefficients] = Simulate.Formulation(CBRs(cbrId), MatDelay, MatNetwork, MapContig, Results(LoopId).Sources);
                        % --------------------- Energy Computation ------------------------------------
                        % Energy consumption on transmission
                        ETX = Simulate.ComputeEnergyTransmission(Energy, Costs.Tx, cbr, Results(LoopId).Sources, Metric, DOD_ant, obj);
                        
                        % Energy consumption on reception
                        ERX = Simulate.ComputeEnergyReception(Energy, Costs.Rx, cbr, Results(LoopId).Sources, Metric, DOD_ant, obj);
                        
                        % Energy consumption on idle mode
                        EON = Simulate.ComputeEnergyIdle(Energy, Costs.Idle, Results(LoopId).Sources, Metric, obj);
                        
                        % Energy consumption for harvesting
                        EHarvesting = Simulate.ComputeEnergyHarvesting(Energy, obj.PowerHarvesting, Results(LoopId).Sources, Metric, obj);
                        
                        % --------------------- Latency and Reliability Computation --------------------
                        [Results(LoopId).(Metric).Reliability(:, cbrId), Results(LoopId).(Metric).Latency(:, cbrId)] = ...
                            Simulate.ComputeLatencyReliability(Costs.Latency, Coefficients.Latency, Energy, cbr, Results(LoopId).Sources);
                        
                        % --------------------- Energy Total Calculation ------------------------------
                        Results(LoopId).(Metric).EnergyConsumption(:, cbrId) = EON + ETX + ERX + EHarvesting;
                        
                        % --------------------- Battery End-of-Life Calculation ------------------------
                        LifeCycle = Simulate.ComputeLifeCycle(Energy, EnergyConsumption, cbr, Results(LoopId).Sources, Metric, obj);
                        
                        % --------------------- End-of-Simulation Battery Status -----------------------
                        Results(LoopId).(Metric).FinalEnergy(:, cbrId) = InitialEnergy - EnergyConsumption;
                        Results(LoopId).(Metric).DOD(:, cbrId) = DOD_ant + EnergyConsumption ./ InitialEnergy;
                        Results(LoopId).(Metric).Lifetime(:, cbrId) = LifeCycle;
                        Results(LoopId).(Metric).GD(:, cbrId) = GD;
                    end
                end
                
                LoopId = LoopId + 1;
            end
        end
    end
    
    methods(Static)
        
        % Formulation of link costs and coefficients
        function [Costs, Coefficients] = Formulation(CBR, MatDelay, MatNetwork, MapContig, Sources)
            % Code for formulation...
        end
        
        % Compute energy consumption on transmission
        function ETX = ComputeEnergyTransmission(Energy, CostTX, CBR, Sources, Metric, DOD_ant, obj)
            % Code for computation...
        end
        
        % Compute energy consumption on reception
        function ERX = ComputeEnergyReception(Energy, CostRX, CBR, Sources, Metric, DOD_ant, obj)
            % Code for computation...
        end
        
        % Compute energy consumption on idle mode
        function EON = ComputeEnergyIdle(Energy, CostON, Sources, Metric, obj)
            % Code for computation...
        end
        
        % Compute energy consumption for harvesting
        function EHarvesting = ComputeEnergyHarvesting(Energy, PowerHarvesting, Sources, Metric, obj)
            % Code for computation...
        end
        
        % Compute latency and reliability
        function [Reliability, Latency] = ComputeLatencyReliability(CostLatency, CoeffLatency, Energy, CBR, Sources)
            % Code for computation...
        end
        
        % Compute battery end-of-life
        function LifeCycle = ComputeLifeCycle(Energy, EnergyConsumption, CBR, Sources, Metric, obj)
            % Code for computation...
        end
    end
end
