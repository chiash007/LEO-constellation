% initial matrix definer
classdef Metrics
    properties
        Directory;
        TotalSimulations;
        Simulation;
        Estimators;
    end
    
    methods   
        
        function obj = Consolidate(Directory, TotalSimulations)
            obj.Directory = Directory;
            obj.TotalSimulations = TotalSimulations;
            load(strcat(Directory, 'Simulation.mat'));
            obj.Simulation = Simulation;
            obj.Estimators.Samples = obj.CalculateSamples();
        end
        
        function obj = set.Estimators(obj, val)
            obj.Estimators = val;
        end       

        function Estimators = CalculateSamples(obj)    
            for simulation_id = 1:obj.TotalSimulations              
                % traffic
                Results ={};
                fprintf('Consolidate Samples (%d/%d)\n', simulation_id, obj.TotalSimulations);
                load(strcat(obj.Directory, 'Results_', num2str(simulation_id), '.mat'));  
                for metric_id = 1:size(obj.Simulation.Metrics, 1)  
                    Metric = num2str(cell2mat(obj.Simulation.Metrics(metric_id)));   
                    for loopId = 1:size(Results, 2)  
                        % Served Demand
                        ServedDemand.(Metric)(loopId, :) ...
                            = mean(Results(loopId).(Metric)...
                            .Sources.serveddemand)./obj.Simulation.CBRs;   
                         % Total Blocked Sources
                         BlockedSources.(Metric)(loopId, :)  ...
                             = sum(Results(loopId).(Metric).Sources.serveddemand  == 0);
                         % Mean Delay
                         Delay.(Metric)(loopId, :) = sum(Results(loopId)...
                             .(Metric).Sources.propagationtime)./sum((Results(loopId)...
                             .(Metric).Sources.propagationtime > 0));          
                         % Mean Hops 
                         Hops.(Metric)(loopId, :) = sum(Results(loopId)...
                             .(Metric).Sources.hops)./sum((Results(loopId)....
                             .(Metric).Sources.hops > 0));                
                          % Total Saturated Links
                          SaturatedLinks.(Metric)(loopId, :)  = sum(Results(loopId)...
                              .(Metric).Sources.saturatedlinks)./size(find(obj.Simulation...
                              .Scenario.Propagation(loopId).Links > 0), 1);     
                          % Route Eclipse
                          RouteEclipse.(Metric)(loopId, :)  = sum(Results(loopId)...
                               .(Metric).Sources.routeeclipse)...
                               ./sum(Results(loopId).(Metric).Sources.routeeclipse > 0) ;                           
                    end
                    Estimators.(Metric).Delay(simulation_id, :) = mean(Delay.(Metric));
                    Estimators.(Metric).ServedDemand(simulation_id, :) = mean(ServedDemand.(Metric));
                    Estimators.(Metric).Hops(simulation_id, :) = mean(Hops.(Metric));
                    Estimators.(Metric).BlockedSources(simulation_id, :) = mean(BlockedSources.(Metric));
                    Estimators.(Metric).SaturatedLinks(simulation_id, :) = mean(SaturatedLinks.(Metric));
                    Estimators.(Metric).RouteEclipse(simulation_id, :) = mean(RouteEclipse.(Metric));    
                end  
            end        
        end         
    end
end
