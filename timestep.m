%models the orbital dynamics
classdef timestep
    properties (Constant)
        mu = 398600.5;              % Earth gravitational constant
        rp = 6371;                  % Planet Radius - km
        rs = 695500;                % Sun Radius - km
    end

    properties (Access = public)
        TotalSimulationTime;        % total simulation time in seconds
        SimulationInterval;         % simulation interval in seconds
        TotalSatellites = 66;       % Total number of satellites in the constellation
        TotalOrbitalPlanes = 6;     % Total number of orbital planes
        TotalSatellitesPerPlane = 11;% Total satellites per plane
        PolarBoundary = 75;         % Polar boundary latitude
        Zones;                      % Earth division into zones
        Constellation;              % Satellite constellation
        Propagation;                % Time propagation
        Eclipse;                    % Eclipse data
    end

    methods
        %% Constructor
        function obj = Scenario(Data, TotalSimulationTime, SimulationInterval)
            % Satellite constellation from TLE
            obj.Constellation = obj.ReadTLE();
            % Division of the Earth into zones
            obj.Zones = Zones();
            % Total simulation time in seconds
            obj.TotalSimulationTime = TotalSimulationTime;
            % Simulation interval in seconds
            obj.SimulationInterval = SimulationInterval;
            % Satellite propagation over time
            obj.Propagation = obj.Propagate(datevec(Data));
            % Eclipse data
            obj.Eclipse = obj.EclipseData();
        end

        %% Extracts data from the TLE file
        function SatelliteData = ReadTLE(obj)
            fd = fopen('TLE.txt', 'rb');
            A0 = fgetl(fd);
            A1 = fgetl(fd);
            A2 = fgetl(fd);
            SatelliteData = {};
            cont = 1;
            while ischar(A2)
                % Epoch of the data
                SatelliteData(cont).epoch = obj.EpochToData(str2num(A1(19:32)));
                % Eccentricity
                SatelliteData(cont).e = str2double(A2(27:33)) / (1e7);
                % Inclination
                SatelliteData(cont).i = str2double(A2(9:16));
                % Right Ascension of Ascending Node
                SatelliteData(cont).raan = str2double(A2(18:25));
                % Argument of Perigee
                SatelliteData(cont).w = str2double(A2(35:42));
                % Mean Anomaly
                SatelliteData(cont).M = str2double(A2(44:51));
                % Mean Motion
                n = str2double(A2(53:63));
                SatelliteData(cont).n = n;
                % Semi-Major Axis
                SatelliteData(cont).a = (obj.mu / (n * 2 * pi / (24 * 3600))^2)^(1/3);
                cont = cont + 1;
                A0 = fgetl(fd);
                A1 = fgetl(fd);
                A2 = fgetl(fd);
            end
            fclose(fd);
        end

        %% Returns date/time vector from TLE epoch
        function date_vector = EpochToData(obj, tle_epoch)
            ymd = floor(tle_epoch);
            yr = fix(ymd / 1000);
            dofyr = mod(ymd, 1000);
            if (yr < 57)
                year = yr + 2000;
            else
                year = yr + 1900;
            end
            decidy = round((tle_epoch - ymd) * 10^8) / 10^8;
            temp = decidy * 24;
            hh = fix(temp);
            temp = (temp - hh) * 60;
            mm = fix(temp);
            temp = (temp - mm) * 60;
            ss = floor(temp);
            nd = eomday(year, 1:12);
            temp = cumsum(nd);
            month = find(temp >= dofyr, 1);
            temp = temp(month) - dofyr;
            date = nd(month) - temp;
            date_vector = [year, month, date, hh, mm, ss];
        end

        %% Returns position vector in ECI from orbital elements
        function rECI = OrbitalElementsToECI(obj, SatNum, Data)
            DeltaT = etime(Data, obj.Constellation(SatNum).epoch);
            j2 = 1.0826359e-3; % J2 perturbation (Earth flattening)
            a = obj.Constellation(SatNum).a; % Semi-major axis
            e = obj.Constellation(SatNum).e; % Eccentricity
            i = obj.Constellation(SatNum).i; % Inclination
            raan = obj.Constellation(SatNum).raan; % Right Ascension of Ascending Node
            w = obj.Constellation(SatNum).w; % Argument of Perigee
            M = obj.Constellation(SatNum).M; % Mean Anomaly
            p = a * (1 - e^2); % Semi-latus rectum parameter
            T = 2 * pi * sqrt(a^3 / obj.mu); % Orbital period in seconds
            Tp = (M / 360 * T); % Time to pass through perigee
            ti = Tp + DeltaT; % Initial time of simulation considering perigee
            j2_raan = -(3/2) * j2 * (obj.rp / p)^2 * sqrt(obj.mu / (a^3)) * cosd(i); % J2 perturbation on RAAN
            j2_w = (3/4) * j2 * (obj.rp / p)^2 * sqrt(obj.mu / (a^3)) * (5 * cosd(i)^2 - 1); % J2 perturbation on argument of perigee
            % Keplerian motion
            n = sqrt(obj.mu / a^3);
            Mt = M + n * ti; % Mean anomaly at time t
            Et = obj.EccentricAnomaly(Mt, e); % Eccentric anomaly at time t
            theta = 2 * atan2(sqrt(1 + e) * sind(Et / 2), sqrt(1 - e) * cosd(Et / 2)); % True anomaly at time t
            % Position vector in perifocal coordinates
            rPQW = [p * cosd(Et) - a; p * sqrt(1 - e^2) * sind(Et); 0];
            % Transformation matrix from perifocal to ECI coordinates
            PQWtoECI = [cos(raan) * cos(w) - sind(raan) * sind(w) * cos(i), -cos(raan) * sind(w) - sind(raan) * cos(w) * cos(i), sind(raan) * sind(i);
                sind(raan) * cos(w) + cos(raan) * sind(w) * cos(i), -sind(raan) * sind(w) + cos(raan) * cos(w) * cos(i), -cos(raan) * sind(i);
                sind(w) * sind(i), cos(w) * sind(i), cos(i)];
            % Position vector in ECI coordinates
            rECI = PQWtoECI * rPQW;
            % Adding J2 perturbation
            rECI(1:2) = rECI(1:2) * (1 + j2 * obj.rp^2 / norm(rECI)^2 * (3/2 * sind(i)^2 - 1)); % Correcting for RAAN and argument of perigee
            rECI(1:2) = rECI(1:2) + j2_raan * ti; % RAAN perturbation
            rECI(1:2) = rECI(1:2) + j2_w * ti; % Argument of perigee perturbation
        end

        %% Propagate satellite position over time
        function PropagationData = Propagate(obj, TimeVector)
            TotalTimeSteps = floor(obj.TotalSimulationTime / obj.SimulationInterval);
            PropagationData = zeros(TotalTimeSteps, 3, obj.TotalSatellites);
            for t = 1:TotalTimeSteps
                for s = 1:obj.TotalSatellites
                    PropagationData(t, :, s) = obj.OrbitalElementsToECI(s, TimeVector);
                end
                TimeVector(6) = TimeVector(6) + obj.SimulationInterval;
                if TimeVector(6) >= 60
                    TimeVector(6) = TimeVector(6) - 60;
                    TimeVector(5) = TimeVector(5) + 1;
                    if TimeVector(5) >= 60
                        TimeVector(5) = TimeVector(5) - 60;
                        TimeVector(4) = TimeVector(4) + 1;
                        if TimeVector(4) >= 24
                            TimeVector(4) = TimeVector(4) - 24;
                            TimeVector(3) = TimeVector(3) + 1;
                        end
                    end
                end
            end
        end

        %% Calculates the eclipse data for each satellite
        function EclipseData = EclipseData(obj)
            TotalTimeSteps = floor(obj.TotalSimulationTime / obj.SimulationInterval);
            EclipseData = zeros(TotalTimeSteps, obj.TotalSatellites);
            for t = 1:TotalTimeSteps
                for s = 1:obj.TotalSatellites
                    EclipseData(t, s) = obj.IsSatelliteInEclipse(s, t);
                end
            end
        end

        %% Check if a satellite is in eclipse at a given time step
        function inEclipse = IsSatelliteInEclipse(obj, SatNum, TimeStep)
            % Calculate satellite position
            SatPos = squeeze(obj.Propagation(TimeStep, :, SatNum));
            % Check if the satellite is in Earth's shadow (umbra)
            inUmbra = obj.IsInUmbra(SatPos);
            if inUmbra
                inEclipse = true;
            else
                % Check if the satellite is in Earth's penumbra
                inPenumbra = obj.IsInPenumbra(SatPos);
                inEclipse = inPenumbra;
            end
        end

        %% Check if a point is in Earth's umbra
        function inUmbra = IsInUmbra(obj, Point)
            % Calculate distances from the point to the Earth and Sun centers
            DistEarth = norm(Point);
            DistSun = norm(Point - [0; obj.rs; 0]);
            % Check if the point is in Earth's umbra (fully in shadow)
            inUmbra = DistEarth < DistSun;
        end

        %% Check if a point is in Earth's penumbra
        function inPenumbra = IsInPenumbra(obj, Point)
            % Calculate distances from the point to the Earth and Sun centers
            DistEarth = norm(Point);
            DistSun = norm(Point - [0; obj.rs; 0]);
            % Check if the point is in Earth's penumbra (partially in shadow)
            inPenumbra = DistEarth > DistSun && DistEarth < DistSun + obj.rs;
        end

        %% Solves Kepler's equation for eccentric anomaly
        function EccentricAnomaly = EccentricAnomaly(obj, MeanAnomaly, Eccentricity)
            % Initial guess for Eccentric Anomaly
            EccentricAnomaly = MeanAnomaly;
            % Iteratively solve Kepler's equation
            for i = 1:100
                f = EccentricAnomaly - Eccentricity * sind(EccentricAnomaly) - MeanAnomaly;
                f_prime = 1 - Eccentricity * cosd(EccentricAnomaly);
                EccentricAnomaly = EccentricAnomaly - f / f_prime;
                % Break if convergence is reached
                if abs(f) < 1e-10
                    break;
                end
            end
        end
    end
end
