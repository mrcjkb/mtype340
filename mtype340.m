classdef mtype340 < handle
    % MTYPE340: Model of a thermal storage tank based on the "MULTIPORT Store-Model" Type 340 (TRNSYS) 
    %          by H. Drueck
    %
    % Source: H. Drueck - MULTIPORT Store-Model for TRNSYS
    %         - Stratified fluid storage tank with four internal heat exchangers, 
    %           ten connections for direct charge and discharge and an internal electrical heater
    %           Type 340
    %           Version 1.99F, March 2006
    %           URL: http://www.trnsys.de/download/en/ts_type_340_en.pdf (Accessed 19.12.2016)
    %
    % For the purpose of simplification of the use, some of the input
    % arguments are set as constant properties in this model.
    %
    % Author:   Marc Jakobi, December 2016
    
    properties
        % Function handle for the ODE solver to be used. The
        % fastest solver depends on the operation type and on the
        % external time step size and can be taken from the
        % following list:
        % The operation types are:
        % [DP = double ports, HX = heat exchangers, AUX0 = auxiliary heater with var == 0]
        %  -> The useage of an auxiliary heater with ver == 1 has no influence on the operation type and speed
        %           
        %  Recommended solver function handles for various operation
        %  scenarios:
        %
        %  - only DP:
        %            @ode23
        %
        %  - DP and HX:
        %            @ode23 - for delt <= 650
        %            @ode23tb - for delt > 650
        %
        %  - DP and AUX0 (or only AUX0):
        %            @ode45
        %
        %  - DP, HX and AUX0:
        %            @ode45 - for delt < 120
        %            @ode113 - for 120 <= delt < 170
        %            @ode23 - for delt >= 170
        solver;
        %   Struct with the following fields:
        %       - aux.var = 1 for variable power, 0 for constant power
        %       - aux.pos = Relative height of the auxialiary heater -> 0 < aux.pos < 1
        %       - aux.T = Set temperature [°C] for the controller (in the case of constant power) 
        %                 [No dead-band temperatrue difference is set in this case]
        aux; 
        zs; % Kx1 vector with the relative positions of up to Nmax temperature sensors (NaN for sensors that are not used)
    end
    properties (SetAccess = 'protected')
        Tho; % 4x1 vector with the temperatures of the heat exchangers [°C]
        Qls; % Heat loss of the entire storage tank to the ambience [W]
        Qlbot; % Heat loss of the storage floor to the ambience [W]
        Qltop; % Heat loss of the storage lid to the ambience [W]
        Qlsx; % Nmax x 1 - vector with the ambient heat losses ofthe Nmax volume segments [W]
        Tdo; % Nx1 -  vector with the temperatures of the N double ports [°C]
        Qd; % Nx1 - vector with the heat flows through the N double ports [W]
        Ts; % Ox1 - vector with the temperatures at the O temperature sensor positions [°C]
        Taux; % Temperature at the position of the temperature controller [°C] (Taux = NaN, if aux(1) == 1)
        Qh; % 4x1 - vector with the thermal output through the M heat exchangers [W] (NaN for unused HX)
        Qhs; % 4x1 - vector with the thermal output between the M heat exchangers [W] (NaN for unused HX)
        Thm; % 4x1 - vector with the mean temperatures of the M heat exchangers [°C] (NaN for unused HX)
        % Utilized auxiliary power [W, el.]
        % (can be lower than the set power for aux.var == 0 if the
        %  auxiliary heater is repeatedly turned on and off)
        Qaux;
    end
    properties (SetAccess = 'immutable')
        Aq; % Cross-section area [m²] of the storage tank
        Hs; % Height [m] of the storage tank
        Vs; % Volume [m³] of the storage tank
        UA_a; % Nominal heat loss rate [W/K] of the storage tank
        UA_u; % Heat loss rate at the bottom of the storage tank [W/K]
        UA_o; % Heat loss rate at the top of the storage tank [W/K]
        UA_h; % Heat loss rates between storage tank and heat exchangers [W/K] (4x1-vector)
        Nmax; % Number of volume segments in the storage tank
        % Relative heights of the inputs of double port 1..N (Nx1-vector, zdi(i) <= 1)
        % N is the number of double ports
        zdi;
        % Relative heights of the outputs of double port 1..N (Nx1-vector, zdo(i) <= 1)
        % N is the number of double ports
        zdo;
        % Nx1-vector according to zdi & zdo to indicate stratified charging and non-stratified charging.
        % scdp(i) == 1 for stratified charging of the double port dp_i and
        % scdp(i) == 0 for non-stratified charging with fixed return flow positions.
        % If no double ports are used, set scdp(1:4) to NaN.
        scdp;
        % 4x1-vector with the relative heights of the HX-inputs of the HX 1..4 [0 <= zhi(i) <= 1]
        %       --> If one of the 4 possible HX is not used, set the
        %           respective positions zhi(i) and zho(i) to a negative value
        zhi; 
        % 4x1-vector with the relative heights of the HX-outputs of HX 1..4 [0 <= zho(i) <= 1]
        %       --> If one of the 4 possible HX is not used, set the
        %           respective positions zhi(i) and zho(i) to a negative value
        zho;
        % 4x1-vector with the volumina [m³] of the four HX
        % (zero for unused HX)
        Vh;
        % Time step size of the simulation [s]
        % e. g. 60 for a time step size of 1 minute
        delt; 
    end
    properties (Hidden, SetAccess = 'protected', GetAccess = 'protected')
        ini = true; % set to false after the initialization of Tz, hxdat, storedat and auxon
        Tz; % Nmax x 1 - vector with the temperatures in the zones x = 1..Nmax
        hxdat; % Struct containing the relevant information of the heat exchangers
        storedat;  % Struct containing the relevant information of the storage tank
        auxon; % true, when auxiliary heater is turned on or variable power is selected; otherwise false.
        ksiak;
        ksiav;
        auxpos;
    end
    properties (Hidden, Constant)
        cp_s = 4179.5; % Specific heat capacity of water at 40 °C (mean storage tank temperature) [Ws/(kg*K)] (= 4.1795 kJ/(kg*K))
        cp_h = [4179.5; 4179.5; 4179.5; 4179.5]; % Specific heat capacity of water in the heat exchangers
        rho_s = 992.21; % Density of water at 40 °C [kg/m³]
        % lambda_eff is determined via an experimental approach and varies
        % only slightly depending on the storage tank type and its
        % operation.
        % Determined values vary between 0.68 W/mK (very rarely) and 2.03 W/mK.
        % Common variations are between 1.3 and 1.8 W/mK. For this model, a
        % 1.52 W/mK was chosen.
        % SOURCES:
        %         - Bundesverband Solarwirtschaft e.V. - Mitarbeit der deutschen
        %           Solarindustrie bei der Überarbeitung der europäischen Normen für
        %           thermische Solaranlagen (Eurosol) [Tabelle 4.2]
        %         - H. Drück, E. Hahne - Kombispeicher auf dem Prüfstand [Tabelle 1]
        lambda = 1.52; % (5.46 kJ/(m*h*K))
        rho_h = 992.21; % Density of the fluid (water) in the heat exchangers [kg/m³] (Drinking water is assumed here)
        %Parameters
        % according to eq. 2.2 ("Heat transfer capacity rate (UA)*hx", p. 6 of the MULTIPORT Store-Model documentation
        %b1: Dependence of the heat loss rate between the HX and the
        %    storage tank on the mass flow (WVR_HX,s)
        %b2: Dependence of the WVR_HX,s on the temperature difference
        %    between the HX and the storage tank  
        %    (can be neglected - Drueck, Bachmann, Mueller-Steinhagen: 
        %      Testing of solar hot water stores by means of up- and down-scaling algorithms)
        %b3: Dependence of the WVR_HX,s on the mean temperature of the HX
        %b3: Abhängigkeit der WVR_WT,s von der mittleren Temperatur des WT
        %Mean values of the above journal article's determined values were chosen for b1 and b3
        %(b1 varies between 0.205 and 0.266 and b3 varies between 0.413 and
        %0.473 for volumes between 300 and 500 l).
        b_h = [0.237; 0; 0.509];
        % Factor for the start-up behaviour of the heat exchangers 
        % (Typical values are between 0.001.. 0.005) (see eq. 2.3 - Type 340 documentation)
        Shx = 0.025; 
    end
    
    methods
        % Constructor
        function t = mtype340(Hs, Vs, UA_a, UA_u, UA_o, UA_h, zdi, zdo, scdp, aux, zhi, zho, Vh, zs, delt, varargin)
            %MTYPE340: Initiates a model of a thermal storage tank based on the "MULTIPORT Store-Model" Type 340 (TRNSYS) 
            %          by H. Drueck
            %
            %   Syntax: s = mtype340(Hs, Vs, UA_a, UA_u, UA_o, UA_h, ...
            %                       zdi, zdo, scdp, aux, zhi, zho, Vh, zs, delt);
            %
            %           s = mtype340(Hs, Vs, UA_a, UA_u, UA_o, UA_h, ...
            %                       zdi, zdo, scdp, aux, zhi, zho, Vh, zs, delt, ...
            %                       'OptionName', 'OptionValue');
            %
            %   Input arguments:
            %
            %   Hs:     Height [m] of the storage tank
            %   Vs:     Volume [m³] of the storage tank
            %   UA_a:   Nominal heat loss rate [W/K] of the storage tank
            %   UA_u:   Heat loss rate at the bottom of the storage tank [W/K]
            %   UA_o:   Heat loss rate at the top of the storage tank [W/K]
            %   UA_h:   Heat loss rates between storage tank and heat exchangers [W/K] (4x1-vector)
            %   zdi:    Relative heights of the outputs of double port 1..N (Nx1-vector, zdo(i) <= 1)
            %           (N is the number of double ports)
            %   zdo:    Relative heights of the outputs of double port 1..N (Nx1-vector, zdo(i) <= 1)
            %   scdp:   Nx1-vector according to zdi & zdo to indicate stratified charging and non-stratified charging.
            %           scdp(i) == 1 for stratified charging of the double port dp_i and
            %           scdp(i) == 0 for non-stratified charging with fixed return flow positions.
            %           If no double ports are used, set scdp(1:4) to NaN.
            %   aux:    Struct with the following fields:
            %           - aux.var = 1 for variable power, 0 for constant power
            %           - aux.pos = Relative height of the auxialiary heater -> 0 < aux.pos < 1
            %           - aux.T = Set temperature [°C] for the controller (in the case of constant power) 
            %                 [No dead-band temperatrue difference is set in this case]
            %   zhi:    4x1-vector with the relative heights of the HX-outputs of HX 1..4 [0 <= zho(i) <= 1]
            %       --> If one of the 4 possible HX is not used, set the
            %           respective positions zhi(i) and zho(i) to a negative value
            %   zho:    4x1-vector with the relative heights of the HX-outputs of HX 1..4 [0 <= zho(i) <= 1]
            %       --> If one of the 4 possible HX is not used, set the
            %           respective positions zhi(i) and zho(i) to a negative value
            %   Vh:     4x1-vector with the volumina [m³] of the four HX
            %           (zero for unused HX)
            %   zs:     Kx1 vector with the relative positions of up to Nmax temperature sensors (NaN for sensors that are not used)
            %   delt:   Time step size of the simulation [s]
            %           e. g. 60 for a time step size of 1 minute
            %
            %
            %   Optional inputs ('OptionName', 'OptionValue' pairs)
            %
            %   'Nmax'   Number of volume segments in the storage tank
            %            (default: 10)
            %   'Solver' Function handle for the ODE-solver 
            %            (default: @ode23)
            
            % optional inputs
            p = inputParser;
            addOptional(p, 'Solver', @ode23)
            addOptional(p, 'Nmax', 10, @(x) isnumeric(x) & round(x) == x & floor(x) == x)
            parse(p, varargin{:})
            t.solver = p.Results.Solver;
            t.Nmax = p.Results.Nmax;
            if ~any([isequal(t.solver, @ode23); ...
                    isequal(t.solver, @ode45); ...
                    isequal(t.solver, @ode23tb); ...
                    isequal(t.solver, @ode113)])
                error('The solver must be one of the following: @ode23, @ode23tb, @ode45 or @ode113')
            end
            t.zdi = zdi; t.zdo = zdo; t.scdp = scdp;
            chk = isnan(zdi) + isnan(zdo) + isnan(scdp);
            if any(chk < 3 & chk > 0)
                error('zdi, zdo and scdp do not match')
            end
            if any([zdi > 1; zdo > 1; zdi < 0; zdo < 0])
                error('zdi and zdo must be between 0..1')
            end
            chk = zhi == 0 + zho == 0;
            if any(chk < 2 & chk > 0)
                error('zhi and zho do not match')
            end
            t.zhi = zhi; t.zho = zho;
            t.Hs = Hs; t.Vs = Vs; t.Aq = Vs./Hs;
            t.Nmax = Nmax;
            t.UA_a = UA_a; t.UA_u = UA_u; t.UA_o = UA_o; t.UA_h = UA_h;
            t.aux = aux; t.zhi = zhi; t.zho = zho;
            t.Vh = Vh; t.zs = zs; t.delt = delt;
            % Additional properties
            t.ksiak = zeros(Nmax,1); 
            t.ksiav = ksiak;
            t.auxpos = int32(ksiak);
        end
        
        function t = simulate(t, Tdi, mdotd, Thi, mdoth, Paux, Tamb)
            % SIMULATE:
        end
    end
    
end

