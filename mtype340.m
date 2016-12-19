classdef mtype340 < handle
    % MTYPE340: Model of a thermal storage tank based on the "MULTIPORT Store-Model" Type 340 (TRNSYS)
    %           by H. Drueck
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
    % mtype340 Methods:
    %
    %       mtype340 - Initiates a model of a thermal storage tank based on the "MULTIPORT Store-Model" Type 340 (TRNSYS)
    %       simulate - Pass simulation variables to MTYPE340 object and update.
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
        function ty = mtype340(Hs, Vs, UA_a, UA_u, UA_o, UA_h, zdi, zdo, scdp, aux, zhi, zho, Vh, zs, delt, varargin)
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
            ty.solver = p.Results.Solver;
            ty.Nmax = p.Results.Nmax;
            if ~any([isequal(ty.solver, @ode23); ...
                    isequal(ty.solver, @ode45); ...
                    isequal(ty.solver, @ode23tb); ...
                    isequal(ty.solver, @ode113)])
                error('The solver must be one of the following: @ode23, @ode23tb, @ode45 or @ode113')
            end
            ty.zdi = zdi; ty.zdo = zdo; ty.scdp = scdp;
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
            ty.zhi = zhi; ty.zho = zho;
            ty.Hs = Hs; ty.Vs = Vs; ty.Aq = Vs./Hs;
            ty.Nmax = Nmax;
            ty.UA_a = UA_a; ty.UA_u = UA_u; ty.UA_o = UA_o; ty.UA_h = UA_h;
            ty.aux = aux; ty.zhi = zhi; ty.zho = zho;
            ty.Vh = Vh; ty.zs = zs; ty.delt = delt;
            % Additional properties
            ty.ksiak = zeros(Nmax,1);
            ty.ksiav = ksiak;
            ty.auxpos = int32(ksiak);
            % initialization of hxdat & storedat, etc.
            ty.hxdat.nh14 = [];
            ty.hxdat.nh23 = [];
            ty.Qlsx = zeros(ty.Nmax,1); % ambient losses
            ty.Qlbot = 0; % Bottom losses
            ty.Qltop = 0; % Top losses
        end
        
        function ty = simulate(ty, Tdi, mdotd, Thi, mdoth, Tamb, Paux)
            % SIMULATE: Pass simulation variables to mtype340 object and update
            %
            % Input arguments:
            %
            %       - Tdi:      Nx1-vector with the input temperatures of the N used double ports [°C]
            %                   --> NaN for unused double ports
            %       - mdotd:    Nx1-vector with the mass flow rates of the N used double ports [kg/s]
            %                   --> [according to zdi & zdo] (positive values)
            %       - Thi:      4x1-vector with the input temperatures of the heat exchangers [°C]
            %                   --> NaN for unused heat exchangers
            %       - mdoth:    4x1-vector with the mass flow rates of the heat exchangers [kg/s]
            %                   --> NaN for unused heat exchangers
            %       - Tamb:     Ambient temperature [°C]
            %       - Paux:     Electrical power of auxiliary heater [W]
            %                   --> NaN if no AUX is used. An efficiency of
            %                   1 is assumed for the AUX in this model.
            
            if nargin == 6
                ty.Tdi = Tdi;
                ty.mdotd = mdotd;
                ty.Thi = Thi;
                ty.mdoth = mdoth;
                ty.Tamb = Tamb;
                ty.Paux = Paux;
            elseif nargin == 0
                Tdi = ty.Tdi;
                mdotd = ty.mdotd;
                Thi = ty.Thi;
                mdoth = ty.mdoth;
                Tamb = ty.Tamb;
                Paux = ty.Paux;
            else
                error('Wrong number of input arguments.')
            end
            irel_bndr = linspace(0,1,ty.Nmax+1); % relative heights of zone boundaries i
            % Auxiliary heater
            switch ty.aux.var
                case 0 % Constant power according to controller
                    % Determination of node that contains AUX
                    ty.auxpos = find(irel_bndr > ty.aux.pos,1)-1;
                    ty.ksiak(ty.auxpos) = 1; ty.ksiav(ty.auxpos) = 0;
                case 1 % variable AUX power
                    ty.auxpos = find(irel_bndr > ty.aux.pos,1)-1;
                    ty.ksiav(ty.auxpos) = 1; ty.ksiak(ty.auxpos) = 0;
                    ty.aux.T = nan;
                case 2 % no auxiliary heater
                    Paux = 0;
                    ty.aux.T = nan;
                    ty.auxpos = 1; % to avoid indexing errors
            end
            % Mathematical description [Eq. 4.1 in TRNSYS doc] -
            % pre-initialisation
            % Conversion of input-positions for the case of stratified DP dis/charging
            zdni = nan(size(ty.zdi));
            zdno = zdni;
            ks = (1:ty.Nmax)';
            msink = zeros(ty.Nmax,length(mdotd)); %log. vector for active mdot in each node
            for i = 1:length(ty.scdp)
                zdni(i) = find(irel_bndr > ty.zdi(i),1)-1;
                zdno(i) = find(irel_bndr > ty.zdo(i),1)-1;
                % find new virtual input node
                if ty.scdp(i) && mdotd(i) ~= 0
                    dTd = ty.Tz - Tdi(i);
                    chk = find(dTd > 0);
                    if ~isempty(chk)
                        zdni(i) = chk(1);
                    else
                        zdni(i) = ty.Nmax;
                    end
                end
                if zdni(i) < zdno(i)
                    mdotd(i) = abs(mdotd(i)); % upward flow direction
                    msink(:,i) = ks<=zdno(i) & ks>=zdni(i);
                else
                    mdotd(i) = 0-abs(mdotd(i)); % downward flow direction
                    msink(:,i) = ks>=zdno(i) & ks<=zdni(i);
                end
            end
            msink = logical(msink);
            [tf,idx] = ismember(ks,zdni);
            %______________logical switches ksi_____________________________________________
            % ksi1 = 1 for mdot_dp > 0, otherwise ksi1 = 0
            ksi1 = true(size(mdotd)); ksi1(mdotd<=0) = 0;
            % ksi2 = 1 for mdot_dp < 0, otherwise ksi2 = 0
            ksi2 = true(size(mdotd)); ksi2(mdotd>=0) = 0;
            % ksi3 = 1 für zones i, that are in contact with HX 1 oder 4 sind,
            % otherwise ksi3 = 0
            if ty.ini || ty.hxdat.nh14 ~= length(find([ty.zho(1),ty.zho(4)]>=0)) || ty.hxdat.nh23 ~= length(find([ty.zho(2),ty.zho(3)]>=0))
                % Find out if hx1 & hx4 are used alone, together or not at all
                numhx14 = length(find([ty.zho(1),ty.zho(4)]>=0));
                % Determine parallel temperature layers of HX 1 & 4
                Th14 = zeros(size(ty.Tz));
                switch numhx14
                    case 0 % neither are used
                        mm1 = []; mm4 = []; % Indexes for respective nodes
                        ksi3 = zeros(1,ty.Nmax);
                        nh14 = 1; % set to 1 to avoid dividing by zero
                        Th14 = ty.Tz;
                        h1nodes = []; h4nodes = [];
                    case 1 % only HX1 is used
                        ksi3 = zeros(1, ty.Nmax);
                        if ty.zho(1) >= 0
                            mm4 = [];
                            set = [find(irel_bndr > ty.zho(1),1)-1; find(irel_bndr > ty.zhi(1),1)-1];
                            % Used for charging or discharging?
                            if ty.zho(1) > ty.zhi(1) % Discharging
                                mm1 = set(2):set(1); l = length(mm1);
                                h1nodes = linspace(Thi(1),ty.Tz(set(1))-5,l+2);
                            else % Charging
                                mm1 = set(1):set(2); l = length(mm1);
                                h1nodes = linspace(ty.Tz(set(1))+5,Thi(1),l+2);
                            end
                            % Apply temperature of HX1 (idealised
                            % assumption of its course first)
                            Th14(mm1) = h1nodes(2:end-1);
                            Th14(mm1) = repmat(Thi(1),size(mm1)); % Init with constant HX temperature
                            ksi3(mm1) = 1;
                            h4nodes = [];
                        else % only HX4 is used
                            mm1 = [];
                            set = [find(irel_bndr > ty.zho(4),1)-1; find(irel_bndr > ty.zhi(4),1)-1];
                            if ty.zho(4) > ty.zhi(4) % Discharging
                                mm4 = set(2):set(1); l = length(mm4);
                                h4nodes = linspace(Thi(4),ty.Tz(set(1))-5,l+2);
                            else % Charging
                                mm4 = set(1):set(2); l = length(mm4);
                                h4nodes = linspace(ty.Tz(set(1))+5,Thi(4),l+2);
                            end
                            Th14(mm4) = h4nodes(2:end-1);
                            Th14(mm4) = repmat(Thi(4),size(mm4));
                            ksi3(mm4) = 1;
                            h1nodes = [];
                        end
                        nh14 = length(find(ksi3)); % Number of nodes that are occupied by HX1 & HX4
                    case 2 % HX1 & HX4 are both used
                        set = [find(irel_bndr > ty.zho(1),1)-1; find(irel_bndr > ty.zhi(1),1)-1;...
                            find(irel_bndr > ty.zho(4),1)-1; find(irel_bndr > ty.zhi(4),1)-1];
                        if ty.zho(1) > ty.zhi(1) % Discharging HX1
                            mm1 = set(2):set(1); l = length(mm1);
                            h1nodes = linspace(Thi(1),ty.Tz(set(1))-5,l+2);
                        else % Charging HX1
                            mm1 = set(1):set(2); l = length(mm1);
                            h1nodes = linspace(ty.Tz(set(1))+5,Thi(1),l+2);
                        end
                        Th14(mm1) = h1nodes(2:end-1);
                        Th14(mm1) = repmat(Thi(1),size(mm1));
                        if ty.zho(4) > ty.zhi(4) % Discharging HX4
                            mm4 = set(4):set(3); l = length(mm4);
                            h4nodes = linspace(Thi(4),ty.Tz(set(3))-5,l+2);
                        else % Charging HX4
                            mm4 = set(3):set(4); l = length(mm4);
                            h4nodes = linspace(ty.Tz(set(3))+5,Thi(4),l+2);
                        end
                        Th14(mm4) =  h4nodes(2:end-1);
                        Th14(mm4) = repmat(Thi(4),size(mm4));
                        set = [mm1,mm4];
                        ksi3 = zeros(1,ty.Nmax); ksi3(set) = 1;
                        nh14 = length(find(ksi3)); % Number of nodes that are occupied;
                end
                %ksi4 = 1 for zones i that are in contact with HX 2 or 3
                %otherwise ksi4 = 0
                %Determine if hx2 and hx3, are used alone, together or not
                %at all
                numhx23 = length(find([ty.zho(2),ty.zho(3)]>=0));
                Th23 = zeros(size(ty.Tz)); % Init
                switch numhx23
                    case 0 % Neither are used
                        mm2 = []; mm3 = [];
                        ksi4 = zeros(1,ty.Nmax);
                        nh23 = 1; % to avoid division by zero
                        Th23 = ty.Tz;
                        h2nodes = []; h3nodes = [];
                    case 1
                        ksi4 = zeros(1,ty.Nmax);
                        if ty.zho(2) >= 0 % Only HX2
                            mm3 = [];
                            set = [find(irel_bndr > ty.zho(2),1)-1; find(irel_bndr > ty.zhi(2),1)-1];
                            if ty.zho(2) > ty.zhi(2) % Discharge
                                mm2 = set(2):set(1); l = length(mm2);
                                h2nodes = linspace(Thi(2), ty.Tz(set(1))-5,l+2);
                            else % Charge
                                mm2 = set(1):set(2); l = length(mm2);
                                h2nodes = linspace(ty.Tz(set(1))+5,Thi(2),l+2);
                            end
                            % Idealised assumption HX2
                            Th23(mm2) = h2nodes(2:end-1);
                            Th23(mm2) = repmat(Thi(2),size(mm2));
                            ksi4(mm2) = 1;
                            h3nodes = [];
                        else %falls nur WT3 verwendet wird
                            mm2 = [];
                            set = [find(irel_bndr > ty.zho(3),1)-1; find(irel_bndr > ty.zhi(3),1)-1];
                            if ty.zho(3) > ty.zhi(3) % Discharge
                                mm3 = set(2):set(1); l = length(mm3);
                                h3nodes = linspace(Thi(3),ty.Tz(set(1))-5,l+2);
                            else %Charge
                                mm3 = set(1):set(2); l = length(mm3);
                                h3nodes = linspace(ty.Tz(set(1))+5,Thi(3),l+2);
                            end
                            %Idealised assumption HX3
                            Th23(mm3) = h3nodes(2:end-1);
                            Th23(mm3) = repmat(Thi(3),size(mm3));
                            ksi4(mm3) = 1;
                            h2nodes = [];
                        end
                        nh23 = length(find(ksi4)); % Number of nodes occupied by HX2 & HX3
                    case 2 % Both HX are used
                        set = [find(irel_bndr > ty.zho(2),1)-1; find(irel_bndr > ty.zhi(2),1)-1;...
                            find(irel_bndr > ty.zho(3),1)-1; find(irel_bndr > ty.zhi(3),1)-1];
                        if ty.zho(2) > ty.zhi(2) % Discharge HX2
                            mm2 = set(2):set(1); l = length(mm2);
                            h2nodes = linspace(Thi(2),ty.Tz(set(1))-5,l+2);
                        else % Charge HX2
                            mm2 = set(1):set(2); l = length(mm2);
                            h2nodes = linspace(ty.Tz(set(1))+5,Thi(2),l+2);
                        end
                        Th23(mm2) = h2nodes(2:end-1);
                        Th23(mm2) = repmat(Thi(2),size(mm2));
                        if ty.zho(3) > ty.zhi(3) % Discharge HX3
                            mm3 = set(4):set(3); l = length(mm3);
                            h3nodes = linspace(Thi(3),ty.Tz(set(3))-5,l+2);
                        else % Charge HX3
                            mm3 = set(3):set(4); l = length(mm3);
                            h3nodes = linspace(ty.Tz(set(3))+5,Thi(3),l+2);
                        end
                        Th23(mm3) =  h3nodes(2:end-1);
                        Th23(mm3) = repmat(Thi(3),size(mm3));
                        set = [mm2,mm3];
                        ksi4 = zeros(1,ty.Nmax); ksi4(set) = 1;
                        nh23 = length(find(ksi4));
                end
                ksi3 = ksi3'; ksi4 = ksi4';
                ty.hxdat.ksi3 = ksi3; ty.hxdat.ksi4 = ksi4;
                ty.hxdat.mm1 = mm1; ty.hxdat.mm2 = mm2; ty.hxdat.mm3 = mm3; ty.hxdat.mm4 = mm4;
                ty.hxdat.nh14 = nh14; ty.hxdat.nh23 = nh23;
                ty.hxdat.h1nodes = h1nodes; ty.hxdat.h2nodes = h2nodes; ty.hxdat.h3nodes = h3nodes; ty.hxdat.h4nodes = h4nodes;
            else
                ksi3 = ty.hxdat.ksi3; ksi4 = ty.hxdat.ksi4;
                mm1 = ty.hxdat.mm1; mm2 = ty.hxdat.mm2; mm3 = ty.hxdat.mm3; mm4 = ty.hxdat.mm4;
                nh14 = ty.hxdat.nh14; nh23 = ty.hxdat.nh23;
                h1nodes = ty.hxdat.h1nodes; h2nodes = ty.hxdat.h2nodes; h3nodes = ty.hxdat.h3nodes; h4nodes = ty.hxdat.h4nodes;
            end
            % Determine start-up of a heat exchanger
            % The inertia during the start-up of a HX is considered in the
            % solver, because it depends on the number of time steps. The
            % following code recognizes the start-up of a HX.
            if ty.ini
                ty.hxdat.mdoth = mdoth;
                ty.hxdat.mdoth(isnan(ty.hxdat.mdoth)) = 0;
                hxini = ones(4,1);
            else
                hxini = zeros(4,1);
                for i = 1:4
                    if ty.hxdat.mdoth(i) == 0 && mdoth(i) ~= 0
                        hxini(i) = 1;
                    end
                end
                ty.hxdat.mdoth = mdoth;
                ty.hxdat.mdoth(isnan(ty.hxdat.mdoth)) = 0;
            end
            % iteration of eq. 4.1 & eq. 4.2 & auxiliary heater for each node
            ksi5 = false(4,1);
            ksi5(mdoth > 0) = 1;
            ksi6 = ~ksi5;
            ksi5(mdoth == 0) = 0;
            ksi6(mdoth == 0) = 0;
            nodevect = zeros(ty.Nmax,1); nodevect([mm1, mm4]) = 1;
            rest = find(nodevect==0);
            Th14(rest) = ty.Tz(rest); % Equate to storage tank temperature so that deltaT = 0
            nodevect = zeros(ty.Nmax,1); nodevect([mm2, mm3]) = 1;
            rest = find(nodevect==0);
            Th23(rest) = ty.Tz(rest);
            if ty.ini
                Qh14_s = zeros(ty.Nmax,1); % Heat flow HX1/4-storage tank
                Qh23_s = Qh14_s; % Heat flow HX2/3-storage tank
            else
                ty.Qlsx = ty.storedat.Qlsx; % MTODO: remove Qlsx, etc from storedat
                Qh14_s = ty.storedat.Qh14_s;
                Qh23_s = ty.storedat.Qh23_s;
                ty.Qlbot = ty.storedat.Qlbot;
                ty.Qltop = ty.storedat.Qltop;
                Th14 = ty.storedat.Th14;
                Th23 = ty.storedat.Th23;
            end
            % Apply solver
            allTs = [ty.Tz; Th14; Th23; ty.Qlsx; Qh14_s; Qh23_s; ty.Qlbot; ty.Qltop; Paux];
            options = odeset('AbsTol',0.1); % Change the absolute error tolerance to 0,1
            [~,dT] = ty.solver(@tempsolver,[0 ty.delt], allTs, options);
            % Evaluation
            ty.Tz = dT(end,1:ty.Nmax)';
            Th14 = dT(end,ty.Nmax+1:2*ty.Nmax)';
            ty.storedat.Th14 = Th14;
            Th23 = dT(end,2*ty.Nmax+1:3*ty.Nmax)';
            ty.storedat.Th23 = Th23;
            % Heat flows through heat exchangers
            Qh1 = zeros(ty.Nmax,1); Qh2 = Qh1; Qh3 = Qh1; Qh4 = Qh1;
            Qh1(mm1) = (ty.Vh(1).*ty.cp_h(1).*ty.rho_h)./length(mm1).*Th14(mm1)./ty.delt;
            Qh2(mm2) = (ty.Vh(2).*ty.cp_h(2).*ty.rho_h)./length(mm2).*Th23(mm2)./ty.delt;
            Qh3(mm3) = (ty.Vh(3).*ty.cp_h(3).*ty.rho_h)./length(mm3).*Th23(mm3)./ty.delt;
            Qh4(mm4) = (ty.Vh(4).*ty.cp_h(4).*ty.rho_h)./length(mm4).*Th14(mm4)./ty.delt;
            ty.Qh = [sum(Qh1); sum(Qh2); sum(Qh3); sum(Qh4)];
            ty.Qh(ty.Qh == 0) = nan;
            ty.Qlsx = dT(end,3*ty.Nmax+1:4*ty.Nmax)'./ty.delt;
            ty.storedat.Qlsx = ty.Qlsx;
            Qh14_s = dT(end,4*ty.Nmax+1:5*ty.Nmax)'./ty.delt;
            ty.storedat.Qh14_s = Qh14_s;
            Qh23_s = dT(end,5*ty.Nmax+1:6*ty.Nmax)'./ty.delt;
            ty.storedat.Qh23_s = Qh23_s;
            ty.Qlbot = dT(end,6*ty.Nmax+1)./ty.delt;
            ty.storedat.Qlbot = ty.Qlbot;
            ty.Qltop = dT(end,6*ty.Nmax+2)./ty.delt;
            ty.storedat.Qltop = ty.Qltop;
            ty.Qaux = dT(end,6*ty.Nmax+3)./ty.delt;
            % Output temperatures of double ports
            ty.Tdo = nan(size(Tdi));
            if ~isnan(Tdi(1))
                set = ty.Tdo;
                for i = 1:length(set)
                    set(i) = find(irel_bndr > ty.zdo(i),1)-1;
                end
                ty.Tdo = ty.Tz(set);
            end
            % Heat flows
            ty.Qls = sum(ty.Qlsx)+ty.Qlbot+ty.Qltop; % total heat loss
            if ~isnan(Tdi(1))
                ty.Qd = abs(mdotd).*ty.cp_s.*(Tdi-ty.Tdo); % Heat flows through double ports
            else
                ty.Qd = nan;
            end
            % HX 1...4 - tank
            if ~isempty(mm1)
                Qh1s = sum(Qh14_s(mm1));
            else
                Qh1s = nan;
            end
            if ~isempty(mm4)
                Qh4s = sum(Qh14_s(mm4));
            else
                Qh4s = nan;
            end
            if ~isempty(mm2)
                Qh2s = sum(Qh23_s(mm2));
            else
                Qh2s = nan;
            end
            if ~isempty(mm3)
                Qh3s = sum(Qh23_s(mm3));
            else
                Qh3s = nan;
            end
            ty.Qhs = [Qh1s; Qh2s; Qh3s; Qh4s];
            % Extraction of the HX output temperatures
            ty.Tho = nan(4,1);
            ty.Thm = ty.Tho;
            if ~isempty(mm1)
                h1nodes(2:end-1) = Th14(mm1);
                chk = h1nodes ~= Thi(1);
                if chk(1) == 0
                    ty.Tho(1) = Thi(1)-abs(ty.Qhs(1)./(mdoth(1).*ty.cp_h(1)));
                    h1nodes(end) = ty.Tho(1);
                else
                    ty.Tho(1) = Thi(1)-abs(ty.Qhs(1)./(mdoth(1).*ty.cp_h(1)));
                    h1nodes(1) = ty.Tho(1);
                end
                ty.Thm(1) = mean(h1nodes);
                ty.hxdat.h1nodes = h1nodes;
            end
            if ~isempty(mm2)
                h2nodes(2:end-1) = Th23(mm2);
                chk = h2nodes ~= Thi(2);
                if chk(1) == 0
                    ty.Tho(2) = Thi(2)-abs(ty.Qhs(2)./(mdoth(2).*ty.cp_h(2)));
                    h2nodes(end) = ty.Tho(2);
                else
                    ty.Tho(2) = Thi(2)-abs(ty.Qhs(2)./(mdoth(2).*ty.cp_h(2)));
                    h2nodes(1) = ty.Tho(2);
                end
                ty.Thm(2) = mean(h2nodes);
                ty.hxdat.h2nodes = h2nodes;
            end
            if ~isempty(mm3)
                h3nodes(2:end-1) = Th23(mm3);
                chk = h3nodes ~= Thi(3);
                if chk(1) == 0
                    ty.Tho(3) = Thi(3)-abs(ty.Qhs(3)./(mdoth(3).*ty.cp_h(3)));
                    h3nodes(end) = ty.Tho(3);
                else
                    ty.Tho(3) = Thi(3)-abs(ty.Qhs(3)./(mdoth(3).*ty.cp_h(3)));
                    h3nodes(1) = ty.Tho(3);
                end
                ty.Thm(3) = mean(h3nodes);
                ty.hxdat.h3nodes = h3nodes;
            end
            if ~isempty(mm4)
                h4nodes(2:end-1) = Th14(mm4);
                chk = h4nodes ~= Thi(4);
                if chk(1) == 0
                    ty.Tho(4) = Thi(4)-abs(ty.Qhs(4)./(mdoth(4).*ty.cp_h(4)));
                    h4nodes(end) = ty.Tho(4);
                else
                    ty.Tho(4) = Thi(4)-abs(ty.Qhs(4)./(mdoth(4).*ty.cp_h(4)));
                    h4nodes(1) = ty.Tho(4);
                end
                ty.Thm(4) = mean(h4nodes);
                ty.hxdat.h4nodes = h4nodes;
            end
            % Temperature sensor outputs
            ty.Ts = nan(size(ty.zs));
            if ~isnan(ty.zs(1))
                set = ty.Ts;
                for i = 1:length(ty.zs)
                    set(i) = find(irel_bndr > ty.zs(i),1)-1;
                end
                ty.Ts = ty.Tz(set);
            end
            if ty.aux.var ~= 2
                ty.Taux = ty.Tz(ty.auxpos);
            else
                ty.Taux = nan;
            end
            function Tdot = tempsolver(t, T)
                % TEMPSOLVER: ODE-solver for the calculation of the temperature changes
                %             The energy balances can be calculated from this.
                Tdot = zeros(6.*ty.Nmax+3,1);
                %Ts: Temperatures in the storage tank
                %Th14: Temperatures in HX 1 & 4
                %Th23: Temperatures in HX 2 & 3
                % Factor for the time dependency of UA_hs has to be calculated 4 times. 1 time for each HX
                if sum(isnan(mdoth)) < 4 && sum(hxini) > 0
                    NITS = zeros(1,4);
                    for x = 1:4
                        NITS(x) = t./ty.delt;
                    end
                    
                    Fhx = ones(4,1);
                    for x = 1:4
                        if ~isnan(mdoth(x)) && hxini(x)
                            Fhxt = linspace(0,1./ty.Shx,1000); % Fhx is between  0 and 1/Shx (temporary crutch)
                            if NITS(x) == 0
                                c = mdoth(x).*ty.delt;
                            else
                                c = mdoth(x).*NITS(x).*ty.delt;
                            end
                            tmp = (1./ty.Shx .* Fhxt) + (c(x) + (1 - Fhxt)); % eq.1: tmp & Fhx not known, bandwidth of Fhx known
                            chk = abs(tmp-Fhxt./ty.Shx); % eq.2: tmp = Fhx./Shx [min(chk) is the approx. intersection between eq.1 & eq.2]
                            tmp = tmp(chk == min(chk)); % Nx1-vector mit tmp for each moment in time
                            % tmp must be <= 1./Shx and >= 0
                            tmp(tmp > 1./ty.Shx) = 1./ty.Shx;
                            tmp(tmp < 0) = 0;
                            tmp = tmp(1);
                            Fhx(x) = tmp.*ty.Shx; % Factos for the time dependency of UA_hx for the current time step
                            
                        end
                    end
                else
                    Fhx = ones(4,1);
                end
                % UA_hs [(UA)*_hx,s - time-, mass flow- and temperature-dependent heat loss rates between HX and storage tank
                % Determination of temperature at positions of HX inputs
                % Differences between T_HXin & T_Snode -->
                dThin_s1 = Thi(1) - T(mm1); dThin_s4 = Thi(4) - T(mm4); dThin_s2 = Thi(2) - T(mm2); dThin_s3 = Thi(3) - T(mm3);
                % Sums of             ----"----------- -->
                sThin_s1 = Thi(1) + T(mm1); sThin_s4 = Thi(4) + T(mm4); sThin_s2 = Thi(2) + T(mm2); sThin_s3 = Thi(3) + T(mm3);
                % Initializations
                F14 = zeros(size(ksi3)); F14(mm1) = Fhx(1); F14(mm4) = Fhx(4);
                dThin_s14 = zeros(size(ksi3)); dThin_s14(mm1) = dThin_s1; dThin_s14(mm4) = dThin_s4;
                sThin_s14 = zeros(size(ksi3)); sThin_s14(mm1) = sThin_s1; sThin_s14(mm4) = sThin_s4;
                mdoth14 = zeros(size(ksi3)); mdoth14(mm1) = mdoth(1); mdoth14(mm4) = mdoth(4);
                % differential equations
                vals = ty.UA_h(1) .* F14 .* abs(mdoth14) .^ ty.b_h(1) .* dThin_s14 .^ ty.b_h(2) .* (sThin_s14./2) .^ ty.b_h(3);
                UA_h14 = zeros(size(ksi3));
                UA_h14(mm1) = vals(mm1);
                vals = ty.UA_h(4) .* F14 .* abs(mdoth14) .^ ty.b_h(1) .* dThin_s14 .^ ty.b_h(2) .* (sThin_s14./2) .^ ty.b_h(3);
                UA_h14(mm4) = vals(mm4);
                F23 = zeros(size(ksi4)); F23(mm2) = Fhx(2); F23(mm3) = Fhx(3);
                dThin_s23 = zeros(size(ksi4)); dThin_s23(mm2) = dThin_s2; dThin_s23(mm3) = dThin_s3;
                sThin_s23 = zeros(size(ksi4)); sThin_s23(mm2) = sThin_s2; sThin_s23(mm3) = sThin_s3;
                mdoth23 = zeros(size(ksi4)); mdoth23(mm2) = mdoth(2); mdoth23(mm3) = mdoth(3);
                vals = ty.UA_h(2) .* F23 .* abs(mdoth23) .^ ty.b_h(1) .* dThin_s23 .^ ty.b_h(2) .* (sThin_s23./2) .^ ty.b_h(3);
                UA_h23 = zeros(size(ksi4));
                UA_h23(mm2) = vals(mm2);
                vals = ty.UA_h(3) .* F23 .* abs(mdoth23) .^ ty.b_h(1) .* dThin_s23 .^ ty.b_h(2) .* (sThin_s23./2) .^ ty.b_h(3);
                UA_h23(mm3) = vals(mm3);
                % Energy balance for the mixing in the storage tank
                if tf(1) % all mdotd except for dp one -> add zero to vector in case of empty matrix; DP-mdotd separate
                    Omdot = [0; mdotd(1:idx(1)-1); mdotd(idx(1)+1:end)]; Omsink = msink(1,:); Omsink = [true, Omsink(1:idx(1)-1), Omsink(idx(1)+1:end)];
                    Oksi2 = [0; ksi2(1:idx(1)-1); ksi2(idx(1)+1:end)]; % O: represent other double ports, not connected to input
                    Tdot(1) = ty.Nmax ./ (ty.Vs .* ty.rho_s .* ty.cp_s).*(ty.ksiak(1) .* Paux .* (T(1) <= ty.aux.T) + ty.ksiav(1) .* Paux...
                        + sum(Omdot(Omsink) .* ty.cp_s .* Oksi2(Omsink) .* (T(1) - T(2)))...
                        + mdotd(idx(1)) .* ty.cp_s .* ksi1(idx(1)) .* (Tdi(idx(1)) - T(1))...
                        + ksi3(1) .* UA_h14(1) ./ nh14 .* (T(ty.Nmax+1) - T(1)) + ksi4(1) .* UA_h23(1) ./ nh23 .* (T(2.*ty.Nmax+1) - T(1))...
                        + ty.lambda .* ty.Aq./ty.Hs .* ty.Nmax .* (T(2) - T(1)) - (ty.UA_a ./ ty.Nmax + ty.UA_u).*(T(1) - Tamb)); % Paux = auxiliary part
                else
                    Tdot(1) = ty.Nmax ./ (ty.Vs .* ty.rho_s .* ty.cp_s) .* (ty.ksiak(1) .* Paux .* (T(1) <= ty.aux.T) + ty.ksiav(1) .* Paux...
                        + sum(mdotd(msink(1,:)) .* ty.cp_s .* (ksi2(msink(1,:)) .* (T(1) - T(2))))...
                        + ksi3(1) .* UA_h14(1) ./ nh14 .* (T(ty.Nmax+1) - T(1)) + ksi4(1) .* UA_h23(1) ./ nh23 .* (T(2.*ty.Nmax+1) - T(1))...
                        + ty.lambda .* ty.Aq ./ ty.Hs .* ty.Nmax .* (T(2) - T(1)) - (ty.UA_a ./ ty.Nmax + ty.UA_u).*(T(1) - Tamb));
                end
                for k = 2:ty.Nmax-1
                    if tf(k)
                        Omdot = [0; mdotd(1:idx(k)-1); mdotd(idx(k)+1:end)]; Omsink = msink(k,:); Omsink = [true, Omsink(1:idx(k)-1), Omsink(idx(k)+1:end)];
                        Oksi1 = [0; ksi1(1:idx(k)-1); ksi1(idx(k)+1:end)]; Oksi2 = [0; ksi2(1:idx(k)-1); ksi2(idx(k)+1:end)];
                        Tdot(k) = ty.Nmax ./ (ty.Vs .* ty.rho_s .* ty.cp_s) .* (ty.ksiak(k).*Paux.*(T(k) <= ty.aux.T) + ty.ksiav(k) .* Paux...
                            + sum(Omdot(Omsink) .* ty.cp_s .* (Oksi1(Omsink) .* (T(k-1) - T(k)) + Oksi2(Omsink) .* (T(k) - T(k+1))))...
                            + mdotd(idx(k)) .* ty.cp_s .* (ksi1(idx(k)) .* (Tdi(idx(k)) - T(k)) + ksi2(idx(k)) .* (T(k) - Tdi(idx(k))))...
                            + ksi3(k) .* UA_h14(k) ./ nh14 .* (T(ty.Nmax+k) - T(k)) + ksi4(k) .* UA_h23(k) ./ nh23 .* (T(2.*ty.Nmax+k) - T(k))...
                            + ty.lambda .* ty.Aq ./ ty.Hs .* ty.Nmax .* ((T(k+1) - T(k)) + (T(k-1) - T(k))) - ty.UA_a ./ ty.Nmax .* (T(k) - Tamb));
                    else
                        Tdot(k) = ty.Nmax ./ (ty.Vs .* ty.rho_s .* ty.cp_s) .* (ty.ksiak(k) .* Paux.*(T(k) <= ty.aux.T) + ty.ksiav(k) .* Paux...
                            + sum(mdotd(msink(k,:)) .* ty.cp_s .* (ksi1(msink(k,:)) .* (T(k-1) - T(k)) + ksi2(msink(k,:)) .* (T(k) - T(k+1))))...
                            + ksi3(k) .* UA_h14(k) ./ nh14 .* (T(ty.Nmax+k) - T(k)) + ksi4(k) .* UA_h23(k) ./ nh23 .* (T(2.*ty.Nmax+k) - T(k))...
                            + ty.lambda .* ty.Aq./ty.Hs .* ty.Nmax .* ((T(k+1) - T(k)) + (T(k-1) - T(k))) - ty.UA_a ./ ty.Nmax .* (T(k) - Tamb));
                    end
                end
                if tf(ty.Nmax)
                    Omdot = [0; mdotd(1:idx(ty.Nmax)-1); mdotd(idx(ty.Nmax)+1:end)]; Omsink = msink(ty.Nmax,:); Omsink = [true, Omsink(1:idx(ty.Nmax)-1), Omsink(idx(ty.Nmax)+1:end)];
                    Oksi1 = [0; ksi1(1:idx(ty.Nmax)-1);ksi1(idx(ty.Nmax)+1:end)];
                    Tdot(ty.Nmax) = ty.Nmax ./ (ty.Vs .* ty.rho_s .* ty.cp_s) .* (ty.ksiak(ty.Nmax) .* Paux .* (T(ty.Nmax) <= ty.aux.T) + ty.ksiav(ty.Nmax) .* Paux...
                        + sum(Omdot(Omsink) .* ty.cp_s .* Oksi1(Omsink) .* (T(ty.Nmax-1) - T(ty.Nmax)))...
                        + mdotd(idx(ty.Nmax)) .* ty.cp_s .* ksi2(idx(ty.Nmax)) .* (T(ty.Nmax) - Tdi(idx(ty.Nmax)))...
                        + ksi3(ty.Nmax) .* UA_h14(ty.Nmax) ./ nh14 .* (T(ty.Nmax+ty.Nmax) - T(ty.Nmax)) + ksi4(ty.Nmax) .* UA_h23(ty.Nmax) ./ nh23 .* (T(2.*ty.Nmax+ty.Nmax) - T(ty.Nmax))...
                        + ty.lambda .* ty.Aq ./ ty.Hs .* ty.Nmax .* (T(ty.Nmax-1) - T(ty.Nmax)) - (ty.UA_a ./ ty.Nmax + ty.UA_o) .* (T(ty.Nmax) - Tamb));
                else
                    Tdot(ty.Nmax) = ty.Nmax ./ (ty.Vs .* ty.rho_s .* ty.cp_s).*(ty.ksiak(ty.Nmax) .* Paux .* (T(ty.Nmax) <= ty.aux.T) + ty.ksiav(ty.Nmax) .* Paux...
                        + sum(mdotd(msink(k,:)) .* ty.cp_s .* ksi1(msink(k,:)) .* (T(ty.Nmax-1) - T(ty.Nmax)))...
                        + ksi3(ty.Nmax) .* UA_h14(ty.Nmax) ./ nh14 .* (T(ty.Nmax+ty.Nmax) - T(ty.Nmax)) + ksi4(ty.Nmax) .* UA_h23(ty.Nmax) ./ nh23 .* (T(2.*ty.Nmax+ty.Nmax) - T(ty.Nmax))...
                        + ty.lambda .* ty.Aq./ty.Hs .* ty.Nmax .* (T(ty.Nmax-1) - T(ty.Nmax)) - (ty.UA_a ./ ty.Nmax + ty.UA_o) .* (T(ty.Nmax) - Tamb));
                end
                % Energy balance for heat exchanger 1
                if ~isempty(mm1)
                    Tdot(ty.Nmax+mm1(1)) = length(mm1)./(ty.Vh(1) .* ty.rho_h .* ty.cp_h(1))...
                        .* (ksi5(1) .* mdoth(1) .* ty.cp_h(1) .* (Thi(1) - T(ty.Nmax+mm1(1))) + ksi6(1) .* mdoth(1) .* ty.cp_h(1) .* (T(ty.Nmax+mm1(1)) - T(ty.Nmax+mm1(2)))...
                        + UA_h14(mm1(1)) ./ length(mm1) .* (T(mm1(1)) - T(ty.Nmax+mm1(1))));
                    for k = 2:length(mm1)-1
                        Tdot(ty.Nmax+mm1(k)) = length(mm1) ./ (ty.Vh(1) .* ty.rho_h .* ty.cp_h(1))...
                            .* (ksi5(1) .* mdoth(1) .* ty.cp_h(1) .* (T(ty.Nmax+mm1(k-1)) - T(ty.Nmax+mm1(k))) + ksi6(1) .* mdoth(1) .* ty.cp_h(1) .* (T(ty.Nmax+mm1(k)) - T(ty.Nmax+mm1(k+1)))...
                            + UA_h14(mm1(k)) ./ length(mm1) .* (T(mm1(k)) - T(ty.Nmax+mm1(k))));
                    end
                    Tdot(ty.Nmax+mm1(end)) = length(mm1) ./ (ty.Vh(1) .* ty.rho_h .* ty.cp_h(1))...
                        .* (ksi5(1) .* mdoth(1) .* ty.cp_h(1) .* (T(ty.Nmax+mm1(end-1)) - T(ty.Nmax+mm1(end))) + ksi6(1) .* mdoth(1) .* ty.cp_h(1) .* (T(ty.Nmax+mm1(end)) - Thi(1))...
                        + UA_h14(mm1(end)) ./ length(mm1) .* (T(mm1(end)) - T(ty.Nmax+mm1(end))));
                end
                if ~isempty(mm4)
                    % Energy balance for heat exchanger 4
                    Tdot(ty.Nmax+mm4(1)) = length(mm4) ./ (ty.Vh(4) .* ty.rho_h .* ty.cp_h(4))...
                        .* (ksi5(4) .* mdoth(4) .* ty.cp_h(4) .* (Thi(4) - T(ty.Nmax+mm4(1))) + ksi6(4) .* mdoth(4) .* ty.cp_h(4) .* (T(ty.Nmax+mm4(1)) - T(ty.Nmax+mm4(2)))...
                        + UA_h14(mm4(1)) ./ length(mm4) .* (T(mm4(1)) - T(ty.Nmax+mm4(1))));
                    for k = 2:length(mm4)-1
                        Tdot(ty.Nmax+mm4(k)) = length(mm4) ./ (ty.Vh(4) .* ty.rho_h .* ty.cp_h(4))...
                            .* (ksi5(4) .* mdoth(4) .* ty.cp_h(4) .* (T(ty.Nmax+mm4(k-1)) - T(ty.Nmax+mm4(k))) + ksi6(4) .* mdoth(4) .* ty.cp_h(4) .* (T(ty.Nmax+mm4(k)) - T(ty.Nmax+mm4(k+1)))...
                            + UA_h14(mm4(k)) ./ length(mm4) .* (T(mm4(k)) - T(ty.Nmax+mm4(k))));
                    end
                    Tdot(ty.Nmax+mm4(end)) = length(mm4) ./ (ty.Vh(4) .* ty.rho_h .* ty.cp_h(4))...
                        .* (ksi5(4) .* mdoth(4) .* ty.cp_h(4) .* (T(ty.Nmax+mm4(end-1)) - T(ty.Nmax+mm4(end))) + ksi6(4) .* mdoth(4) .* ty.cp_h(4) .* (T(ty.Nmax+mm4(end)) - Thi(4))...
                        + UA_h14(mm4(end)) ./ length(mm4) .* (T(mm4(end)) - T(ty.Nmax+mm4(end))));
                end
                if ~isempty(mm2)
                    % Energy balances for heat exchanger 2
                    Tdot(2.*ty.Nmax+mm2(1)) = length(mm2) ./ (ty.Vh(2) .* ty.rho_h .* ty.cp_h(2))...
                        .* (ksi5(2) .* mdoth(2) .* ty.cp_h(2) .* (Thi(2) - T(2.*ty.Nmax+mm2(1))) + ksi6(2) .* mdoth(2) .* ty.cp_h(2) .* (T(2.*ty.Nmax+mm2(1)) - T(2.*ty.Nmax+mm2(2)))...
                        + UA_h23(mm2(1)) ./ length(mm2) .* (T(mm2(1)) - T(2 .* ty.Nmax+mm2(1))));
                    for k = 2:length(mm2)-1
                        Tdot(2 .* ty.Nmax+mm2(k)) = length(mm2) ./ (ty.Vh(2) .* ty.rho_h .* ty.cp_h(2))...
                            .* (ksi5(2) .* mdoth(2) .* ty.cp_h(2) .* (T(2.*ty.Nmax+mm2(k-1)) - T(2.*ty.Nmax+mm2(k))) + ksi6(2) .* mdoth(2) .* ty.cp_h(2) .* (T(2.*ty.Nmax+mm2(k)) - T(2.*ty.Nmax+mm2(k+1)))...
                            + UA_h23(mm2(k)) ./ length(mm2) .* (T(mm2(k)) - T(2.*ty.Nmax+mm2(k))));
                    end
                    Tdot(2.*ty.Nmax+mm2(end)) = length(mm2) ./ (ty.Vh(2) .* ty.rho_h .* ty.cp_h(2))...
                        .* (ksi5(2).*mdoth(2) .* ty.cp_h(2) .* (T(2.*ty.Nmax+mm2(end-1)) - T(2.*ty.Nmax+mm2(end))) + ksi6(2) .* mdoth(2) .* ty.cp_h(2) .* (T(2.*ty.Nmax+mm2(end)) - Thi(2))...
                        + UA_h23(mm2(end)) ./ length(mm2) .* (T(mm2(end)) - T(2.*ty.Nmax+mm2(end))));
                end
                if ~isempty(mm3)
                    % Energy balances for heat exchanger 3
                    Tdot(2.*ty.Nmax+mm3(1)) = length(mm3) ./ (ty.Vh(3) .* ty.rho_h .* ty.cp_h(3))...
                        .* (ksi5(3) .* mdoth(3) .* ty.cp_h(3) .* (Thi(3) - T(2 .* ty.Nmax + mm3(1))) + ksi6(3) .* mdoth(3) .* ty.cp_h(3) .* (T(2.*ty.Nmax+mm3(1)) - T(2.*ty.Nmax+mm3(2)))...
                        + UA_h23(mm3(1)) ./ length(mm3) .* (T(mm3(1)) - T(2.*ty.Nmax+mm3(1))));
                    for k = 2:length(mm3)-1
                        Tdot(2.*ty.Nmax+mm3(k)) = length(mm3) ./ (ty.Vh(3) .* ty.rho_h .* ty.cp_h(3))...
                            .* (ksi5(3) .* mdoth(3) .* ty.cp_h(3) .* (T(2.*ty.Nmax+mm3(k-1)) - T(2.*ty.Nmax+mm3(k))) + ksi6(3) .* mdoth(3) .* ty.cp_h(3) .* (T(2.*ty.Nmax+mm3(k)) - T(2.*ty.Nmax+mm3(k+1)))...
                            + UA_h23(mm3(k)) ./ length(mm3) .* (T(mm3(k)) - T(2.*ty.Nmax+mm3(k))));
                    end
                    Tdot(2.*ty.Nmax+mm3(end)) = length(mm3) ./ (ty.Vh(3) .* ty.rho_h .* ty.cp_h(3))...
                        .* (ksi5(3) .* mdoth(3) .* ty.cp_h(3) .* (T(2 .* ty.Nmax+mm3(end-1)) - T(2.*ty.Nmax+mm3(end)))+ksi6(3) .* mdoth(3) .* ty.cp_h(3).*(T(2.*ty.Nmax+mm3(end)) - Thi(3))...
                        + UA_h23(mm3(end)) ./ length(mm3) .* (T(mm3(end)) - T(2.*ty.Nmax+mm3(end))));
                end
                % Heat flow rates: Losses storage tank - ambience
                Tdot(3.*ty.Nmax+1) = -(ty.UA_a ./ ty.Nmax) .* (T(1) - Tamb);
                for k = 2:ty.Nmax-1
                    Tdot(3 .* ty.Nmax+k) = -ty.UA_a ./ ty.Nmax .* (T(k) - Tamb);
                end
                Tdot(3.*ty.Nmax + ty.Nmax) = -(ty.UA_a ./ ty.Nmax) .* (T(ty.Nmax) -Tamb);
                % Heat flow rates between HX1..4 and storage tank
                if ~isempty(mm1)
                    for k = 1:length(mm1)
                        Tdot(4.*ty.Nmax+mm1(k)) = UA_h14(mm1(k)) ./ length(mm1) .* (T(mm1(k)) - T(ty.Nmax+mm1(k)));
                    end
                end
                if ~isempty(mm4)
                    for k = 1:length(mm4)
                        Tdot(4.*ty.Nmax+mm4(k)) = UA_h14(mm4(k)) ./ length(mm4) .* (T(mm4(k)) - T(ty.Nmax+mm4(k)));
                    end
                end
                if ~isempty(mm2)
                    for k = 1:length(mm2)
                        Tdot(5.*ty.Nmax+mm2(k)) = UA_h23(mm2(k)) ./ length(mm2) .* (T(mm2(k)) - T(2.*ty.Nmax+mm2(k)));
                    end
                end
                if ~isempty(mm3)
                    for k = 1:length(mm3)
                        Tdot(5.*ty.Nmax+mm3(k)) = UA_h23(mm3(k)) ./ length(mm3) .* (T(mm3(k)) - T(2.*ty.Nmax+mm3(k)));
                    end
                end
                % Heat losses tank floor/-top to ambient
                Tdot(6.*ty.Nmax+1) = -ty.UA_u .* (T(1) - Tamb);
                Tdot(6.*ty.Nmax+2) = -ty.UA_o .* (T(ty.Nmax) - Tamb);
                % Consumption auxiliary heater
                Tdot(6 .* ty.Nmax+3) = sum(ty.ksiak) .* Paux .* (T(ty.auxpos) <= ty.aux.T) + sum(ty.ksiav) .* Paux;
            end
        end
    end
end

