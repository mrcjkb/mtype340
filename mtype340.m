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
        end
        
        function ty = simulate(ty, Tdi, mdotd, Thi, mdoth, Tamb, Paux)
            % SIMULATE: Pass simulation variables to mtype340 object and update
            
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
            ksi1 = ones(size(mdotd)); ksi1(mdotd<=0) = 0;
            % ksi2 = 1 for mdot_dp < 0, otherwise ksi2 = 0
            ksi2 = ones(size(mdotd)); ksi2(mdotd>=0) = 0;
            % ksi3 = 1 für zones i, that are in contact with HX 1 oder 4 sind,
            % otherwise ksi3 = 0
            if ty.ini
                ty.hxdat.nh14 = [];
                ty.hxdat.nh23 = [];
            end
            if ty.ini || ty.hxdat.nh14 ~= length(find([ty.zho(1),ty.zho(4)]>=0)) || ty.hxdat.nh23 ~= length(find([ty.zho(2),ty.zho(3)]>=0))
                % Find out if hx1 & hx4 are used alone, together or not at all
                numhx14 = length(find([ty.zho(1),ty.zho(4)]>=0));
                % Determine parallel temperature layers of HX 1 & 4
                Th14 = zeros(size(ty.Tz));
                switch numhx14
                    case 0 % neither are used
                        mm1 = []; mm4 = []; %Indizes für die entsprechenden Nodes
                        ksi3 = zeros(1,Nmax);
                        nh14 = 1; %wird 1 gesetzt, um Teilen durch Null zu vermeiden
                        Th14 = Tz;
                        h1nodes = []; h4nodes = [];
                    case 1
                        ksi3 = zeros(1,Nmax);
                        if zho(1) >= 0 %falls nur WT1 verwendet wird
                            mm4 = [];
                            set = [find(irel_bndr > zho(1),1)-1; find(irel_bndr > zhi(1),1)-1];
                            %Unterscheidung zw. B- und Entladen des Speichers
                            if zho(1) > zhi(1) %Entladen
                                mm1 = set(2):set(1); l = length(mm1);
                                h1nodes = linspace(Thi(1),Tz(set(1))-5,l+2);
                            else %Beladen
                                mm1 = set(1):set(2); l = length(mm1);
                                h1nodes = linspace(Tz(set(1))+5,Thi(1),l+2);
                            end
                            %Einsetzen der Temperatur des WT1 (zun. idealisierte Annahme des Verlaufs)
                            Th14(mm1) = h1nodes(2:end-1);
                            Th14(mm1) = repmat(Thi(1),size(mm1)); %Initialisierung mit konstanter hx-Temperatur
                            ksi3(mm1) = 1;
                            h4nodes = [];
                        else %falls nur WT4 verwendet wird
                            mm1 = [];
                            set = [find(irel_bndr > zho(4),1)-1; find(irel_bndr > zhi(4),1)-1];
                            if zho(4) > zhi(4) %Entladen
                                mm4 = set(2):set(1); l = length(mm4);
                                h4nodes = linspace(Thi(4),Tz(set(1))-5,l+2);
                            else %Beladen
                                mm4 = set(1):set(2); l = length(mm4);
                                h4nodes = linspace(Tz(set(1))+5,Thi(4),l+2);
                            end
                            Th14(mm4) = h4nodes(2:end-1);
                            Th14(mm4) = repmat(Thi(4),size(mm4));
                            ksi3(mm4) = 1;
                            h1nodes = [];
                        end
                        nh14 = length(find(ksi3)); %Anzahl Nodes, die von WT1 & WT4 besetzt sind;
                        %ksi3 = repmat(ksi3,length(time),1);
                    case 2 %falls WT1 und WT4 beide verwendet werden
                        set = [find(irel_bndr > zho(1),1)-1; find(irel_bndr > zhi(1),1)-1;...
                            find(irel_bndr > zho(4),1)-1; find(irel_bndr > zhi(4),1)-1];
                        %Einsetzen der Temperaturen der WT 1 & 4 entsprechend ihrer Nodes (zunächst idealisierte Annahme des Verlaufs)
                        %Unterscheidung zwischen Be- und Entladen des Speichers
                        if zho(1) > zhi(1) %Entladen WT1
                            mm1 = set(2):set(1); l = length(mm1);
                            h1nodes = linspace(Thi(1),Tz(set(1))-5,l+2);
                        else %Beladen WT1
                            mm1 = set(1):set(2); l = length(mm1);
                            h1nodes = linspace(Tz(set(1))+5,Thi(1),l+2);
                        end
                        Th14(mm1) = h1nodes(2:end-1);
                        Th14(mm1) = repmat(Thi(1),size(mm1));
                        if zho(4) > zhi(4) %Entladen WT4
                            mm4 = set(4):set(3); l = length(mm4);
                            h4nodes = linspace(Thi(4),Tz(set(3))-5,l+2);
                        else %Beladen WT4
                            mm4 = set(3):set(4); l = length(mm4);
                            h4nodes = linspace(Tz(set(3))+5,Thi(4),l+2);
                        end
                        Th14(mm4) =  h4nodes(2:end-1);
                        Th14(mm4) = repmat(Thi(4),size(mm4));
                        set = [mm1,mm4];
                        ksi3 = zeros(1,Nmax); ksi3(set) = 1;
                        nh14 = length(find(ksi3)); %Anzahl Nodes, die von WT1 & WT4 besetzt sind;
                        %ksi3 = repmat(ksi3,length(time),1);
                end
                %ksi4 = 1 für Zonen i, die in Kontakt mit Wärmetauscher 2 oder 3 sind,
                %ansonsten ksi4 = 0
                %Ermittlung, ob hx1, hx4, alleine, zusammen oder gar nicht genutzt werden:
                numhx23 = length(find([zho(2),zho(3)]>=0));
                Th23 = zeros(size(Tz)); %Initialisierung der entspr. Temperaturschichten der WT 2 & 3
                switch numhx23
                    case 0 %falls WT2 und WT3 beide nicht verwendet werden
                        mm2 = []; mm3 = [];
                        ksi4 = zeros(1,Nmax);
                        nh23 = 1; %wird 1 gesetzt, um Teilen durch Null zu vermeiden
                        Th23 = Tz;
                        h2nodes = []; h3nodes = [];
                    case 1
                        ksi4 = zeros(1,Nmax);
                        if zho(2) >= 0 %falls nur WT2 verwendet wird
                            mm3 = [];
                            set = [find(irel_bndr > zho(2),1)-1; find(irel_bndr > zhi(2),1)-1];
                            if zho(2) > zhi(2) %Entladen
                                mm2 = set(2):set(1); l = length(mm2);
                                h2nodes = linspace(Thi(2),Tz(set(1))-5,l+2);
                            else %Beladen
                                mm2 = set(1):set(2); l = length(mm2);
                                h2nodes = linspace(Tz(set(1))+5,Thi(2),l+2);
                            end
                            %Einsetzen der Temperatur des WT2 (zun. idealisierte Annahme des Verlaufs)
                            Th23(mm2) = h2nodes(2:end-1);
                            Th23(mm2) = repmat(Thi(2),size(mm2));
                            ksi4(mm2) = 1;
                            h3nodes = [];
                        else %falls nur WT3 verwendet wird
                            mm2 = [];
                            set = [find(irel_bndr > zho(3),1)-1; find(irel_bndr > zhi(3),1)-1];
                            if zho(3) > zhi(3) %Entladen
                                mm3 = set(2):set(1); l = length(mm3);
                                h3nodes = linspace(Thi(3),Tz(set(1))-5,l+2);
                            else %Beladen
                                mm3 = set(1):set(2); l = length(mm3);
                                h3nodes = linspace(Tz(set(1))+5,Thi(3),l+2);
                            end
                            %Einsetzen der Temperatur des WT3 (zun. idealisierte Annahme des Verlaufs)
                            Th23(mm3) = h3nodes(2:end-1);
                            Th23(mm3) = repmat(Thi(3),size(mm3));
                            ksi4(mm3) = 1;
                            h2nodes = [];
                        end
                        nh23 = length(find(ksi4)); %Anzahl Nodes, die von WT2 & WT3 besetzt sind
                        %ksi4 = repmat(ksi4,length(time),1);
                    case 2 %falls WT2 und WT3 beide verwendet werden
                        set = [find(irel_bndr > zho(2),1)-1; find(irel_bndr > zhi(2),1)-1;...
                            find(irel_bndr > zho(3),1)-1; find(irel_bndr > zhi(3),1)-1];
                        if zho(2) > zhi(2) %Entladen WT2
                            mm2 = set(2):set(1); l = length(mm2);
                            h2nodes = linspace(Thi(2),Tz(set(1))-5,l+2);
                        else %Beladen WT2
                            mm2 = set(1):set(2); l = length(mm2);
                            h2nodes = linspace(Tz(set(1))+5,Thi(2),l+2);
                        end
                        Th23(mm2) = h2nodes(2:end-1);
                        Th23(mm2) = repmat(Thi(2),size(mm2));
                        if zho(3) > zhi(3) %Entladen WT3
                            mm3 = set(4):set(3); l = length(mm3);
                            h3nodes = linspace(Thi(3),Tz(set(3))-5,l+2);
                        else %Beladen WT3
                            mm3 = set(3):set(4); l = length(mm3);
                            h3nodes = linspace(Tz(set(3))+5,Thi(3),l+2);
                        end
                        Th23(mm3) =  h3nodes(2:end-1);
                        Th23(mm3) = repmat(Thi(3),size(mm3));
                        set = [mm2,mm3];
                        ksi4 = zeros(1,Nmax); ksi4(set) = 1;
                        nh23 = length(find(ksi4));
                        
                end
                ksi3 = ksi3'; ksi4 = ksi4';
                hxdat.ksi3 = ksi3; hxdat.ksi4 = ksi4;
                hxdat.mm1 = mm1; hxdat.mm2 = mm2; hxdat.mm3 = mm3; hxdat.mm4 = mm4;
                hxdat.nh14 = nh14; hxdat.nh23 = nh23;
                hxdat.h1nodes = h1nodes; hxdat.h2nodes = h2nodes; hxdat.h3nodes = h3nodes; hxdat.h4nodes = h4nodes;
            else
                ksi3 = hxdat.ksi3; ksi4 = hxdat.ksi4;
                mm1 = hxdat.mm1; mm2 = hxdat.mm2; mm3 = hxdat.mm3; mm4 = hxdat.mm4;
                nh14 = hxdat.nh14; nh23 = hxdat.nh23;
                h1nodes = hxdat.h1nodes; h2nodes = hxdat.h2nodes; h3nodes = hxdat.h3nodes; h4nodes = hxdat.h4nodes;
            end
            
            
            %% Feststellen eines Starts der Wärmetauscher
            %Die Massenträgheit beim Start eines WT wird im Solver berücksichtigt, da
            %diese von der Anzahl interner Zeitschritte abhängt.
            %In diesem Abschnitt wird ein Anfahren eines der 4 WT erkannt
            if ini
                hxdat.mdoth = mdoth;
                hxdat.mdoth(isnan(hxdat.mdoth)) = 0;
                hxini = ones(4,1);
            else
                hxini = zeros(4,1);
                for i = 1:4
                    if hxdat.mdoth(i) == 0 && mdoth(i) ~= 0
                        hxini(i) = 1;
                    end
                end
                hxdat.mdoth = mdoth;
                hxdat.mdoth(isnan(hxdat.mdoth)) = 0;
            end
            
            
            %____Ermittlung des Vektors ndzk__________________________________________________________
            %ndzk: Vektor mit Anzahl Nodes in der zum jew. Node gehörigen Zone k (k = 1...4)
            %{
ndz = Nmax*dz;
%Ermittlung der Anzahl Zonen konstanter Wärmeverlustrate
for i = 1:4
    if sum(dz(1:i)) == 1
        numz = i;
        break
    end
end
switch numz %problem: Generates 12x1 matrix with example data
    case 1 %Fall 1: Es gibt nur 1 Zone --> Zone z1 besetzt alle Nmax Nodes
        ndzk = repmat(Nmax,1,Nmax)'; %Anzahl zonen für die WVR auf Speicherschichten skaliert
        UAs_ak = repmat(UA_a(1),1,Nmax)'; %WVR der entspr. Zonen auf Speicherschichten skaliert
    case 2 %Fall 2: Es gibt 2 Zonen --> Nodes werden zw. z1 & z2 aufgeteilt
        ndzk = [repmat(ndz(1),1,uint8(ndz(1))), repmat(ndz(2),1,uint8(ndz(2)))]';
        UAs_ak = [repmat(UA_a(1),1,uint8(ndz(1))); repmat(UA_a(2),1,uint8(ndz(2)))]';
    case 3 %Fall 3: Es gibt 3 Zonen
        ndzk = [repmat(ndz(1),1,uint8(ndz(1))), repmat(ndz(2),1,uint8(ndz(2))), repmat(ndz(3),1,uint8(ndz(3)))]';
        UAs_ak = [repmat(UA_a(1),1,uint8(ndz(1))), repmat(UA_a(2),1,uint8(ndz(2))), repmat(UA_a(3),1,uint8(ndz(3)))]';
    case 4 %Fall 4: Es gibt 4 Zonen
        ndzk = [repmat(ndz(1),1,uint8(ndz(1))), repmat(ndz(2),1,uint8(ndz(2))), repmat(ndz(3),1,uint8(ndz(3))), repmat(ndz(4),1,uint8(ndz(4)))]';
        UAs_ak = [repmat(UA_a(1),1,uint8(ndz(1))), repmat(UA_a(2),1,uint8(ndz(2))), repmat(UA_a(3),1,uint8(ndz(3))), repmat(UA_a(4),1,uint8(ndz(4)))]';
end
            %}
            
            
            %% iteration of Gl. 4.1 & Gl. 4.2 & auxiliary heater for each node
            
            ksi5 = zeros(4,1);
            ksi5(mdoth > 0) = 1;
            ksi6 = ~ksi5; %ksi5 = 0 für negtive WT-Massenströme, ksi6 = 0 für positive WT-Massenströme
            ksi5(mdoth == 0) = 0;
            ksi6(mdoth == 0) = 0;
            
            nodevect = zeros(Nmax,1); nodevect([mm1, mm4]) = 1;
            rest = find(nodevect==0);
            Th14(rest) = Tz(rest); %Gleichsetzen mit Speichertemperatur, damit deltaT = 0
            nodevect = zeros(Nmax,1); nodevect([mm2, mm3]) = 1;
            rest = find(nodevect==0);
            Th23(rest) = Tz(rest);
            
            if ini
                Qlsx = zeros(Nmax,1); %Verluste Speicher-Umgebung
                Qh14_s = zeros(Nmax,1); %Wärmefluss HX1/4-Speicher
                Qh23_s = Qh14_s; %Wärmefluss HX2/3-Speicher
                Qlbot = 0; %Verluste Speicherboden
                Qltop = 0; %Verluste Speicherdeckel
            else
                Qlsx = storedat.Qlsx;
                Qh14_s = storedat.Qh14_s;
                Qh23_s = storedat.Qh23_s;
                Qlbot = storedat.Qlbot;
                Qltop = storedat.Qltop;
                Th14 = storedat.Th14;
                Th23 = storedat.Th23;
            end
            %mdotd = abs(mdotd);
            
            %% Solver-Anwendung
            allTs = [Tz; Th14; Th23; Qlsx; Qh14_s; Qh23_s; Qlbot; Qltop; Paux];
            
            options = odeset('AbsTol',0.1); %Ändern der absoluten Fehlertoleranz auf 0,1
            [~,dT] = solver(@tempsolver,[0 delt], allTs, options);
            
            
            %% % für AUTOMATISCHE SOLVER-AUSWAHL diesen Abschnitt aktivieren
            %{

if sum(isnan(mdoth)) == 0 && aux.var ~= 2 %es werden nur Doppelports (oder Heizstab variabler Leistung) verwendet
    [t,dT] = ode23(@tempsolver,[0 delt], allTs);
elseif aux.var ~= 2 %Es werden Doppelports und Wärmetauscher verwendet
    if delt <= 650
        [t,dT] = ode23(@tempsolver,[0 delt], allTs);
    else
        [t,dT] = ode23tb(@tempsolver,[0 delt], allTs);
    end
elseif sum(isnan(mdoth)) == 0 && aux.var == 2 %Es werden Doppelports & Heizstab (konstanter Leistung) verwendet
    [t,dT] = ode45(@tempsolver,[0 delt], allTs);
else %Es werden DP, WT und Heizstab konstanter Leistung verwendet
    if delt < 120
        [t,dT] = ode45(@tempsolver,[0 delt], allTs);
    elseif delt >= 120 && delt < 170
        [t,dT] = ode113(@tempsolver,[0 delt], allTs);
    else
        [t,dT] = ode23(@tempsolver,[0 delt], allTs);
    end
end
            %}
            %%
            Tz = dT(end,1:Nmax)';
            Th14 = dT(end,Nmax+1:2*Nmax)';
            storedat.Th14 = Th14;
            Th23 = dT(end,2*Nmax+1:3*Nmax)';
            storedat.Th23 = Th23;
            
            %Wärmeströme durch die Wärmetauscher
            Qh1 = zeros(Nmax,1); Qh2 = Qh1; Qh3 = Qh1; Qh4 = Qh1;
            Qh1(mm1) = (Vh(1).*cp_h(1).*rho_h)./length(mm1).*Th14(mm1)./delt;
            Qh2(mm2) = (Vh(2).*cp_h(2).*rho_h)./length(mm2).*Th23(mm2)./delt;
            Qh3(mm3) = (Vh(3).*cp_h(3).*rho_h)./length(mm3).*Th23(mm3)./delt;
            Qh4(mm4) = (Vh(4).*cp_h(4).*rho_h)./length(mm4).*Th14(mm4)./delt;
            %Qh1 = Qh1(Qh1~=0); Qh2 = Qh2(Qh2~=0); Qh3 = Qh3(Qh3~=0); Qh4 = Qh4(Qh4~=0);
            Qh = [sum(Qh1); sum(Qh2); sum(Qh3); sum(Qh4)];
            Qh(Qh == 0) = nan;
            
            %Qlsx = mean(dT(:,3*Nmax+1:4*Nmax))./delt; Qlsx = Qlsx';
            Qlsx = dT(end,3*Nmax+1:4*Nmax)'./delt;
            storedat.Qlsx = Qlsx;
            %Qh14_s = mean(dT(:,4*Nmax+1:5*Nmax))./delt; Qh14_s = Qh14_s';
            Qh14_s = dT(end,4*Nmax+1:5*Nmax)'./delt;
            storedat.Qh14_s = Qh14_s;
            %Qh23_s = mean(dT(:,5*Nmax+1:6*Nmax))./delt; Qh23_s = Qh23_s';
            Qh23_s = dT(end,5*Nmax+1:6*Nmax)'./delt;
            storedat.Qh23_s = Qh23_s;
            %Qlbot = mean(dT(:,6*Nmax+1))./delt;
            Qlbot = dT(end,6*Nmax+1)./delt;
            
            storedat.Qlbot = Qlbot;
            %Qltop = mean(dT(:,6*Nmax+2))./delt;
            Qltop = dT(end,6*Nmax+2)./delt;
            storedat.Qltop = Qltop;
            %Qaux = sum(dT(:,6*Nmax+3))./delt;
            %stag = Qaux./Paux.*length(t);
            %Qaux = [t', dT(:,6*Nmax+3)'];
            
            %Qaux = sum((dT(:,auxpos) >= aux.T).*[0; diff(t)])./sum((dT(:,auxpos) <= aux.T).*[0; diff(t)]).*Paux;
            %Qaux = mean(dT(:,6*Nmax+3))./delt;
            Qaux = dT(end,6*Nmax+3)./delt;
            
            %Qh14 = (dT(end,Nmax+1:2*Nmax)'./(Vh(1).*rho_h.*cp_h))./delt + abs(storedat.Qh14_s);
            
            
            %% Ermittlung der Ausgangstemperaturen der Doppelports
            Tdo = nan(size(Tdi));
            if ~isnan(Tdi(1))
                set = Tdo;
                for i = 1:length(set)
                    set(i) = find(irel_bndr > zdo(i),1)-1;
                end
                Tdo = Tz(set);
            end
            %% Ausgabe der Wärmeströme
            Qls = sum(Qlsx)+Qlbot+Qltop; %ges. Wärmeverluste des Speichers
            if ~isnan(Tdi(1))
                Qd = abs(mdotd).*cp_s.*(Tdi-Tdo); %Wärmeströme durch die Doppelports
            else
                Qd = nan;
            end
            %Qh = abs(mdoth).*cp_h.*(Thi-Tho); %Wärmeströme durch die Wärmetauscher
            %Wärmetauscher 1...4 - Speicher
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
            Qhs = [Qh1s; Qh2s; Qh3s; Qh4s];
            %QhDiff = Qh + Qhs; %Restwärme nach der Speicherbeladung
            %QhDiff
            %Qhs
            %Tho = Thi-Qh./(cp_h.*abs(mdoth));
            %Q_WT = cp*mdot*(Tein-Taus)
            
            %% Extraktion der Ausgangstemperaturen der Wärmetauscher
            Tho = nan(4,1);
            Thm = Tho;
            if ~isempty(mm1)
                h1nodes(2:end-1) = Th14(mm1);
                chk = h1nodes ~= Thi(1);
                if chk(1) == 0
                    %Tho(1) = h1nodes(end-1);
                    Tho(1) = Thi(1)-abs(Qhs(1)./(mdoth(1).*cp_h(1)));
                    h1nodes(end) = Tho(1);
                else
                    % Tho(1) = h1nodes(2);
                    Tho(1) = Thi(1)-abs(Qhs(1)./(mdoth(1).*cp_h(1)));
                    h1nodes(1) = Tho(1);
                end
                Thm(1) = mean(h1nodes);
                hxdat.h1nodes = h1nodes;
                %cp_h(1) = storedat.cpw(round(Thm(1)));
            end
            if ~isempty(mm2)
                h2nodes(2:end-1) = Th23(mm2);
                chk = h2nodes ~= Thi(2);
                if chk(1) == 0
                    %Tho(2) = h2nodes(end-1);
                    Tho(2) = Thi(2)-abs(Qhs(2)./(mdoth(2).*cp_h(2)));
                    h2nodes(end) = Tho(2);
                else
                    %Tho(2) = h2nodes(2);
                    Tho(2) = Thi(2)-abs(Qhs(2)./(mdoth(2).*cp_h(2)));
                    h2nodes(1) = Tho(2);
                end
                Thm(2) = mean(h2nodes);
                hxdat.h2nodes = h2nodes;
                %cp_h(2) = storedat.cpw(round(Thm(2)));
            end
            if ~isempty(mm3)
                h3nodes(2:end-1) = Th23(mm3);
                chk = h3nodes ~= Thi(3);
                if chk(1) == 0
                    %Tho(3) = h3nodes(end-1);
                    Tho(3) = Thi(3)-abs(Qhs(3)./(mdoth(3).*cp_h(3)));
                    h3nodes(end) = Tho(3);
                else
                    %Tho(3) = h3nodes(2);
                    Tho(3) = Thi(3)-abs(Qhs(3)./(mdoth(3).*cp_h(3)));
                    h3nodes(1) = Tho(3);
                end
                Thm(3) = mean(h3nodes);
                hxdat.h3nodes = h3nodes;
                %cp_h(3) = storedat.cpw(round(Thm(3)));
            end
            if ~isempty(mm4)
                h4nodes(2:end-1) = Th14(mm4);
                chk = h4nodes ~= Thi(4);
                if chk(1) == 0
                    %Tho(4) = h4nodes(end-1);
                    Tho(4) = Thi(4)-abs(Qhs(4)./(mdoth(4).*cp_h(4)));
                    h4nodes(end) = Tho(4);
                else
                    %Tho(4) = h4nodes(2);
                    Tho(4) = Thi(4)-abs(Qhs(4)./(mdoth(4).*cp_h(4)));
                    h4nodes(1) = Tho(4);
                end
                Thm(4) = mean(h4nodes);
                hxdat.h4nodes = h4nodes;
                %cp_h(4) = storedat.cpw(round(Thm(4)));
            end
            %hxdat.cp_h = cp_h;
            
            
            
            
            
            
            %% Temperatursensoren-Ausgabe + Temperatur an der Stelle des Heizstabes
            Ts = nan(size(zs));
            if ~isnan(zs(1))
                set = Ts;
                for i = 1:length(zs)
                    set(i) = find(irel_bndr > zs(i),1)-1;
                end
                Ts = Tz(set);
            end
            if aux.var ~= 2
                Taux = Tz(auxpos);
            else
                Taux = nan;
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

