%% Test simulation of the mtype340: charging/discharging scenario (stratified) with a double port each
%--> constant power




%% fixed inputs:
Vs = 1.25; % [m�] -> 1250 litre tank
Hs = 1.7; % [m] height of storage tank
UA_a = 4.2691; % heat loss rates of the zones 1..4 [W/K]
UA_u = 0.333; % heat loss rate bottom [W/K]
UA_o = 0.333; % heat loss rate top [W/K]
% Heat exchangers not used - HX variables are set to nans:
zhi = -ones(4,1); zho = zhi; Vh = nan(4,1); UA_h = zhi; Thi = zhi; mdoth = zhi;
zs = linspace(0.05,0.95,10); % 10 temperature sensors (relative positions)
aux.var = 2; Paux = nan; % no auxiliary heater
delt = 5; %60;%300;%300; % time step size [s]
nsth = 3600./delt; % Number of passes per hour
Tamb = 15; % ambient temperature [�C]
scdp = 1; % stratified charging
Nmax = 10; % number of nodes

%% Variable parameters (change in the course of the simulation)
% double ports
zdi = 0.05; % Charging from below (relative position)
zdo = 0.04; % rel. position of the DP output
% Tdi = 90; %DP-input temperature [�C]
Tz = 15.*ones(Nmax,1); % Start: 15 �C in the whole storage tank

%% Initialize mtype340 object
ty = mtype340(Hs, Vs, UA_a, UA_u, UA_o, UA_h, zdi, zdo, scdp, aux, zhi, zho, Vh, zs, delt);


% Init temperature-matrix
Tmat = nan(Nmax,14.*nsth-1);
Tmat(:,1) = Tz;
SOC = nan(1,14.*nsth-1);
SOC(1) = (mean(Tz)-15)./75.*100;
Qout = zeros(1,14.*nsth-1);
eta = Qout; eta(1) = 1;
mdots = zeros(1,15.*nsth-1); %Massenstr�me
mmdots = zeros(1,6);%mittlerer Massenstrom bei jedem Abschnitt
Qloss = zeros(size(SOC));
Qinput = Qloss;


%% [1] Heating to stand-by temperature, 30 kW, 1 hour
%Assumption DP-output. Tdo = 15 �C
%with cpw = 4179,5 Ws/(kg*K)
cpw = 4179.5;
Qin = 30000; %Input-power [W]

% Set properties
ty.Tdi = 90; % DP-input temperature [�C]
ty.mdotd = -Qin./(cpw.*(ty.Tdi-15)); % mass flow rate [kg/s]
ty.Tamb = 15; % ambient temperature [�C]

% Initialize simulation
ty.simulate;

% Save values
Qinput(2) = Qin;
Tmat(:,2) = ty.Ts;
SOC(2) = (mean(ty.Tz)-15)./75.*100;
Qout(2) = ty.Qd;
eta(2) = ty.Qd./Qin;
Qloss(2) = ty.Qls;
t = ty.delt./3600; % simulated time in hours
% Loop for further simulations in section [1]
for i = 2:nsth-1 %#ok<*BDSCI>
    % Determine new mass flow rate for constant power
    ty.mdotd = -Qin./(cpw.*(ty.Tdi-ty.Tdo));
    mdots(i) = ty.mdotd;
    ty.simulate;
    % Save properties in arrays
    Tmat(:,i) = ty.Ts;
    SOC(i) = (mean(ty.Tz)-15)./75.*100;
    t = [t; t(end) + ty.delt./3600]; %#ok<*AGROW>
    Qout(i) = ty.Qd;
    eta(i) = ty.Qd./Qin;
    Qloss(i) = ty.Qls;
    Qinput(i) = Qin;
end
mmdots(1) = mean(mdots(1:i));
%{
fig = figure;
plot(t,Tmat(:,1:nsth-1)'); hold on
%}

%% [2] down time, 1 hour
Qin = 0;
ty.mdotd = 0;
for i = nsth:2*nsth-1
    ty.simulate;
    % Save properties in array
    Tmat(:,i) = ty.Ts;
    t = [t; t(end) + delt./3600];
    SOC(i) = (mean(ty.Tz)-15)./75.*100;
    Qout(i) = ty.Qd;
    eta(i) = ty.Qd./Qin;
    Qloss(i) = ty.Qls;
    Qinput(i) = Qin;
end
%{
fig = figure;
plot(t,Tmat(:,1:2*nsth-1)'); hold on
%}

%% [3] Charge with 20 kW, 4 hours
Qin = 20000;
for i = 2*nsth:6*nsth-1
    % Calculate mass flow rate required for constant power
    ty.mdotd = -Qin./(cpw.*(ty.Tdi-ty.Tdo));
    if isnan(ty.mdotd) || ty.mdotd == -inf % correction
        ty.mdotd = mdots(i-1);
    end
    mdots(i) = ty.mdotd;
    ty.simulate;
    % Save properties
    Tmat(:,i) = ty.Ts;
    t = [t; t(end) + ty.delt./3600];
    SOC(i) = (mean(ty.Tz)-15)./75.*100;
    Qout(i) = ty.Qd;
    eta(i) = ty.Qd./Qin;
    Qloss(i) = ty.Qls;
    Qinput(i) = Qin;
end
mmdots(3) = mean(mdots(2*nsth:i));
%{
fig = figure;
plot(t,Tmat(:,1:6*nsth-1)'); hold on
%}


%% [4] down time, 1 hour
Qin = 0;
ty.mdotd = 0;
for i = 6*nsth:7*nsth-1
    ty.simulate;
    Tmat(:,i) = ty.Ts;
    t = [t; t(end) + ty.delt./3600];
    SOC(i) = (mean(ty.Tz)-15)./75.*100;
    Qout(i) = ty.Qd;
    eta(i) = ty.Qd./Qin;
    Qloss(i) = ty.Qls;
    Qinput(i) = Qin;
end
%{
fig = figure;
plot(t,Tmat(:,1:7*nsth-1)'); hold on
%}

%% [5] Discharging with 20 kW, 3.5 hours
Qin = -20000;
ty.Tdi = 15;
ty.zdo = 0.95;
ty.aux.var = 0; % constant auxiliary heater AUX power according to controller
ty.aux.pos = 0.55; % relative height of AUX (6th node)
ty.aux.T = 60; % Minimum temperature of 60 �C at respective position
ty.aux.var = 2;
ty.Paux = 0;
% Backup heating with heat exchangers
%ty.zhi(1) = 0.75;
%ty.zho(1) = 0.05;
%ty.Vh(1) = 0.01;
%ty.Thi(1) = 90;
%ty.UA_h(1) = 542;
%ty.mdoth(1) = ty.Paux./(cpw.*(ty.Thi(1)-ty.Tz(1)));

for i = 7*nsth:10.5*nsth-1
    ty.mdotd = Qin./(cpw.*(ty.Tdi-ty.Tdo));
    mdots(i) = ty.mdotd;
    ty.simulate;
    Tmat(:,i) = ty.Ts;
    t = [t; t(end) + ty.delt./3600];
    SOC(i) = (mean(ty.Tz)-15)./75.*100;
    Qloss(i) = ty.Qls;
    Qinput(i) = Qin+ty.Qaux;
end
mmdots(5) = mean(mdots(7*nsth:i));
%{
fig = figure;
plot(t(1:10.5*nsth-1),Tmat(:,1:10.5*nsth-1)'); hold on
%}

%% [6] Discharging with 15 kW and securing the stand-by temperature, 3.5 hours
Qin = -15000;
ty.zdo = 0.75; % Withdrawal of water - 60�C
%ty.zdo = 0.95; % water withdrawal (as hot as possible)

% % Back-up heating with auxiliary heater
% ty.aux.var = 0;
% %ty.aux.pos = 0.65; % AUX in 7th node
% %ty.Paux = 0;
% ty.Paux = 15000;
% ty.aux.T = 60; % AUX control temperature: 60 �C

% Back-up heating with heat exchangers
ty.zhi(1) = 0.75;
ty.zho(1) = 0.01;
ty.Vh(1) = 0.01;
ty.Thi(1) = 90;
ty.UA_h(1) = 542;
ty.Tho(1) = Tz(2); 

% Back-up heating with double ports
% ty.zdi = [ty.zdi; ty.zdi];
% ty.zdo = [ty.zdo; 0.04];
% ty.Tdi = [15; 80]; ty.scdp = [1; 1];
% ty.Tdo = [ty.Tdo; ty.Tz(1)];
% Qhs = 0; Qh = 0;
Qhxin = [];
auxon = true; % Auxiliary heater on
Qhxtran = []; thout = [];

for i = 10.5*nsth:14*nsth-1
    Qhi = -(Qin + ty.Qltop + ty.Qltop + ty.Qls)+800; % backup with HX
    ty.mdoth(1) = -Qhi./(cpw.*(ty.Thi(1)-ty.Tho(1))); % backup with HX
    if isnan(ty.mdoth(1))
       ty.mdoth(1) = lastmdoth;
    end
    ty.mdotd = Qin./(cpw.*(ty.Tdi-ty.Tdo));
    if isnan(ty.mdotd(1)) || ty.mdotd(1) == -inf
        ty.mdotd(1) = lastmdotd;
    end
    ty.simulate;
    lastmdotd = ty.mdotd(1);
    lastmdoth = ty.mdoth(1);
    Tmat(:,i) = ty.Ts;
    t = [t; t(end) + ty.delt./3600];
    SOC(i) = (mean(ty.Tz)-15)./75.*100;
    Qout(i) = ty.Qd;
    eta(i) = ty.Qd./(Qin+ty.Paux.*auxon);
    Qinput(i) = Qin+ty.Qh(1);
    Qloss(i) = ty.Qls;
    Qhxin = [Qhxin; ty.Qh(1)];
    Qhxtran = [Qhxtran; ty.Qhs(1)];
    thout = [thout; ty.Tho(1)];
end
mmdots(6) = mean(mdots(10.5*nsth:i));

%% Plot results
fig = figure;
fig.Position(3) = 2.*fig.Position(3);
plot(t,flipud(Tmat(:,1:14*nsth-1))','LineWidth',2); hold on
plot(t,SOC,'k--','LineWidth',2);
legend('\vartheta_{top}','\vartheta_9','\vartheta_8','\vartheta_7',...
'\vartheta_6','\vartheta_5','\vartheta_4','\vartheta_3',...
'\vartheta_2','\vartheta_{bottom}','SOC');
lh=findall(gcf,'tag','legend');
set(lh,'location','northeastoutside');
xlabel('time in h')
ylabel('Temperature in �C and SOC in %')
hold on
y = get(gca,'ylim');
plot([1,1],y,'k'); hold on
plot([2,2],y,'k'); hold on
plot([6,6],y,'k'); hold on
plot([7,7],y,'k');
plot([10.5,10.5],y,'k')
text(0.45,95,'[1]')
text(1.45,95,'[2]')
text(3.95,95,'[3]')
text(6.45,95,'[4]')
text(8.75,95,'[5]')
text(12.4,95,'[6]')
title('Temperatures in the storage tank nodes')
