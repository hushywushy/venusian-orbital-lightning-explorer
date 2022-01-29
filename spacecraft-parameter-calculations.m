format compact;
%% Inputs
%Venus
R_v = 6034*1000;%Venus radius, units = m 
M_v = 4.87*10^(24);  %Venus mass, units = kg
G = 6.67408*10^(-11); %gravitational constant %units = m3 kg-1 s-2
g = 8.87;  %gravity, units = m/s2
rho = 9.81*10^(-10); %atmospheric density

% Mothership
apogee_ms = 60000*1000; %Mothership apogee altitude, units = m
    a_ms = apogee_ms + R_v;
perigee_ms = 500*1000; %Mothership perigee altitude, units = m
    p_ms = perigee_ms + R_v;
semimajor_ms = apogee_ms + perigee_ms + R_v; %Mothership Semi-Major Axis, units = m
    e_ms = (semimajor_ms/2-p_ms)/(semimajor_ms/2); 
    
%CubeSat
m_s = 9; %CubeSat mass, units = kg
A_ram = (10*10)*0.0001; %CubeSat area, units = m2
C_d = 2.2; %CubeSat coefficient of drag 
apogee_s = 60000*1000; %CubeSat apogee altitude, units = m
	a_s = apogee_s + R_v;
perigee_s = 500*1000; %CubeSat perigee altitude, units = m
    p_s = perigee_s + R_v;
alt_initial = apogee_s + perigee_s; %CubeSat initial "circular" altitude, units = m
    a_initial = alt_initial + R_v; %CubeSat initial "circular" orbital radius, units = m
        h_initial = alt_initial/1000; %CubeSat initial "circular" orbital radius, units = km
alt_final = 500*1000; %CubeSat final "circular" altitude, units = m
    a_final = alt_final + R_v; %CubeSat final "circular" orbital radius, units = m
        h_final = alt_final/1000; %CubeSat final "circular" orbital radius, units = km
f10_a = 300;
Ap_a = 0;
t = 0.5; %units = year
P_thermal = 14.4; %units = W
P_camera = 240*10^(-3); %units = W
P_radio = 3; %includes radio % units = W
P_OBC = 3*550*10^(-3); %units = W
P_ACS = 2; %units = W
P_comm = 4.3; %units = W
KE = 10.8; %Energy of CubeSat release, units = J
alt_maneuver = 500*1000; %altitude of maneuver, units = m 
max_altitude = 3200; %units = km
g2 = g*(R_v/(R_v+alt_maneuver))^2;  %gravity at altitude of maneuver, units = m/s2
    
%% Mothership Orbital Calculations

v_ms_perigee = sqrt(G*M_v*(1+e_ms)/p_ms); %Mothership Perigee velocity, units = m/s
v_ms_apogee = sqrt(G*M_v*(1+e_ms)/a_ms);  %Mothership Apogee velocity, units = m/s
T_ms = 2*pi*sqrt(semimajor_ms^3/(G*M_v)); %CubeSat Period, units = s

fprintf('The Mothership travels at a speed of %f km/s at apogee and at %f km/s at perigee with a period of %f years.\n',v_ms_apogee/1000,v_ms_perigee/1000,T_ms*(3.17098*10^(-8)))

%% CubeSat Orbital Calculations

v_s_perigee = v_ms_perigee;
v_s_apogee = v_ms_apogee;
T_s = T_ms;
%v_s = sqrt(g*R_v^2/a_a); %CubeSat velocity, units = m/s
%T_s = 2*pi*sqrt(a_a^3/(G*M_v)); %CubeSat Period, units = s

fprintf('The CubeSat travels at a speed of %f km/s at apogee and at %f km/s at perigee with a period of %f years.\n',v_ms_apogee/1000,v_ms_perigee/1000,T_ms*(3.17098*10^(-8)))

%% Drag 

D = C_d*(A_ram/m_s)*rho_s*v_s_perigee^2;

fprintf('The CubeSat experiences a drag of %f at perigee.\n', D)

counter_a = 0;
alt_a = alt_initial;
a_initial = alt_a + R_v; %CubeSat orbital radius, units = m
    h_a = alt_a/1000; %CubeSat orbital radius, units = km
while alt_a >= alt_final %units = m
   
P_a =sqrt(4*pi^2*a_initial^3/(G*M_v));
    T_a = 900 + 2.5*(f10_a-70)+1.5*Ap_a;
    m_a = 27 - 0.012*(h_a-200);
    H_a = T_a/m_a;
    
    if alt_a > 50*10^3
        rho_a = rho;
    end
dPdt_a = -3*pi*a_initial*rho_a*C_d*A_ram/m_s;
    
% Calculating new P 
    P_a = P_a + dPdt_a;
    counter_a = counter_a +1;
    
    a_initial = ((P_a^2*G*M_v)/(4*pi^2))^(1/3);
    alt_a = a_initial-R_v; h_a = alt_a/1000;
end

fprintf('The CubeSat stays in orbit for %f years at an altitude of %f km.\n', counter_a*(3.17098*10^(-8)),alt_a/1000)

%% POSITION CALCULATIONS
% Reading Data Files
dataVOLE = xlsread('Positions_Vole.xls');

% Assigning Data Values to Variables
Vole_X = dataVOLE(1:13393,2);
    Vole_X = [Vole_X; Vole_X; Vole_X; Vole_X];
Vole_Y = dataVOLE(1:13393,3);
    Vole_Y = [Vole_Y; Vole_Y; Vole_Y;  Vole_Y];
Vole_Z = dataVOLE(1:13393,4);
    Vole_Z = [Vole_Z; Vole_Z; Vole_Z;  Vole_Z];

% One Orbit
% Initializing Matrices
time = zeros(1260,1); dist = zeros(1260,1); slantrange = zeros(1260,1);

slantrange = zeros(1260,1)+ max_altitude;
counter = 0;

% Computing Position
for k = 1:1260
	time(k) = (k-1)*10;
	dist(k) = (sqrt((Vole_X(k))^2+(Vole_Y(k))^2 +(Vole_Z(k))^2))-R_v/1000;

    if dist(k)<= slantrange(k)
    	counter = counter +10;
    end
end

% One Year
% Initializing Matrices
time_full = zeros(52560,1); dist_full = zeros(52560,1); slantrange_full= zeros(52560,1);

slantrange_full = zeros(52560,1)+ max_altitude;
counter_full = 0;

% Computing Position
for k = 1:52560
	time_full(k) = (k-1)*10;
	dist_full(k) = (sqrt((Vole_X(k))^2+(Vole_Y(k))^2 +(Vole_Z(k))^2))-R_v/1000;

    if dist_full(k)<= slantrange_full(k)
    	counter_full = counter_full +10;
    end
end

% Plotting 
fig =  figure(1);
    grid on      
    set(fig,'color','w');
subplot(1,2,1)
plot(time,smooth(dist),'-b')
hold on 
plot(time, slantrange,'--r')
xlim([0 1260]);
xlabel('Time [min]')
ylabel('Altitude [km]')
        set(gca,'FontSize',15)
title('Altitude of VOLE, One Orbit')

subplot(1,2,2)
plot(time_full,smooth(dist_full),'-b')
hold on 
plot(time_full, slantrange_full,'--r')
xlabel('Time [min]')
    xlim([0 length(Vole_X)])
ylabel('Altitude [km]')
        set(gca,'FontSize',15)
title('Altitude of VOLE, One Year')

fprintf('We have %f minutes of data collection oppurtunities in one orbit.\n',counter)
fprintf('We have %f minutes of data collection oppurtunities in one year.\n',counter_full)

%% POWER
P_tot = P_thermal + P_camera + P_OBC + P_radio + P_ACS + P_comm

S_V = 2601.3; %units = W/m^2

eff = 0.30;
A = 14*10^(-2);  %units = m^2

Psa = S_V*eff*A

%% PROPULSION
del_v = sqrt(2*KE/m_f); % units = m/s
fprintf('The CubeSat will be deployed at a del_v of %0.1f m/s.\n',del_v)

Isp = 2640; %units = s 
m_dot_prop = 48*10^(-9); %units = kg/s

thrust = Isp*m_dot_prop*g2; %units = kg*m/s2
del_t = del_v/ (thrust/ m_f); %units = s
fprintf('The propulsion will take %0.2f hours to achieve the del_v of %0.1f m/s.\n',del_t*0.000277778,del_v)
  
dist = del_v*del_t;
fprintf('In the time that it takes to achieve the del_v, VOLE will be at a distance of %0.1f m from the ISRO craft.\n',dist/2)
