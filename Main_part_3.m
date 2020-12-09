clear;clc;close all

m_aircraft = 1600;              % kg; mass of aircraft without fuel & equipment (aka dry mass)
eff_actual_cycle = 0.430218;    % efficiency of actual cycle
m_fuel = linspace(600,0,601);   % kg; mass of fuel as it is consumed
H_ethane = 47484;               % kJ/kg; heating value for ethane
g = 9.8;                        % m/s^2
v = linspace(210,450,256);      % km/hr; cruising speed of drone
rho_air = 1.007;                % kg/m^3; density of air at 2000 m altitude
A = 10;                         % m^2; planform area of UAV wing

load coefficients_chart
alpha = coefficients_chart(:,1);   % selects every value in the column for angle of attack alpha
C_L = coefficients_chart(:,2);     % selects every value in the column for coeff of lift
C_D = coefficients_chart(:,3);     % selects every value in the column for coeff of drag

% Declaring variables before the loop to save memory
range = zeros(1,length(v));
drag = zeros(1,length(v));
coeff_D = zeros(1,length(v));
alpha_actual = zeros(1,length(v));
coeff_L = zeros(1,length(v));

for i = 1:length(m_fuel)
    Q_released = m_fuel(end-1)*(H_ethane*1000);        % J; or kg*m^2/s^2; amount of energy from combustion
    lift = (m_fuel(i) + m_aircraft)*g;             % Newtons; or kg*m/s^2
    
    for j = 1:length(v)
        coeff_L(1,j) = (2*lift)/(rho_air*A*(v(j)/3.6).^(2));  % units cancel - checked;
        alpha_actual(1,j) = interp1(C_L,alpha,coeff_L(1,j));    % degrees; corresponds to the alpha for the given C_L in the table
        coeff_D(1,j) = interp1(alpha,C_D,alpha_actual(1,j));    % corresponds to the C_D for the given alpha in the table
        drag(1,j) = 0.5*rho_air*A*(v(j)^2)*coeff_D(1,j);        % Newtons; or kg*m/s^2
        range(1,j) = range(1,j) + ((Q_released*eff_actual_cycle)./drag(1,j))/1000;  % kilometers
    end
end
%%
hold on
plot(v,range(1,:),'b')  %   plotted seperately only so i could choose the colors
title('Drone Range w.r.t its Cruising Speed for non-const fuel mass')
xlabel('Cruising Speed (km/hr)')
ylabel('Range (km)')
legend('600 kg')
leg = legend('show');       %   adds a title to the legend
title(leg,'Mass of Fuel')   %   adds a title to the legend
