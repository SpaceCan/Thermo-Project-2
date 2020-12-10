%% Declare Variables
clear;clc;close all

m_aircraft = 1600;              % kg; mass of aircraft without fuel & equipment
eff_actual_cycle = 0.430218;    % efficiency of actual cycle
m_fuel = [200 400 600];         % kg; mass of fuel
H_ethane = 47484;               % kJ/kg; heating value for ethane
g = 9.8;                        % m/s^2
v = linspace(210,450,256);                  % km/hr; cruising speed of drone
rho_air = 1.007;                % kg/m^3; density of air at 2000 m altitude
A = 10;                         % m^2; planform area of UAV wing

load coefficients_chart
alpha = coefficients_chart(:,1);   % selects every value in the column for angle of attack alpha
C_L = coefficients_chart(:,2);     % selects every value in the column for coeff of lift
C_D = coefficients_chart(:,3);     % selects every value in the column for coeff of drag
%% Calculate Ranges
for i = 1:length(m_fuel)
    Q_released = m_fuel(i)*(H_ethane*1000);        % J; or kg*m^2/s^2; amount of energy from combustion
    lift = (m_fuel(i) + m_aircraft)*g;             % Newtons; or kg*m/s^2
    
    for j = 1:length(v)
        coeff_L(i,j) = (2*lift)/(rho_air*A*(v(j)/3.6).^(2));  % units cancel - checked;
        alpha_actual(i,j) = interp1(C_L,alpha,coeff_L(i,j));    % degrees; corresponds to the alpha for the given C_L in the table
        coeff_D(i,j) = interp1(alpha,C_D,alpha_actual(i,j));    % corresponds to the C_D for the given alpha in the table
        drag(i,j) = 0.5*rho_air*A*(v(j)^2)*coeff_D(i,j);        % Newtons; or kg*m/s^2
        range(i,j) = ((Q_released*eff_actual_cycle)./drag(i,j))/1000;  % kilometers
    end
end
%% Plotting
figure('Units','inches','Position',[2 1 3.56 2.7].*1.5)
hold on
plot(v,range(1,:),'Color','#001275','LineWidth',1.5) %   plotted seperately only so i could choose the colors
plot(v,range(2,:),'Color','#1DDDBA','LineWidth',1.5) %   could be plotted easier with plot(v,range)
plot(v,range(3,:),'Color','#EC0955','LineWidth',1.5)
% title('Drone Range w.r.t its Cruising Speed and Amount of Fuel')
xlabel('Cruising Speed $$\left[\frac{km}{hr}\right]$$','Interpreter','latex')
ylabel('Range (km)','Interpreter','latex')
legend('200 kg','400 kg','600 kg','Interpreter','latex')
leg = legend('show');       %   adds a title to the legend
title(leg,'Mass of Fuel')   %   adds a title to the legend
print('Graph-1','-r300','-djpeg') % Auto-export figure1 at a crisp 300 dpi
