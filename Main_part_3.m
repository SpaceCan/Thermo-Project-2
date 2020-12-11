%% Variable Declarations
clear;clc;close all

m_aircraft = 1600;              % kg; mass of aircraft without fuel & equipment (aka dry mass)
eff_actual_cycle = 0.430218;    % efficiency of actual cycle
m_fuel = linspace(600,0,601);   % kg; mass of fuel as it is consumed

H_methane = 50010;              % kJ/kg; heating value for methane
M_mass_methane = 16.043;        % g/mol; molecular mass of methane
carbon_atoms_methane = 1;       % # of carbon atoms in methanes chemical formula
M_mass_co2 = 44.0095;           % g/mol; molar mass of carbon dioxide

g = 9.8;                        % m/s^2
v = linspace(210,450,256);      % km/hr; cruising speed of drone
rho_air = 1.007;                % kg/m^3; density of air at 2000 m altitude
A = 10;                         % m^2; planform area of UAV wing

load coefficients_chart
alpha = coefficients_chart(:,1);   % selects every value in the column for angle of attack alpha
C_L = coefficients_chart(:,2);     % selects every value in the column for coeff of lift
C_D = coefficients_chart(:,3);     % selects every value in the column for coeff of drag
%% Range Calculation
% Declaring variables before the loop to save memory
range = zeros(1,length(v));
drag = zeros(1,length(v));
coeff_D = zeros(1,length(v));
alpha_actual = zeros(1,length(v));
coeff_L = zeros(1,length(v));

for i = 1:length(m_fuel) % Integrating fuel as it decreases
    Q_released = m_fuel(end-1)*(H_methane*1000);        % J; or kg*m^2/s^2; amount of energy from combustion
    lift = (m_fuel(i) + m_aircraft)*g;             % Newtons; or kg*m/s^2
    
    for j = 1:length(v)
        coeff_L(1,j) = (2*lift)/(rho_air*A*(v(j)/3.6).^(2));  % units cancel - checked;
        alpha_actual(1,j) = interp1(C_L,alpha,coeff_L(1,j));    % degrees; corresponds to the alpha for the given C_L in the table
        coeff_D(1,j) = interp1(alpha,C_D,alpha_actual(1,j));    % corresponds to the C_D for the given alpha in the table
        drag(1,j) = 0.5*rho_air*A*(v(j)^2)*coeff_D(1,j);        % Newtons; or kg*m/s^2
        range(1,j) = range(1,j) + ((Q_released*eff_actual_cycle)./drag(1,j))/1000;  % kilometers
    end
end
%% Calulating CO2 usage vs range
moles_ethane = max(m_fuel)*1000/M_mass_methane;    %   calculates moles of ethane from max fuel
moles_co2_methane = moles_ethane*(carbon_atoms_methane);   %  calculates moles of co2 using ratio
emissions_methane = (moles_co2_methane*M_mass_co2)/1000;   %   kg; carbon dioxide emissions
range_emissions = emissions_methane./(range(1,:));    %   kg/km carbon emissions per kilometer vs cruising speed

%% Plotting
figure('Units','inches','Position',[2 1 3.56 2.7].*1.5)
plot(v,range(1,:),'Color','#EC0955','LineWidth',1.5)  %   plotted seperately only so i could choose the colors
%title('Drone Range w.r.t its Cruising Speed for non-const fuel mass')
xlabel('Cruising Speed $$\left[\frac{km}{hr}\right]$$','Interpreter','latex')
ylabel('Range (km)','Interpreter','latex')
legend('600 kg','Interpreter','latex')
leg = legend('show');       %   adds a title to the legend
title(leg,'Mass of Fuel')   %   adds a title to the legend
print('Graph-3','-r300','-djpeg') % Auto-export figure1 at a crisp 300 dpi

figure('Units','inches','Position',[2 1 3.56 2.7].*1.5)
plot(v,range_emissions,':','Color','k','LineWidth',2)
xlabel('Cruising Speed $$\left[\frac{km}{hr}\right]$$','Interpreter','latex')
ylabel('$CO_{2}$ emissions $$\left[\frac{kg}{km}\right]$$','Interpreter','latex')
print('Graph-4','-r300','-djpeg') % Auto-export figure2 at a crisp 300 dpi
