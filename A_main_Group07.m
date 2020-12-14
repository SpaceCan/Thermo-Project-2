%% Part 1
%% Declare Variables
clear;clc

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

%% Part 2
%% Declare Variables
clc;clear

m_aircraft = 1600;              % kg; mass of aircraft without fuel & equipment
eff_actual_cycle = 0.430218;    % efficiency of actual cycle
m_fuel = 400;                   % kg; mass of fuel

H_ethane = 47484;               % kJ/kg; heating value for ethane
M_mass_ethane = 30.070;         % g/mol; molecular mass of ethane
carbon_atoms_ethane = 2;        % # of carbon atoms in methanes chemical formula

H_methane = 50010;              % kJ/kg; heating value for methane
M_mass_methane = 16.043;        % g/mol; molecular mass of methane
carbon_atoms_methane = 1;        % # of carbon atoms in methanes chemical formula

H_propane = 46352;              % kJ/kg; heating value for propane
M_mass_propane = 44.094;        % g/mol; molecular mass of propane
carbon_atoms_propane = 3;        % # of carbon atoms in methanes chemical formula

H = [H_ethane H_methane H_propane];
M_mass_co2 = 44.0095;           % g/mol; molar mass of carbon dioxide
g = 9.8;                        % m/s^2
v = linspace(210,450,256);                  % km/hr; cruising speed of drone
rho_air = 1.007;                % kg/m^3; density of air at 2000 m altitude
A = 10;                         % m^2; planform area of UAV wing

load coefficients_chart
alpha = coefficients_chart(:,1);   % selects every value in the column for angle of attack alpha
C_L = coefficients_chart(:,2);     % selects every value in the column for coeff of lift
C_D = coefficients_chart(:,3);     % selects every value in the column for coeff of drag
%% Calculate Ranges
for i = 1:length(H)
    Q_released = m_fuel*(H(i)*1000);        % J; or kg*m^2/s^2; amount of energy from combustion
    lift = (m_fuel + m_aircraft)*g;             % Newtons; or kg*m/s^2
    
    for j = 1:length(v)
        coeff_L(i,j) = (2*lift)/(rho_air*A*(v(j)/3.6).^(2));  % units cancel - checked;
        alpha_actual(i,j) = interp1(C_L,alpha,coeff_L(i,j));    % degrees; corresponds to the alpha for the given C_L in the table
        coeff_D(i,j) = interp1(alpha,C_D,alpha_actual(i,j));    % corresponds to the C_D for the given alpha in the table
        drag(i,j) = 0.5*rho_air*A*(v(j)^2)*coeff_D(i,j);        % Newtons; or kg*m/s^2
        range(i,j) = ((Q_released*eff_actual_cycle)./drag(i,j))/1000;  % kilometers
    end
end
%% Calculating Emissions
moles_ethane = m_fuel*1000/M_mass_ethane;    %   calculates moles of ethane from given am't of fuel
moles_co2_ethane = moles_ethane*(carbon_atoms_ethane);   %  calculates moles of co2 using ratio
emissions_ethane = (moles_co2_ethane*M_mass_co2)/1000;   %   kg; carbon dioxide emissions

moles_methane = m_fuel*1000/M_mass_methane;    %   calculates moles of methane from given am't of fuel
moles_co2_methane = moles_methane*(carbon_atoms_methane);   %  calculates moles of co2 using ratio
emissions_methane = (moles_co2_methane*M_mass_co2)/1000;   %   kg; carbon dioxide emissions

moles_propane = m_fuel*1000/M_mass_propane;    %   calculates moles of propane from given am't of fuel
moles_co2_propane = moles_propane*(carbon_atoms_propane);   %  calculates moles of co2 using ratio
emissions_propane = (moles_co2_propane*M_mass_co2)/1000;   %   kg; carbon dioxide emissions

%% Plotting
figure('Units','inches','Position',[2 1 3.56 2.7].*1.5)
hold on
plot(v,range(2,:),'Color','#001275','LineWidth',1.5,'DisplayName','Methane')  %   could be plotted easier with plot(v,range)
plot(v,range(1,:),'Color','#F0A808','LineWidth',1.5,'DisplayName','Ethane')  %   plotted seperately only so i could choose the colors
plot(v,range(3,:),'Color','#EC0955','LineWidth',1.5,'DisplayName','Propane')
% title('UAV Range w.r.t its Cruising Speed and Fuel Type')
xlabel('Cruising Speed $$\left[\frac{km}{hr}\right]$$','Interpreter','latex')
ylabel('Range (km)','Interpreter','latex')
legend('Interpreter','latex')
leg = legend('show');       %   adds a title to the legend
title(leg,'Fuel Types')     %   adds a title to the legend
print('Graph-2','-r300','-djpeg') % Auto-export figure1 at a crisp 300 dpi
%% Part 3
%% Variable Declarations
clear;clc

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
