function xdot=calc_xdot(t,x,g)
%Calculates the derivative of state variables

global V0 u0 C m0 tau

h = g(1);
w = g(2);

%Constant:
u0 = 1.2566e-6; %Wb/Am
%Resistance
Rc = 0.03; %ohms
Re = 0; %ohms
Rpe = 0; %ohms
%Inductance:
Lc = 3.5e-8; %Henry
Le = 0; %Henry
%Plasma, based on graphs
%1eV = 1.60218e-19 Joules
%k = 1.3807e-23 J/K
Te = (1.5 * 1.60218e-19)/1.3807e-23 ; %Kelvin
ne = 1e21; %1/m^3
%Current sheet mass
m0 = 1e-8; %kg
%Geometry:
% h = 0.02; %meter
% w = 0.01; %meter
%t_el = 0.002; %meter

Lpe = u0 *(h/w);
%Lpe = (u0/pi) * (1.5 + log (h/(w+t_el)));
%Where:
%m0 = current sheet mass at t=0
%u0 = permeability of free space
%h = electrode gap
%w = electrode width
Lpe_t = Lpe*x(1);
L_tot_t = Lc + Le + Lpe_t;
%Where:
%L_tot_t = total inductance of the circuit as function of time
%Lc = capacitor inductance
%Le = leads and wires inductance
%Lpe_t = parallel electrodes inductance as function of time
Rp = 8.08 * (h/(w*(Te^0.75)))  * sqrt((u0 * log (1.24e7 * ((Te^3) / ne )^0.5 ))/tau);
%Rp = .0130392;
%Where:
%Te = plasma temperature, use the graphical data
%ne=  electrone number density, use the graphical data
%tau = charateristic pulse time equals to 1/4 of the PPT ringing period
R_tot = Rc + Re + Rpe + Rp;
%Where:
%Rc = capacitor resistance
%Re =  leads and wires resistance
%Rpe = parallel electrodes resitance
%Rp = plasma resistance

xdot(1) = x(3); %Current Sheet Position
xdot(2) = x(4); %Capacitor Charge
xdot(3) = (0.5/m0)*(u0*h/w)*(x(4))^2; %Velocity
xdot(4) = (V0 - (x(2)/C) - (x(4)*R_tot) - (Lpe*x(3)*x(4)))/L_tot_t; %Current

%Where:
%C = capacitance
%V0 = initial discharge voltage of capacitor at time t=0

xdot = xdot';

end




%##########################
%Validation Data

% %Input: LES-6 PPT
% 
% u0 = 1.2566e-6; %Wb/Am
% V0 = 1360; %Volts
% C = 2e-6; %F
% %Resistance
% Rc = 0.03; %ohms
% Re = 0; %ohms
% Rpe = 0; %ohms
% %Inductance:
% Lc = 34e-9; %Henry
% Le = 0; %Henry
% %Plasma, based on graphs
% %1eV = 1.60218e-19 Joules
% %k = 1.3807e-23 J/K
% Te = (1.5 * 1.60218e-19)/1.3807e-23 ; %Kelvin
% ne = 1e21; %1/m^3
% tau = 4e-7; %sec
% %Current sheet mass
% m0 = 1e-8; %kg
% %Geometry:
% h = 0.03; %meter
% w = 0.01; %meter
% t_el = 0.002; %meter

% %Input: LES 8/9
% 
% u0 = 1.2566e-6; %Wb/Am
% V0 = 1538; %Volts
% C = 17e-6; %F
% %Resistance
% Rc = 0.03; %ohms
% Re = 0; %ohms
% Rpe = 0; %ohms
% %Inductance:
% Lc = 35e-9; %Henry
% Le = 0; %Henry
% %Plasma, based on graphs
% %1eV = 1.60218e-19 Joules
% %k = 1.3807e-23 J/K
% Te = (5 * 1.60218e-19)/1.3807e-23 ; %Kelvin
% ne = 1e21; %1/m^3
% tau = 4e-7; %sec
% %Current sheet mass
% m0 = 2.85e-8; %kg
% %Geometry:
% h = 0.0254; %meter
% w = 0.0254; %meter
% t_el = 0.002; %meter
