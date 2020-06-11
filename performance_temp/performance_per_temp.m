function [ISP Ibit] = performance_per_temp(cap)

g=[2 1];
global V0 u0 C m0 tau

%g for the geometry input

%Input:
V0 = 1500; %Volts
C =cap; %F
tau = 4e-7; %Pulse time (sec)

%Solver parameters
finalTime = 1.8e-5;
timesteps = 1000;
tolerance = 1e-6;
tt=linspace(0,finalTime,timesteps);

options = odeset('AbsTol', tolerance);
[t,x] = ode45(@(t,x) calc_xdot(t,x,g), tt, [0 0 0 0],options);

%Results as function of time
position = x(:,1);
voltage = V0-(x(:,2)./17e-6);
velocity = x(:,3);
current = x(:,4);

%Performance
channelExitTime = interp1(position,t,0.025);
channelExitVelocity = interp1(t,velocity,channelExitTime);
ISP = channelExitVelocity/9.81;
Ibit = channelExitVelocity*m0;
thrust_effiency = (m0*(channelExitVelocity^2) / (C*V0^2))*100;


E = 0.5*C*(voltage.^2);
% E_consumed = 0.5*C*(V0^2) - interp1(t,E,channelExitTime);
% eff = (9.81*Isp*Ibit) / (E_consumed*2) *100
% TP = Ibit/E_consumed

%plot(t,E)

end

