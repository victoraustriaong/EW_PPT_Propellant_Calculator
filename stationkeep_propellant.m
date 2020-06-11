%This program calculates the propellant required from the delta-v required for East-West
%Stationkeeping mission.
%By: Victor Ong
clear, clc, close all;

mCube = input('Enter nanoSat mass(kg): '); %kg
Isp = input('Enter Isp (s): '); %s

% LEO East-West Station Keeping Delta-V requirements according to 
%Ref: 
DeltaVstation = 0.1162; %m/s every 21.21 days
DeltaVstation = DeltaVstation * (365/21.21);%Every year

mPropellant = (mCube*DeltaVstation)/(9.81*Isp); %kg

fprintf('Propellant required = %.2f grams \n', mPropellant*1000)
rhoTef = 2.2; %g/cm^3
volPropellant = mPropellant*1000/rhoTef; %cm^3
fprintf('For PPT propulsion with Teflon propellant, propellant volume is %.2f cm^3 \n', volPropellant)