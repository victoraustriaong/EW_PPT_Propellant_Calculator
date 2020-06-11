%This program calculates the number of capacitors in parallel (given the 
%number of capacitors in series and capacitance of 1 capacitor) to achieve 
%the required capacitance.
%By: Victor Ong

%Parallel capacitor equivalent formula
%Ct = C1 + C2 + C3 + ... + Cn
%Series capacitor equivalent formula
%1/Ct = 1/C1 + 1/C2 + 1/C3 + ... + 1/Cn



req_cap = input('Enter required capacitance (Farad): ');               % Farad
cap_of1cap = input('Enter capacitance of 1 capacitor (Farad): ');      % Farad

num_caps_inSeries = input('Enter number capacitors in series (Farad): ');

series_equiv = 1 / ( (1/cap_of1cap*num_caps_inSeries));         % Farad

num_caps_inParallel = ceil(req_cap/series_equiv);               % Rounded up

%Check the actual capacitance
actual_cap = num_caps_inParallel * series_equiv                % Farad
total_num_capacitor = num_caps_inParallel * num_caps_inSeries