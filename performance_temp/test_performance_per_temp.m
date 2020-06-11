clear
clc
close


percent = -15:1:15;
N=length(percent);
ISP=[];
Ibit=[];
nominal_cap = 3e-6;


for i=1:N
    
    cap = (percent(i)/100 +1)*nominal_cap;
    
    [ISP(i) Ibit(i)]= performance_per_temp(cap);
    
end

ISP=ISP';
Ibit=Ibit';