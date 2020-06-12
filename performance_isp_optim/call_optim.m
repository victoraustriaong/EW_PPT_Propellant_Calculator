clear;
close;
clc;

LB = [1;1];
UB = [2;2];
Aeq=[1 1];
Beq=[0];

X0=[1;1];


options = optimoptions(@fmincon, 'ConstraintTolerance',1e-12);

%[X isp] = fmincon(@ISP_function,X0,[],[],Aeq,Beq,LB,UB,[],options);
[X isp]= fmincon(@ISP_function,X0,[],[],[],[],LB,UB,[],options);

