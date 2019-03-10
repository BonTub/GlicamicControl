clear all
close all
clc

%% Initialization
param = zeros(15,1);
param(1) = 0.066;
param(2) = 0.006;
param(3) = 0.06;
param(4) = 0.03;
param(5) = 0.138;
param(6) = 0.12;
param(7) = 0.16;
param(8) = 0.8;
param(9) = 40;
param(10) = 51.2*10^-4;
param(11) = 8.2*10^-4;
param(12) = 520*10^-4;
param(13) = 0.0161;
param(14) = 0.0097;
param(15) = 55;
    %Defining parameters
    k12 = param(1);
    ka1 = param(2);
    ka2 = param(3);
    ka3 = param(4);
    ke = param(5);
    VI = param(6);
    Vg = param(7);
    Ag = param(8);
    tmaxG = param(9);
    kb1 = param(10)*ka1;
    kb2 = param(11)*ka2;
    kb3 = param(12)*ka3;
    EGP0 = param(13);
    F01 = param(14);
    tmaxI = param(15);
    Dg=2;
    
%     %Formulized parameters
%     G = x(1)/Vg;
%     Ug = (Dg*Ag*t*exp(-t/tmaxG))/tmaxG^2;
%     UI = x(4)/tmaxI;
%     
%     if G >= 4.5
%         F01c = F01;
%     else
%         F01c = F01*G/4.5;
%     end
%     
%     if G >= 9
%         Fr = 0.003*(G-9)*Vg;
%     else
%         Fr = 0;
%     end
%% Matrices
%Operating point
Q1bar = 6*Vg;
Q2bar = 5*Vg;
x1bar = 0.04;
x2bar = 0.02;

A = [-x1bar k12 0 0 0 -Q1bar 0 -EGP0;
    x1bar -k12-x2bar 0 0 0 Q1bar -Q2bar 0;
    0 0 -1/tmaxI 0 0 0 0 0;
    0 0 1/tmaxI -1/tmaxI 0 0 0 0;
    0 0 0 1/(VI*tmaxI) -ke 0 0 0;
    0 0 0 0 kb1 -ka1 0 0;
    0 0 0 0 kb2 0 -ka2 0;
    0 0 0 0 kb3 0 0 -ka3];
B = [0 1 1; 0 0 0; 1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0];
C = [1/Vg 0 0 0 0 0 0 0];
D = [0 0 0];

System = ss(A,B,C,D);
%% Simulation
tspan = 0:1:2000;
Dg = 5;
delay = 0;
for i = 1:length(tspan)
    t = tspan(i);
    Ug(i) = max(0,(Dg*Ag*(t-delay)*exp(-(t-delay)/tmaxG))/tmaxG^2);
end
u = 0.1*ones(size(Ug));
F01c = F01*ones(size(Ug));
Fr = zeros(size(Ug));
EGP0_vectorized = EGP0*ones(size(Ug));
Constant_input = F01c+Fr+EGP0_vectorized;
x0 = [1.04
    0.485
    0
    0
    0
    0
    0
    0];

[x,t] = lsim(System,[u; Ug; Constant_input], tspan, x0);
figure
plot(t,x)