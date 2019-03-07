clear all
close all
clc

%% Defining parameters;
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

%% Model simulation
tspan=0:1:2000;
u=11;
f = @(t,x)GlucoseModel(t,x,u,param);
x0 = [5*param(7)
    4*param(7)
    0.01
    0.01
    0
    param(10)
    param(11)
    param(12)];

[t,x]=ode45(f, tspan, x0);
y = x(:,1)./param(7);
%% Generating figures
figure
hold on
plot(t,y)
legend('Stable Dead')