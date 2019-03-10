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
tspan=0:1:3000;
u=0.1;
f = @(t,x)GlucoseModel(t,x,u,param);
x0 = [1.04
    0.485
    0
    0
    0
    0
    0
    0];

[t,x]=ode45(f, tspan, x0);
figure
hold on
for i = 1:1
    y(:,i)= x(:,i)/param(7);
plot(t,y(:,i))
end
%% Generating figures
