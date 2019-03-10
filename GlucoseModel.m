function [xdot] = GlucoseModel(t,x,u,param)
    
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
    Dg=0;
    %u = u*exp(-0.1*t);
    
    %Formulized parameters
    G = x(1)/Vg;
    Ug = (Dg*Ag*t*exp(-t/tmaxG))/tmaxG^2;
    UI = x(4)/tmaxI;
    
    if G >= 4.5
        F01c = F01;
    else
        F01c = F01*G/4.5;
    end
    
    if G >= 9
        Fr = 0.003*(G-9)*Vg;
    else
        Fr = 0;
    end

    %System equations
    xdot=zeros(8,1);
    xdot(1) = -F01c-x(6)*x(1)+k12*x(2)-Fr+Ug+EGP0*(1-x(8)); %dQ1/dt
    xdot(2) = x(6)*x(1)-(k12+x(7))*x(2); %dQ2/dt
    xdot(3) = u-x(3)/tmaxI; %dS1/dt;
    xdot(4) = x(3)/tmaxI - x(4)/tmaxI; %dS2/dt
    xdot(5) = UI/VI-ke*x(5);
    xdot(6) = -ka1*x(6)+kb1*x(5);
    xdot(7) = -ka2*x(7)+kb2*x(5);
    xdot(8) = -ka3*x(8)+kb3*x(5);
    
    
end