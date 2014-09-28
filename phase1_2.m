%Timecourse - Simulation time with timesteps 
simulT = 100; %miliseconds
dT = .01;
t= 0:dT:simulT;



%Injected Current ( affects the change in voltage) - Step Pulse
I = zeros(1, length(t))
%gives step function at I
for j = 1:50:(length(t)-1)./2
    I(1, (2.*j)-1:(2.*j)+50) = 5  %indexing 
end


%Parameter Constants
gbar_K = 36;
gbar_Na = 120;
gbar_L = 0.3;
E_K = -12;
E_Na = 115;
E_L = 10.6;
C_m = 1;

%Initial Condition & Gating Variables  
V_m = 0;
alpha_m = 0.1 * ( (25-V_m) / (exp((25-V_m)/10)-1) ); 
beta_m = 4 * exp(-V_m/18);
alpha_n = 0.01 * ( (10-V_m) / (exp((10-V_m)/10)-1) );
beta_n = 0.125 * exp(-V_m/80) ;
alpha_h = 0.07 * exp(-V_m/20) ; 
beta_h = 1/ (exp((30-V_m)/10) +1);

m(1) = alpha_m / (alpha_m + beta_m);
n(1) = alpha_n / (alpha_n + beta_n);
h(1) = alpha_h / (alpha_h + beta_h);

for i = 1:numel(t)-1
    
    %Gating variables at each timestep
    alpha_m(i) = 0.1 * ( (25-V_m(i)) / (exp((25-V_m(i))/10)-1) ); 
    beta_m(i) = 4 * exp(-V_m(i)/18);
    alpha_n(i) = 0.01 * ( (10-V_m(i)) / (exp((10-V_m(i))/10)-1) );
    beta_n(i) = 0.125 * exp(-V_m(i)/80) ;
    alpha_h(i) = 0.07 * exp(-V_m(i)/20) ; 
    beta_h(i) = 1/ (exp((30-V_m(i))/10) +1);
    
    %Currents
    I_Na = (m(i)^3 * gbar_Na * h(i) * (V_m(i) - E_Na));
    I_K = (n(i)^4 * gbar_K * (V_m(i) - E_K));
    I_L = gbar_L * (V_m(i) - E_L);
    I_ion = I(i) - I_K - I_Na - I_L;
    
    %First order approximation using Euler's Method
    V_m(i+1) = V_m(i) + dT * I_ion/C_m;
    n(i+1) = n(i) + dT * (alpha_n(i) * (1-n(i)) - beta_n(i) * n(i));
    m(i+1) = m(i) + dT * (alpha_m(i) * (1-m(i)) - beta_m(i) * m(i));
    h(i+1) = h(i) + dT * (alpha_h(i) * (1-h(i)) - beta_h(i) * h(i));
    
    
end

V_m = V_m - 70;

%Plot Membrane Potential.
subplot(1, 2, 1)
plot (t, V_m, 'LineWidth', 3) 
hold on
legend ({'Voltage'})
xlabel('Time (ms)')
ylabel('Voltage (mv)')
title('Membrane Potential')

%Plot Conductances.
subplot(1, 2, 2)
P1 = plot(t,gbar_K * n.^4, 'Linewidth', 2 );
hold on 
P2 = plot(t, gbar_Na * (m.^3).*h, 'r', 'LineWidth', 2);
legend([P1,P2], 'gK', 'gNa')
xlabel('Time (ms)')
ylabel('Conductance (S/cm^2)')
title ('gK and gNa')