%Timecourse - Simulation time with timesteps 
simulT = 100; %miliseconds
dT = .01; %Time Interval
t= 0:dT:simulT;  %Starting from 0 to 100 at timestep of dT


%Injected Current ( affects the change in voltage) - Constant Current
I = zeros(1, length(t))    
I(1:end) = 5.*pi.*(1e-8).*(1e3); %Setting from 1 to end; 
%the corresponding values of the current input for each timestep is set to five.
%the given unit of microAmp/cm^2 is converted to microAmps by multiplying
%the given value by cell surface area (sphere: 4*pi*R^2) with diameter of 1 micron^2 (Diameter
%of human nerve cell as referenced by Wikipedia)    


%Parameter Constants (Provided in the Modeling Phase1 document)
gbar_K = 36;   %maximum Conductance of K+
gbar_Na = 120; %maximum Conductance of Na+
gbar_L = 0.3; %maximum Conductance of Leakage
E_K = -12; %Nernst Voltage of K+
E_Na = 115; %Nernst Voltage of NA+
E_L = 10.6; %Nernst Voltage of Leakage
C_m = 1; %Membrane Capacitance

%Initial Condition & Gating Variables (Provided in the Modeling Phase1 document)
V_m = 0;    %initial voltage set to zero

%gating variables have no units
alpha_m = 0.1 * ( (25-V_m) / (exp((25-V_m)/10)-1) ); %probability of alpha_m 
beta_m = 4 * exp(-V_m/18); %probability of beta_m 
alpha_n = 0.01 * ( (10-V_m) / (exp((10-V_m)/10)-1) ); %probability of alpha_n 
beta_n = 0.125 * exp(-V_m/80) ; %probability of beta_n
alpha_h = 0.07 * exp(-V_m/20) ;  %probability of alpha_h 
beta_h = 1/ (exp((30-V_m)/10) +1); %probability of beta_h 

m(1) = alpha_m / (alpha_m + beta_m); %probability of channel open and let ions flow for m
n(1) = alpha_n / (alpha_n + beta_n);%probability of channel open and let ions flow for n
h(1) = alpha_h / (alpha_h + beta_h);%probability of channel inactivation

for i = 1:length(t)-1 %For indexing until 10000, one has to substract 1.
    %for loop is used to take information from previous timestep and also
    %multiplying timestep by given function and constant.
    %Gating variables at each timestep 
    %The same gating variables listed above. The variables are calculated
    %at each timestep.
    alpha_m(i) = 0.1 * ( (25-V_m(i)) / (exp((25-V_m(i))/10)-1) ); 
    beta_m(i) = 4 * exp(-V_m(i)/18);
    alpha_n(i) = 0.01 * ( (10-V_m(i)) / (exp((10-V_m(i))/10)-1) );
    beta_n(i) = 0.125 * exp(-V_m(i)/80) ;
    alpha_h(i) = 0.07 * exp(-V_m(i)/20) ; 
    beta_h(i) = 1/ (exp((30-V_m(i))/10) +1);
    
    %Currents (Provided in the Modeling Phase1 document)
    I_Na = (m(i)^3 * gbar_Na * h(i) * (V_m(i) - E_Na)); %Current for Na+ ions
    I_K = (n(i)^4 * gbar_K * (V_m(i) - E_K));%Current for K+ ions
    I_L = gbar_L * (V_m(i) - E_L);%Current for leakage
    I_ion = I(i) - I_K - I_Na - I_L;%Current for ions
    
    %First order approximation using Euler's Method
    V_m(i+1) = V_m(i) + dT * I_ion/C_m;
    n(i+1) = n(i) + dT * (alpha_n(i) * (1-n(i)) - beta_n(i) * n(i));
    m(i+1) = m(i) + dT * (alpha_m(i) * (1-m(i)) - beta_m(i) * m(i));
    h(i+1) = h(i) + dT * (alpha_h(i) * (1-h(i)) - beta_h(i) * h(i));
    %The values of n(i+1), m(i+1), and h(i+1) are dependent on the previous
    %values of n(i), m(i), and h(i) respectively. 
    %These correspond to the derivates listed on the Modeling Phase1
    %document

    
end

V_m = V_m - 70;   %The resting membrane potential is set by subtracting 70mV from the
%intial voltage value

%Plot Membrane Potential.
subplot(1, 2, 1)  %Refer to the left figure  
plot (t, V_m)  %Plot membrane voltage over duration of time 
legend ({'Voltage'}) % Boxed fonts representing the line of graph
xlabel('Time (ms)')  %X-axis with Time in miliseconds
ylabel('Voltage (mv)') %Y-axis with Voltage in milivolts
title('Membrane Potential') % Title of the graph

%Plot Conductances.
subplot(1, 2, 2) %Refer to the right figure 
P1 = plot(t,gbar_K * n.^4); %Plot K+ conductance over duration of time
hold on %Allow to hold two distinguished lines on the same graph.
P2 = plot(t, gbar_Na * (m.^3).*h, 'r'); %Plot Na+ conductance over duration of time
legend([P1,P2], 'gK', 'gNa') % Boxed fonts representing each line of graph
xlabel('Time (ms)')   %X-axis with Time in miliseconds
ylabel('Conductance (mS/cm^2)') %Y-axis with Conductance in S/cm^2
title ('gK and gNa')  % Title of the graph