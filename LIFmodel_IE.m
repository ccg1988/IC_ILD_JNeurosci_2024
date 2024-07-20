function [spikes, V] = LIFmodel_IE(Ge, Gi, noise_magnitude, Erest)

% modified from Bendor, PLOS Computational Biology, 2015
% this function will be called by function "F_percent_IE" and program "F_demo_IE"

spikes=[]; V=[]; t=1; i=1;
step=.0001;   % 0.1 ms duration (temporal increment for running simulation)
% noise_sigma = 0.01 ; % [V] scale of noise, default is 10mV
noise_sigma = 0.00 ; % [V] scale of noise, default is 10mV
C = 0.25*1e-9;          %0.25 nF (Wehr and Zador); 10 ms time constant(Cm*Rm)
Grest = 25*1e-9;      %25 nS (Wehr and Zador)
Ee = 0; % 0 mV (Xiong 2013)
Ei = -0.070;  %-70 mV (Xiong 2013)
V(1) = Erest;

%avoid negative conductances
Ge=Ge+noise_magnitude*randn(1, length(Ge)); %mean is 0
Gi=Gi+noise_magnitude*randn(1, length(Gi));
Ge(Ge<0)=0;
Gi(Gi<0)=0;

while(t<(length(Ge)-0))
    V(t+1)=(-step*(Ge(t)*(V(t)-Ee)+Gi(t)*(V(t)-Ei)+Grest*(V(t)-Erest))/C)+V(t)+noise_sigma*randn*sqrt(step); % 
    if V(t+1)>(Erest+0.020) %20 mV above Erest as spike threshold; more negative, more spikes
        spikes(i)=step*(t+1);
        i=i+1;
        t=t+1;
        V(t+1)=Erest;
    end
    t=t+1;
% V(end) = [];
end
