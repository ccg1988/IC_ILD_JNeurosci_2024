function [spike_distribution, BW, rate_peak, ILD50]=BW_SYNnum(synapse_num_In,BW_type)
% change the # of inhibitory inputs for "ILD" or "Hz"
% modified from "BW_SYNnum_example"
% by CCG @ 2023-03-14

spike_num_Poi = 7 ; % Number of Possion spikes; use '7' to reduce firing rate
noise_magnitude = 1.5e-8 ; %default noise level in conductance--decide firing rate (>Ge and Gi)
IE_delay = 2.0 ; % + means Ex lead In, longer delay==more spikes (Xiong 2013)
x = 1.0:0.005:5; x(end)=[];
synapse_num_Ex = 10 ;    % N excitatory synaptic inputs

if strcmp(BW_type, 'ILD')
    mean_Ex=3.5; std_Ex=0.75; % mu=Mean of logarithmic values
    mean_In=1.5; std_In=1; % sigma=Scale parameter of logarithmic values
    E_strength = 1.2+flip(mean_Ex*pdf('LogLogistic', x, 1, std_Ex));
    I_strength = 1.2+flip(mean_In*pdf('LogLogistic', x, 1, std_In));
elseif strcmp(BW_type, 'Hz')
    mean_Ex=1.5; std_Ex=0.5;
    mean_In=3.5; std_In=1;%raw
    E_strength = 1+transpose(mean_Ex*normpdf(x,3,std_Ex)); %0.5sd in paper
    I_strength = 1+transpose(mean_In*normpdf(x,3,std_In)); %1.0sd in paper---raw
end
nreps = numel(x);
step=.0001; % [S]
Erest = -0.065 ; % more negative, more spikes---FIXED
kernel_time_constant=.005;  %time constant of 5 ms (Wehr and Zador)
stimulus_duration=0.5;  %half second
spike_distribution=NaN(1, nreps);
t=step:step:(kernel_time_constant*10);
kernel=t.*exp(-t/kernel_time_constant); %1D vector  a=1/kernel_time_constant
kernel=1e-9*kernel/max(kernel); %amplitude of 1 nS (shape of single spike)

input=zeros(size(step:step:stimulus_duration));
stimulus_input_length=length(step:step:stimulus_duration);
for r = 1 : nreps
    E_input = input;
    I_input = input;
    spike_num = poissrnd(spike_num_Poi);
    spike_gap = floor(stimulus_input_length/spike_num) ;
    for i = 1 :  spike_num
        time_window_end = (i-1)*spike_gap + length(kernel);
        time_window_start = 1 + time_window_end-length(kernel);
        time_window = time_window_start : time_window_end ;
        if max(time_window)<=stimulus_input_length
        %assign curve at specific point range     
        E_input (time_window) =  E_input (time_window)+kernel*synapse_num_Ex; 
        I_input (time_window) = I_input (time_window)+kernel*synapse_num_In;
        end
    end
    
    delay=round(abs(IE_delay)/(1000*step));  % delay in steps, single value
    delay_chunk=zeros(1,delay); %1D vector
    Ge=E_input * E_strength(r); %number of kernels * kernel amplitude(A; nS)
    Gi=[delay_chunk I_input(1:(length(I_input)-delay))] * I_strength(r);
 
    %Excitatory and Inhibitory conductance (nS)!!!
    [spikes, ~] = LIFmodel_IE (Ge, Gi, noise_magnitude, Erest); 
    spike_distribution(r)=numel(spikes)/stimulus_duration;
end
rate_sm=smooth(1:nreps, spike_distribution,0.1,'sgolay');
rate_peak=max(rate_sm); rate_hf=rate_peak/2;
[~,ILD50]=min(abs(rate_sm-rate_hf));
BW=(100*numel(find(rate_sm>=rate_hf))/nreps);
