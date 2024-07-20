% modified from program used in my stimulus-specific facilitation paper
% LIFmodel_IE is modified from LIFmodel_basic (Bendor, PLOS Computational Biology, 2015)
% some para. from Xiong 2013, some para. are still using previous one
% it will only call "LIFmodel_IE"

% it include BW for both ILD (LogLogistic) and Hz (normpdf); only change synapse number
% it could shift the postion of In. tuning curve for unmatched cases 
% by CCG @ 2023-03-13

clear; clc; 
% close all

synapse_num_In = 5 ;    %--------------------------------change synapse number
tuning_type='ILD'; % 'Hz' or 'ILD'
match_YN='Y'; % 'match' or 'not match'

spike_num_Poi = 7 ; % Number of Possion spikes
noise_magnitude = 1.5e-8 ; %default noise level in conductance--decide firing rate (>Ge and Gi)
IE_delay = 2.0 ; % + means Ex lead In, longer delay==more spikes (Xiong 2013)
x = 1.0:0.005:5; x(end)=[];
synapse_num_Ex = 10 ;    % N excitatory synaptic inputs
%%%%%%%%%%%%%%%%%%%%%%% ILD %%%%%%%%%%%%%%%%%%%%%%%
if strcmp(tuning_type, 'ILD')
    %mu=Mean of logarithmic values; sigma=Scale parameter of logarithmic values
    mean_Ex=3.5; std_Ex=0.75; 
    mean_In=1.5; std_In=1; %-- synapse strength; KEEP it fixed
    E_strength=1.2+flip(mean_Ex*pdf('LogLogistic',x,1, std_Ex));
if strcmp(match_YN, 'Y')
    I_strength=1.2+flip(mean_In*pdf('LogLogistic',x,1, std_In));%for matched EX vs. IN
elseif strcmp(match_YN, 'N')
    I_strength=1.2+flip(mean_In*pdf('LogLogistic',x,1+0.5, std_In));%--------for unmatched
end
% figure;plot(x,E_strength);hold on; plot(x,I_strength);legend({'Ex','In'}); 
% xlim([min(x) max(x)]);% ylim(1+[0 1.5])

%%%%%%%%%%%%%%%%%%%%%%% Frequency %%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(tuning_type, 'Hz')
    mean_Ex=1.5; std_Ex=0.5; 
    mean_In=3.5; std_In=1;%-- synapse strength; KEEP it fixed
    E_strength = 1+transpose(mean_Ex*normpdf(x,3,std_Ex)); 
if strcmp(match_YN, 'Y')
    I_strength = 1+transpose(mean_In*normpdf(x,3,std_In)); % for matched EX vs. IN
elseif strcmp(match_YN, 'N')
    I_strength = 1+transpose(mean_In*normpdf(x,3+1.5,std_In));% for unmatched EX vs. IN
end
% figure;plot(x,E_strength);hold on; plot(x,I_strength);legend({'Ex','In'}); 
% xlim([min(x) max(x)]);ylim(1+[0 1.5])
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
%%
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
        %time_window(time_window>stimulus_input_length)=[];
        %last spike's kernal is usually outside of the sound duration
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

%%
figure;tiledlayout(2,1);
nexttile
plot(E_strength);hold on;plot(I_strength);legend({'Ex','In'})
ylabel('amplitude'); xlabel('# trial'); 
title(['mean E',num2str(mean_Ex),'  mean In',num2str(mean_In),...
    '  std E',num2str(std_Ex),'  std In',num2str(std_In)])
if strcmp(match_YN, 'Y')
    text(100, 1.8, 'match Ex & In curves'); % for matched EX vs. IN
elseif strcmp(match_YN, 'N')
    text(100, 1.8, 'unmatch Ex & In curves');% for unmatched EX vs. IN
end
% nexttile
% hist_range=min(I_strength):0.05:max(E_strength);
% histogram(E_strength,hist_range);hold on;histogram(I_strength,hist_range);
% legend({'Ex','In'}); 
% title(['# In syn.=',num2str(synapse_num_In),'  # Ex syn.=',num2str(synapse_num_Ex)])
% xlabel('amplitude'); ylabel('# counts')
nexttile
plot(spike_distribution,'Color','k','Marker', '.', 'MarkerSize',5,...
    'LineStyle','none');hold on; % 'LineStyle','none'
rate_sm=smooth(1:nreps, spike_distribution,0.1,'sgolay');
rate_hf=max(rate_sm)/2;
[~,ILD50]=min(abs(rate_sm-rate_hf));
plot(rate_sm,'LineWidth',2); 
plot(1:nreps,rate_hf*ones(nreps,1),'Color',rgb('DarkGreen'))
wid=round(100*numel(find(rate_sm>=rate_hf))/nreps);
title(['Half-peak tuning width=',num2str(wid), ...
    '%  peak firing=',num2str(max(rate_sm)),'spikes/s',...
    '   50% ILD=',num2str(ILD50)])
ylabel('spikes/s'); xlabel('# trial'); ylim([0 inf])