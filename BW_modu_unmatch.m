% change the unmatched inhibitory synapse number&strength
% check the modulation of tuning bandwidth and firing rate
% it will call: BW_SYNnum_unmatch, BW_SYNamp_unmatch, and bonf_holm
% we need to save "unmatch_Hz & unmatch_ILD" in the last cell
% these two values will be used in "plot_match_unmatch_ILD_Hz"
% by CCG @ 2023-03-14
   
clear;clc
BW_type='ILD'; % either 'Hz' or 'ILD'
Rep=10; % We used 50 in the paper; 10 and 50 have similar results but need 1/5 time
%% number of inhibitory synapse changes compute---This cell takes minutes depends on "Rep"
Syn_num_In=1:10; N_syn=numel(Syn_num_In);
BW_num=nan(N_syn,Rep); rate_peak_num=nan(N_syn,Rep); ILD50_num=nan(N_syn,Rep);
% fixed "x"/stimuli range is 1-2-3-4-5
if strcmp(BW_type, 'Hz') % tuning curve center is "3", random between 1.5 to 4.5 
    shift_sign=[1.5*ones(N_syn*Rep/2,1);-1.5*ones(N_syn*Rep/2,1)];
elseif strcmp(BW_type, 'ILD') % tuning curve center is "1", random between 0.5 to 1.5
    shift_sign=[0.5*ones(N_syn*Rep/2,1);-0.5*ones(N_syn*Rep/2,1)];
end    
shift_sign=shift_sign(randperm(N_syn*Rep));
shift_sign=reshape(shift_sign,[N_syn,Rep]);
for syn = 1 : N_syn
    disp(['synapse number=',num2str(Syn_num_In(syn))])
    for r = 1 : Rep
        cf_shift=shift_sign(syn,r)*rand(1);
        [~, BW_num(syn,r), rate_peak_num(syn,r), ILD50_num(syn,r)]...
            =BW_SYNnum_unmatch(Syn_num_In(syn),BW_type,cf_shift);
    end
end    
%% strength of inhibitory synapse changes compute
Syn_amp_In=(1:10)/10; N_amp=numel(Syn_amp_In);
BW_amp=nan(N_amp,Rep); rate_peak_amp=nan(N_amp,Rep); ILD50_amp=nan(N_amp,Rep);
if strcmp(BW_type, 'Hz')
    shift_sign=[1.5*ones(N_amp*Rep/2,1);-1.5*ones(N_amp*Rep/2,1)];
elseif strcmp(BW_type, 'ILD')
    shift_sign=[0.5*ones(N_amp*Rep/2,1);-0.5*ones(N_amp*Rep/2,1)];
end 
shift_sign=shift_sign(randperm(N_amp*Rep));
shift_sign=reshape(shift_sign,[N_amp,Rep]);
for syn = 1 : N_amp
    disp(['synapse strength=',num2str(Syn_amp_In(syn)*100),'%'])
    for r = 1 : Rep
        cf_shift=shift_sign(syn,r)*rand(1);
        [~, BW_amp(syn,r), rate_peak_amp(syn,r), ILD50_amp(syn,r)]...
            =BW_SYNamp_unmatch(Syn_amp_In(syn),BW_type,cf_shift);
    end
end    
%% bandwidth changes plot & p-values; Figure 5g-Top
figure;
errorbar(1:N_amp, mean(BW_amp,2), std(BW_amp,[],2), 'Marker','.','MarkerSize',20);hold on
errorbar(1:N_syn, mean(BW_num,2), std(BW_num,[],2), 'Marker','.','MarkerSize',20)
legend({'In. Syn. strength','In. Syn. number'})
xlabel('Inhibitory synapses strength/number');ylabel('Half-peak bandwidth (%)')
xlim([0 N_amp+1]); 
if strcmp(BW_type, 'Hz'); ylim([5 30]); elseif strcmp(BW_type, 'ILD'); ylim([5 40]);end
xticks(1:N_amp);
xticklabels({'10','20','30','40','50','60','70','80','90','100'})
BW_num_vs_amp_p=nan(N_syn,1);
for syn=1:N_syn
   BW_num_vs_amp_p(syn)=ranksum(BW_num(syn,:),BW_amp(syn,:));
end    
[BW_num_vs_amp_p, ~]=bonf_holm(BW_num_vs_amp_p,0.05);
%% firing rate changes plot & p-values; Figure 5f
figure;
errorbar(1:N_amp, mean(rate_peak_amp,2), std(rate_peak_amp,[],2),'Marker','.','MarkerSize',20);hold on
errorbar(1:N_syn, mean(rate_peak_num,2), std(rate_peak_num,[],2),'Marker','.','MarkerSize',20)
legend({'In. Syn. strength','In. Syn. number'})
xlabel('Inhibitory synapses strength/number');ylabel('Firing rate (spikes/s)')
xlim([0 N_amp+1]); ylim([0 20])
xticks(1:N_amp);
xticklabels({'10','20','30','40','50','60','70','80','90','100'})
rates_num_vs_amp_p=nan(N_syn,1);
for syn=1:N_syn
   rates_num_vs_amp_p(syn)=ranksum(rate_peak_num(syn,:),rate_peak_amp(syn,:));
end    
[rates_num_vs_amp_p, h]=bonf_holm(rates_num_vs_amp_p,0.05);
%% ILD50% changes plot & p-values; Figure 5g-Bottom
if strcmp(BW_type, 'ILD')
figure;
errorbar(1:N_amp, mean(ILD50_amp,2), std(ILD50_amp,[],2),'Marker','.','MarkerSize',20);hold on
errorbar(1:N_syn, mean(ILD50_num,2), std(ILD50_num,[],2),'Marker','.','MarkerSize',20)
legend({'In. Syn. strength','In. Syn. number'})
xlabel('Inhibitory synapses strength/number');ylabel('50% peak ILD')
xlim([0 N_amp+1]); %ylim([0 30])
xticks(1:N_amp);
xticklabels({'10','20','30','40','50','60','70','80','90','100'})
ILD50_num_vs_amp_p=nan(N_syn,1);
for syn=1:N_syn
   ILD50_num_vs_amp_p(syn)=ranksum(ILD50_num(syn,:),ILD50_amp(syn,:));
end    
[ILD50_num_vs_amp_p, h]=bonf_holm(ILD50_num_vs_amp_p,0.05);
end
%% save the data
unmatch_Hz.BW_num=BW_num;
unmatch_Hz.BW_amp=BW_amp;
unmatch_Hz.BW_num_vs_amp_p=BW_num_vs_amp_p; 

% unmatch_Hz.ILD50_num=ILD50_num;
% unmatch_Hz.ILD50_amp=ILD50_amp;
% unmatch_Hz.ILD50_num_vs_amp_p=ILD50_num_vs_amp_p;

unmatch_Hz.rate_peak_num=rate_peak_num;
unmatch_Hz.rate_peak_amp=rate_peak_amp;
unmatch_Hz.rates_num_vs_amp_p=rates_num_vs_amp_p;