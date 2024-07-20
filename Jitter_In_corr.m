% jitter the number of inhibitory synapse between Hz and ILD which is correlated
% 1 to 10+randi([-1 +1]*0) indicate Hz and ILD have the same number (e.g. 5) of In synapse
% 1 to 10+randi([-1 +1]*5) indicate Hz could be 5-5=0 and ILD could be 5+5=10
% notice thise is different from totally random number in second situation
% because [-5 +5] is the maximum jitter, and it's only exist in 2/(5*2+1)=18% of trials

% call: BW_SYNnum + fit_pr + fit_scatter
% CCG @ 2023-03-15
clear;clc

%% Cotuned Inhibitory neurons Hz vs ILD or ILD50%
Ns=10; 
Nr=10; 
rand_range=0:1:5; % no jitter+5 jitters of 2/4/6/8/10
Nsft=numel(rand_range);
Syn_num_In=repmat(1:Ns,1, Nr); % range from 1-10/Ns, and repeat "Nr" times
Syn_num_In=Syn_num_In(randperm(Ns*Nr)); % shuffle the number
Rep=numel(Syn_num_In);
BW_Hz_In=nan(Rep,Nsft); %Hz50_In=nan(Rep,Nsft); %rate_peak_Hz=nan(Rep,1); 
BW_ILD_In=nan(Rep,Nsft);ILD50_In=nan(Rep,Nsft);%rate_peak_ILD=nan(Rep,1);
fit_pr_HzILD=nan(2,Nsft);fit_pr_HzILD50=nan(2,Nsft);
syn_Hz_all=nan(Nsft, Rep);
syn_ILD_all=nan(Nsft, Rep);
for ran= 1 : Nsft % 1 to 6
    ran_val=rand_range(ran); % 0(no-jitter) to 5(100%jitter)
for r = 1 : Rep
    if mod(r,Ns)==0; disp(['randi=',num2str(ran_val), ' rep=',num2str(r/Ns)]);end
    syn_Hz=Syn_num_In(r)+randi([-1 1]*ran_val);
    syn_Hz(syn_Hz<0)=0;
    syn_Hz_all(ran, r)=syn_Hz;
    [~, BW_Hz_In(r,ran), ~, ~]=BW_SYNnum(syn_Hz,'Hz');
    syn_ILD=Syn_num_In(r)+randi([-1 1])*ran_val;
    syn_ILD(syn_ILD<0)=0;
    syn_ILD_all(ran, r)=syn_ILD;
    [~, BW_ILD_In(r,ran), ~, ILD50_In(r,ran)]=BW_SYNnum(syn_ILD,'ILD');
end
    [fit_pr_HzILD(1,ran),fit_pr_HzILD(2,ran)]=fit_pr(BW_Hz_In(:,ran), BW_ILD_In(:,ran));
    [fit_pr_HzILD50(1,ran),fit_pr_HzILD50(2,ran)]=fit_pr(BW_Hz_In(:,ran), ILD50_In(:,ran));
end
%% show the jittered synapse number of ILD and Hz tuning inputs
pos=get(0,'ScreenSize'); X_size=pos(3);Y_size=pos(4);
figure('position',[X_size*0.01 Y_size*0.3 1600 200]);
tiledlayout(1, Nsft);
for stim = 1: Nsft
    nexttile
    plot(syn_Hz_all(stim, :)); 
    hold on
    plot(syn_ILD_all(stim, :));
    title(['Jitter=', num2str(rand_range(stim)*20),'%'])
    if stim==Nsft
    legend({'Freq.', 'ILD'})
    end
end
%% Independent/Unjittered Excitatory neurons Hz vs ILD or ILD50%
Ns=10;
Nr=5;
Syn_num_In=repmat(1:Ns,1,Nr); 
Syn_num_In_Hz=Syn_num_In(randperm(Ns*Nr));
Syn_num_In_ILD=Syn_num_In(randperm(Ns*Nr));
Rep=numel(Syn_num_In_Hz);
BW_Hz_Ex=nan(Rep,1); rate_peak_Hz=nan(Rep,1); Hz50_Ex=nan(Rep,1);
BW_ILD_Ex=nan(Rep,1); rate_peak_ILD=nan(Rep,1); ILD50_Ex=nan(Rep,1);
for r = 1 : Rep
    if mod(r,Ns)==0; disp(['rep=',num2str(r/Ns)]);end
    syn_Hz=Syn_num_In_Hz(r);
    [~, BW_Hz_Ex(r), rate_peak_Hz(r), Hz50_Ex(r)]...
        =BW_SYNnum(syn_Hz,'Hz');
    syn_ILD=Syn_num_In_ILD(r);
    [~, BW_ILD_Ex(r), rate_peak_ILD(r), ILD50_Ex(r)]...
        =BW_SYNnum(syn_ILD,'ILD');
end
lm_Ex_HzILD=   fitlm(BW_Hz_Ex, BW_ILD_Ex); 
lm_Ex_HzILD50= fitlm(BW_Hz_Ex, ILD50_Ex);
Ex_HzILD_p=lm_Ex_HzILD.Coefficients.pValue(2);Ex_HzILD_r=sqrt(lm_Ex_HzILD.Rsquared.Ordinary);
Ex_HzILD50_p=lm_Ex_HzILD50.Coefficients.pValue(2);Ex_HzILD50_r=sqrt(lm_Ex_HzILD50.Rsquared.Ordinary);
%% plot the scatters of BW-Hz vs BW-ILD
figure; x_range='fixed';
colors = turbo(7);
for ran= 1 : Nsft
[x,y]=fit_scatter(BW_Hz_In(:,ran)*2/100, BW_ILD_In(:,ran)/100, x_range); 
plot(x,y,'Color',colors(ran+1,:),'LineWidth',2);hold on
text(0.01,0.01+0.05*ran,['p=',num2str(fit_pr_HzILD(1,ran)),' r=',num2str(fit_pr_HzILD(2,ran)),...
    '  5 +/- ',num2str(rand_range(ran)), '(',num2str(100*(rand_range(ran))/5), '%)'])
end
[x,y]=fit_scatter(BW_Hz_Ex*2/100, BW_ILD_Ex/100, x_range); 
plot(x,y,'Color',rgb('black'),'LineWidth',2);
legend({num2str((0:6)')},'Location','southeast')
text(0.01,0.01+0.05*(ran+1),['Excitatory/black p=',num2str(Ex_HzILD_p),' r=-',num2str(Ex_HzILD_r)])
% axis equal
xlim([0 .6]);xlabel('Bandwith of Hz')
ylim([0 .4]);ylabel('Bandwith of ILD') 
%% plot the scatters of BW-Hz vs ILD50
figure;sz=30; % tiledlayout(1,2);
colors = turbo(7); x_range='fixed' ;% 'fixed' or 'free'
for ran= 1 : Nsft
[x,y]=fit_scatter(BW_Hz_In(:,ran)*2/100, ILD50_In(:,ran)/1000, x_range); 
plot(x,y,'Color',colors(ran+1,:),'LineWidth',2);hold on
text(0.01,0.74-0.03*ran,['p=',num2str(fit_pr_HzILD(1,ran)),' r=-',num2str(fit_pr_HzILD(2,ran)),...
    '  5 +/- ',num2str(rand_range(ran)), '(',num2str(100*(rand_range(ran))/5), '%)'])
end
[x,y]=fit_scatter(BW_Hz_Ex*2/100, ILD50_Ex/1000, x_range); 
plot(x,y,'Color',rgb('black'),'LineWidth',2);
text(0.01,0.74-0.03*(ran+1),['Excitatory/black p=',num2str(Ex_HzILD50_p),' r=',num2str(Ex_HzILD50_r)])
% axis equal
xlim([0 .6]);xlabel('Bandwith of Hz')
ylim([0.45 .75]);ylabel('ILD 50% point') 