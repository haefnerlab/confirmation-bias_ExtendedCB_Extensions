clear all;close all;clc;
boots = 500;
bins = 15;
hpr_ridge = logspace(-1, 5, 7);
hpr_ar1 = 0.0;
hpr_curvature = logspace(-1, 5, 7);
standardize = 0;
folds = 10;
performance_bnd = 0.7;
thresh_boots = 10;
num_frames = 10;
cases = 2;
time(1) = 1/12;
time(2) = 1/12;
expt_type = 2;
dir = 'RawDataFoveaPeripheryStimulus';
subjects = {...
    'bpgFinaltest-subject04', 'bpgFinaltestBiggerStimulus-subject04';%esther
    'bpgFinaltest-subject12', 'bpgFinaltestBiggerStimulus-subject08';%yihe
    'bpgFinaltest-subject08', 'bpgFinaltestBiggerStimulus-subject05';%anabel
    'bpgFinaltest-subject05', 'bpgFinaltestBiggerStimulus-subject03';%maxxwell
    'bpgFinaltest-subject09', 'bpgFinaltestBiggerStimulus-subject01';%anya
    'bpgFinaltest-subject02', 'bpgFinaltestBiggerStimulus-subject02';%brianna
    'bpgFinaltest-subject06', 'bpgFinaltestBiggerStimulus-subject06';%camryn
    'bpgFinaltest-subject10', 'bpgFinaltestBiggerStimulus-subject07'%julie
    };

subjectsID_fovea = {...
    'bpgFinaltest-subject04';%esther
    'bpgFinaltest-subject12';%yihe
    'bpgFinaltest-subject08';%anabel
    'bpgFinaltest-subject05';%maxxwell
    'bpgFinaltest-subject09';%anya
    'bpgFinaltest-subject02';%brianna
    'bpgFinaltest-subject06';%camryn
    'bpgFinaltest-subject10'%julie
    };

subjectsID_periphery = {...
    'bpgFinaltestBiggerStimulus-subject04';%esther
    'bpgFinaltestBiggerStimulus-subject08';%yihe
    'bpgFinaltestBiggerStimulus-subject05';%anabel
    'bpgFinaltestBiggerStimulus-subject03';%maxxwell
    'bpgFinaltestBiggerStimulus-subject01';%anya
    'bpgFinaltestBiggerStimulus-subject02';%brianna
    'bpgFinaltestBiggerStimulus-subject06';%camryn
    'bpgFinaltestBiggerStimulus-subject07'%julie
    };

disp('Preloading big stimulus statistics for faster post computations....');
[num_sub,~] = size(subjects);
big_matched = 1;
preload_use;

disp('Starting to find best hyperparameters for FOVEAL data across subjects....');
[best_hprs_fovea] = CustomRegression.combined_hprs_search_foveal_stimulus(subjectsID_fovea, expt_type, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds, dir);
disp('Starting to find best hyperparameters for PERIPHERAL data across subjects....');
[best_hprs_periphery] = CustomRegression.combined_hprs_search_peripheral_stimulus(subjectsID_periphery, correct_signals_pre, faulty_signals_pre, ratios, contrasts, expt_type, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds, dir);


for i=1:(num_sub)
    tic;
    figure(i);
    for j=1:cases
        if j==1
            disp(['Starting analysis for FOVEAL data of Subject ' num2str(i) ' ...']);
            best_hprs = best_hprs_fovea;
            [params_boot,~,~,sobl,~,~,abbl_exp,~,~,trials,bin_centers,~,~,data,frame_signals,sobl_time_locked,~,~,log_bernoulli_sub{i,j}] = run_analysis_noise_only(subjects{i,j},expt_type,time(j),boots,best_hprs,standardize,dir);
        else
            disp(['Starting analysis for PERIPHERAL data of Subject ' num2str(i) ' ...']);
            best_hprs = best_hprs_periphery;
            [params_boot,sobl,abbl_exp,trials,bin_centers,~,~,data,frame_signals,sobl_time_locked,log_bernoulli_sub{i,j}] = run_analysis_noise_only_big_stim(subjects{i,j},expt_type,time(j),boots,best_hprs,correct_signals_pre,faulty_signals_pre,ratios,contrasts,standardize,dir);
        end
        %         alpha(i,j,:) = [prctile(params_boot(:, end).^2, 50) std(params_boot(:, end).^2)];
        alpha(i,j,:) = [prctile(1e-4+(1-1e-4) * sigmoid(params_boot(:,end)), 50) std(1e-4+(1-1e-4) * sigmoid(params_boot(:,end)))/sqrt(size(params_boot,1))];
        bias(i,j) = prctile(params_boot(:, end-1), 50);
        temporal_kernel(i,j,:) = prctile(params_boot(:, 1:num_frames), 50);
        norm_temporal_kernel(i,j,:) = temporal_kernel(i,j,:)/mean(temporal_kernel(i,j,:));%prctile(params_boot_norm(:, 1:num_frames), 50);
        lo_temporal_kernel(i,j,:) = prctile(params_boot(:, 1:num_frames), 50) - prctile(params_boot(:, 1:num_frames), 16);
        hi_temporal_kernel(i,j,:) = prctile(params_boot(:, 1:num_frames), 84) - prctile(params_boot(:, 1:num_frames), 50);
        num_trials(i,j) = trials;
        all_exp(i,j,:,:) = abbl_exp;
        beta_all(i,j,:) = abbl_exp(:,2);
        beta(i,j) = prctile(squeeze(all_exp(i,j,:,2)),50);
        all_linear(i,j,:,:) = sobl;
        norm_all_linear(i,j,:,:) = [sobl(:,1)/mean(temporal_kernel(i,j,:)) sobl(:,2)/mean(temporal_kernel(i,j,:)) sobl(:,3) sobl(:,4)];%sobl_norm;
        slope(i,j) = prctile(squeeze(all_linear(i,j,:,1)),50);
        slope_all(i,j,:) = sobl(:,1);
        norm_slope_all(i,j,:) = norm_all_linear(i,j,:,1);
        norm_slope(i,j) = prctile(squeeze(norm_all_linear(i,j,:,1)),50);
        hprs_used(i,j,:) = best_hprs;
        data_sub{i,j} = data;
        rng1 = -50;
        rng2 = 50;
        
        subplot(cases,3,1+3*(j-1))
        plot((1:length(data.noise)), data.noise);
        hold on;
        plot(linspace(1,length(data.noise),100), ones(1,100),'r','LineWidth',1);
        xlabel('Trials');
        ylabel('Noise Level');
        axis('tight');
        thresh(i,j,:) = getBootstrapTheshold(data,frame_signals,performance_bnd,thresh_boots);
        
        subplot(cases,3,3+3*(j-1))
        errorbar(1:num_frames,squeeze(temporal_kernel(i,j,1:num_frames)),squeeze(lo_temporal_kernel(i,j,:)),squeeze(hi_temporal_kernel(i,j,:)),'k','LineWidth',1)
        xlabel('Frames');
        ylabel('Weights');
        axis('tight');
        
        subplot(cases,3,2+3*(j-1))
        bins = 10;
        subject_pm_curve =[];
        uniq_vals = linspace(-0.8,0.8,bins);
        tr_kappa = data.sign_noise;
        noise_signal = data.ideal_frame_signals;
        for tt=1:(length(uniq_vals)-1)
            subj_resp(i,j,tt) = mean(data.choice(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)));
            ntrial_subj(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
        end
        vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
        errorbar(vals,squeeze(subj_resp(i,j,:)),squeeze((subj_resp(i,j,:)).*(1-subj_resp(i,j,:))./sqrt(ntrial_subj(i,j,:))),'o');
        subject_pm_curve = (1./(1+exp(-(noise_signal*squeeze(temporal_kernel(i,j,:))+bias(i,j)))))*( 1-(alpha(i,j,1)))+(alpha(i,j,1)/2);
        for tt=1:(length(uniq_vals)-1)
            subj_resp_pred(i,j,tt) = mean(subject_pm_curve(((tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)))));
            ntrial_subj_pred(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
        end
        hold on;
        plot(vals,squeeze(subj_resp_pred(i,j,:)),'Linewidth',2);
        yline(0.5,'--k');
        xline(0.0,'--k');
        xlabel('Signed Kappa');
        ylabel('Percent chose left');
        xlim([-0.8 0.8])
        ylim([0.0 1.0])
    end
    sgtitle(['Top Row: Foveal Frames and Bottom Row: Peripheral Frames for Subject ' num2str(i)])
    disp(['All analysis complete for Subject ' num2str(i) ' !!!!']);
    toc;
    disp('-----------------------------------------------------------------------------------------------------');
end

disp ('Clearing huge preload files....')
clear_preload_huge_files;
%%
bin_num = 15;
figure();
llo_mn = [];
choice_mn = [];
err_ch = [];
err_ch_mn = [];
llo_mean = [];
choice_mean = [];
for sub=1:num_sub
    subplot(2,5,sub)
    %     llo = llo( ~any( isnan( llo ) | isinf( llo ), 2 ),: );
    [sorted_llo,order_llo] = sort(log_bernoulli_sub{sub,1});
    choice_used = data_sub{sub,1}.choice(order_llo);
    bin_edges = linspace(min(sorted_llo),max(sorted_llo),bin_num+1);
    for bn=1:length(bin_edges)-1
        llo_mean(bn) = mean(sorted_llo(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        choice_mean(bn) = mean(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        err_ch(bn) = sqrt(var(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)))/length(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1))));
    end
    [llo_mn, ord] = sort(llo_mean);
    choice_mn = choice_mean(ord);
    err_ch_mn = err_ch(ord);
    errorbar(llo_mn,choice_mn,err_ch_mn,'-or','Linewidth',2);
    %     shadedErrorBar(llo_mn,choice_mn,err_ch_mn);
    if sub==1 || sub==6
        ylabel('Probability chose vertical','Fontsize',20);
    end
    if sub==3 || sub==8
        xlabel('Log likelihood odds','Fontsize',20);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    xlim([-1*max(abs(llo_mn)) max(abs(llo_mn))]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.0,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Check how good weights predict foveal stimulus behavior','fontsize',30);

bin_num = 15;
figure();
llo_mn = [];
choice_mn = [];
err_ch = [];
err_ch_mn = [];
llo_mean = [];
choice_mean = [];
for sub=1:num_sub
    subplot(2,5,sub)
    %     llo = llo( ~any( isnan( llo ) | isinf( llo ), 2 ),: );
    [sorted_llo,order_llo] = sort(log_bernoulli_sub{sub,2});
    choice_used =  data_sub{sub,2}.choice(order_llo);
    bin_edges = linspace(min(sorted_llo),max(sorted_llo),bin_num+1);
    for bn=1:length(bin_edges)-1
        llo_mean(bn) = mean(sorted_llo(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        choice_mean(bn) = mean(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
        err_ch(bn) = sqrt(var(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)))/length(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1))));
    end
    [llo_mn, ord] = sort(llo_mean);
    choice_mn = choice_mean(ord);
    err_ch_mn = err_ch(ord);
    errorbar(llo_mn,choice_mn,err_ch_mn,'-or','Linewidth',2);
    %     shadedErrorBar(llo_mn,choice_mn,err_ch_mn);
    if sub==1 || sub==6
        ylabel('Probability chose vertical','Fontsize',20);
    end
    if sub==3 || sub==8
        xlabel('Log likelihood odds','Fontsize',20);
    end
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.0 1]);
    xlim([-1*max(abs(llo_mn)) max(abs(llo_mn))]);
    yline(0.5,'-k','linewidth',2);
    hold on;
    xline(0.0,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)], 'fontsize',20)
end
sgtitle('Check how good weights predict peripheral stimulus behavior','fontsize',30);

%%

figure()
for i=1:(num_sub)
    subplot(2,5,i);
    errorbar(1:num_frames,squeeze(temporal_kernel(i,1,1:num_frames)),squeeze(lo_temporal_kernel(i,1,:)),squeeze(hi_temporal_kernel(i,1,:)),'b','LineWidth',2);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    errorbar(1:num_frames,squeeze(temporal_kernel(i,2,1:num_frames)),squeeze(lo_temporal_kernel(i,2,:)),squeeze(hi_temporal_kernel(i,2,:)),'g','LineWidth',2);
    legend({['foveal ' '(' num2str(num_trials(i,1)) ' trials)'],['peripheral ' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Foveal Vs Peripheral Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Temporal weights for all subjects','fontsize',30);

figure()
for i=1:(num_sub)
    subplot(2,5,i);
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'b','LineWidth',2);
    xlabel('Frames');
    ylabel('Norm Weights');
    hold('on');
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'g','LineWidth',2);
    legend({['foveal ' '(' num2str(num_trials(i,1)) ' trials)'],['peripheral ' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Foveal Vs Peripheral Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Normalized temporal weights for all subjects','fontsize',30);

figure()
for i=1:(num_sub)
    subplot(2,5,i);
    plot(1:num_frames,prctile(squeeze(all_exp(i,1,:,1)),50) * exp(prctile(squeeze(all_exp(i,1,:,2)),50) * (0:num_frames-1)),'-ob','LineWidth',2);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(all_exp(i,2,:,1)),50) * exp(prctile(squeeze(all_exp(i,2,:,2)),50) * (0:num_frames-1)),'-og','LineWidth',2);
    legend({['foveal ' '(' num2str(num_trials(i,1)) ' trials)'],['peripheral ' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Foveal Vs Peripheral Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Exponential temporal weights for all subjects','fontsize',30);


figure()
for i=1:(num_sub)
    subplot(2,5,i);
    plot(1:num_frames,prctile(squeeze(all_linear(i,1,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,1,:,1)),50),'-ob','LineWidth',2);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(all_linear(i,2,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,2,:,1)),50),'-og','LineWidth',2);
    legend({['foveal ' '(' num2str(num_trials(i,1)) ' trials)'],['peripheral ' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Foveal Vs Peripheral Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Linear temporal weights for all subjects','fontsize',30);

%%
f = figure();
set(f,'defaultLegendAutoUpdate','off');
subplot(2,2,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(temporal_kernel(i,1,1:num_frames)),'b','LineWidth',0.2);
    cs1_sub(i,:) = squeeze(temporal_kernel(i,1,1:num_frames));
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,squeeze(temporal_kernel(i,2,1:num_frames)),'g','LineWidth',0.2);
    cs2_sub(i,:) = squeeze(temporal_kernel(i,2,1:num_frames));
    legend({['foveal trials'],['peripheral trials']});
    title(['Foveal Vs Peripheral Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-ob','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-og','LineWidth',4);
title('Temporal Weights');

subplot(2,2,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'b','LineWidth',0.2);
    cs1_sub(i,:) = squeeze(norm_temporal_kernel(i,1,1:num_frames));
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'g','LineWidth',0.2);
    cs2_sub(i,:) = squeeze(norm_temporal_kernel(i,2,1:num_frames));
    legend({['foveal trials'],['peripheral trials']});
    title(['Foveal Vs Peripheral Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-ob','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-og','LineWidth',4);
title('Normalized Temporal Weights');

subplot(2,2,4)
for i=1:(num_sub)
    plot(1:num_frames,prctile(squeeze(all_linear(i,1,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,1,:,1)),50),'b','LineWidth',0.2);
    cs1_sub(i,:) = prctile(squeeze(all_linear(i,1,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,1,:,1)),50);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(all_linear(i,2,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,2,:,1)),50),'g','LineWidth',0.2);
    cs2_sub(i,:) = prctile(squeeze(all_linear(i,2,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,2,:,1)),50);
    legend({['foveal trials'],['peripheral trials']});
    title(['Foveal Vs Peripheral Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-ob','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-og','LineWidth',4);
title('Linear Weights');

subplot(2,2,3);
for i=1:(num_sub)
    plot(1:num_frames,prctile(squeeze(all_exp(i,1,:,1)),50) * exp(prctile(squeeze(all_exp(i,1,:,2)),50) * (0:num_frames-1)),'b','LineWidth',0.2);
    cs1_sub(i,:) = prctile(squeeze(all_exp(i,1,:,1)),50) * exp(prctile(squeeze(all_exp(i,1,:,2)),50) * (0:num_frames-1));
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(all_exp(i,2,:,1)),50) * exp(prctile(squeeze(all_exp(i,2,:,2)),50) * (0:num_frames-1)),'g','LineWidth',0.2);
    cs2_sub(i,:) = prctile(squeeze(all_exp(i,2,:,1)),50) * exp(prctile(squeeze(all_exp(i,2,:,2)),50) * (0:num_frames-1));
    legend({['foveal trials'],['peripheral trials']});
    title(['Foveal Vs Peripheral Stimulus for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-ob','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-og','LineWidth',4);
title('Exponential Weights');

sgtitle('All types of temporal weights fit across subjects','fontsize',30);

%%
figure(num_sub+3)
ax1 = subplot(2,5,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(temporal_kernel(i,1,1:num_frames)),'b');
    hold('on');
end
hold('on');
axis('tight')
xlabel('Frames');
ylabel('Weights');
hold('on');
plot(1:num_frames,mean(squeeze((temporal_kernel(:,1,1:num_frames))),1),'b','LineWidth',2);
title('Foveal Stimulus Weights')

ax2 = subplot(2,5,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(temporal_kernel(i,2,1:num_frames)),'g');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Weights');
hold('on');
plot(1:num_frames,mean(squeeze((temporal_kernel(:,2,1:num_frames))),1),'g','LineWidth',2);
title('Peripheral Stimulus Weights')
hold('on');
axis('tight')

ax3 = subplot(2,5,3);
for i=1:(num_sub)
    plot((1:2),slope(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, slope(i,1), 'b');
    hold on;
    scatter(2, slope(i,2), 'g');
    hold on;
end
hold('on');
xticklabels({'Foveal Stimulus','','Peripheral Stimulus'});
xlabel('Foveal Stimulus (blue) to Peripheral Stimulus (green)');
ylabel('Slopes');
hold('on');
title('Slope Comparison')
hold('on');


ax4 = subplot(2,5,4);
scatter(slope(:,1),slope(:,2),'r');
mn_s = min(min(slope));
mx_s = max(max(slope));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
scatter(mean(mean(slope_all(:,1,:))),mean(mean(slope_all(:,2,:))),100,'r','filled');
xlabel('Slopes for Foveal Stimulus ');
ylabel('Slopes for Peripheral Stimulus');
hold('on');
title('Slope Comparison')
hold('on');
axis('tight')

ax5 = subplot(2,5,5);
scatter(beta(:,1),beta(:,2),'r');
mn_b = min(min(beta));
mx_b = max(max(beta));
hold on;
plot(linspace(mn_b,mx_b,10),linspace(mn_b,mx_b,10),'k','LineWidth',2);
hold on;
scatter(mean(mean(beta_all(:,1,:))),mean(mean(beta_all(:,2,:))),100,'r','filled');
xlabel('Beta for Foveal Stimulus ');
ylabel('Beta for Peripheral Stimulus');
hold('on');
title('Beta Comparison')
hold('on');
axis('tight')

linkaxes([ax1,ax2],'y')



ax6 = subplot(2,5,6);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'b');
    hold('on');
end
hold('on');
axis('tight')
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'b','LineWidth',2);
title('Foveal Stimulus Norm Weights')

ax7 = subplot(2,5,7);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'g');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'g','LineWidth',2);
title('Peripheral Stimulus Norm Weights')
hold('on');
axis('tight')


ax8 = subplot(2,5,8);
for i=1:(num_sub)
    plot((1:2),norm_slope(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, norm_slope(i,1), 'b');
    hold on;
    scatter(2, norm_slope(i,2), 'g');
    hold on;
end
hold('on');
xticklabels({'Foveal Stimulus','','Peripheral Stimulus'});
xlabel('Foveal Stimulus (blue) to Periphery Stimulus (green)');
ylabel('Norm Slopes');
hold('on');
title('Norm Slope Comparison')
hold('on');


ax9 = subplot(2,5,9);
scatter(norm_slope(:,1),norm_slope(:,2),'r');
mn_s = min(min(norm_slope));
mx_s = max(max(norm_slope));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
scatter(mean(norm_slope(:,1)),mean(norm_slope(:,2)),100,'r','filled');
xlabel('Norm Slopes for Foveal Stimulus ');
ylabel('Norm Slopes for Peripheral Stimulus');
hold('on');
title('Norm Slope Comparison')
hold('on');
axis('tight')

ax10 = subplot(2,5,10);
scatter(beta(:,1),beta(:,2),'r');
mn_b = min(min(beta));
mx_b = max(max(beta));
hold on;
plot(linspace(mn_b,mx_b,10),linspace(mn_b,mx_b,10),'k','LineWidth',2);
hold on;
scatter(mean(beta(:,1)),mean(beta(:,2)),100,'r','filled');
xlabel('Beta for Foveal Stimulus ');
ylabel('Beta for Peripheral Stimulus');
hold('on');
title('Beta Comparison')
hold('on');
axis('tight')
linkaxes([ax6,ax7],'y')

%% Summary Fig
figure(num_sub+4)
ax1 = subplot(2,4,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'b');
    hold('on');
end
hold('on');
axis('tight')
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'b','LineWidth',2);
title('Foveal Stimulus Norm Weights')

ax2 = subplot(2,4,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'g');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'g','LineWidth',2);
title('Peripheral Stimulus Norm Weights')
hold('on');
axis('tight')
 
ax3 = subplot(2,4,3);
for i=1:(num_sub)
    plot((1:2),norm_slope(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, norm_slope(i,1), 'b');
    hold on;
    scatter(2, norm_slope(i,2), 'g');
    hold on;
end
hold('on');
xticklabels({'Foveal Stimulus','','Peripheral Stimulus'});
ylabel('Slopes');
hold('on');
title('Norm Slope Comparison')
hold('on');

ax4 = subplot(2,4,4);
for i=1:(num_sub)
    plot((1:2),beta(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, beta(i,1), 'b');
    hold on;
    scatter(2, beta(i,2), 'g');
    hold on;
end
hold('on');
xticklabels({'Foveal Stimulus','','Peripheral Stimulus'});
ylabel('Beta');
hold('on');
title('Beta Comparison')
hold('on');

ax5 = subplot(2,4,5);
scatter(norm_slope(:,1),norm_slope(:,2),'k');
mn_s = min(min(norm_slope));
mx_s = max(max(norm_slope));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
scatter(mean(norm_slope(:,1)),mean(norm_slope(:,2)),100,'k','filled');
xlabel('Norm Slopes for Foveal Stimulus ');
ylabel('Norm Slopes for Peripheral Stimulus');
hold('on');
title('Norm Slope Comparison')
hold('on');
axis('tight')

ax6 = subplot(2,4,6);
scatter(beta(:,1),beta(:,2),'k');
mn_b = min(min(beta));
mx_b = max(max(beta));
hold on;
plot(linspace(mn_b,mx_b,10),linspace(mn_b,mx_b,10),'k','LineWidth',2);
hold on;
scatter(mean(beta(:,1)),mean(beta(:,2)),100,'k','filled');
xlabel('Beta for Foveal Stimulus ');
ylabel('Beta for Peripheral Stimulus');
hold('on');
title('Beta Comparison')
hold('on');
axis('tight')

linkaxes([ax1,ax2],'y')

subplot(2,4,7)
for i=1:num_sub
    scatter(squeeze(norm_slope_all(i,1,:)),squeeze(norm_slope_all(i,2,:)),'MarkerFaceColor',rand(1,3));
    hold on;
end
hold on;
for i=1:num_sub
    scatter(mean(squeeze(norm_slope(i,1))),mean(squeeze(norm_slope(i,2))),100,'k','filled');
    hold on;
end
nmn_s_all = min(min(min(norm_slope_all)));
nmx_s_all = max(max(max(norm_slope_all)));
hold on;
plot(linspace(nmn_s_all,nmx_s_all,10),linspace(nmn_s_all,nmx_s_all,10),'k','LineWidth',2);

xlabel('Norm Slope for Foveal Stimulus ');
ylabel('Norm Slope for Peripheral Stimulus');
hold on;
scatter(mean(mean(norm_slope(:,1))),mean(mean(norm_slope(:,2))),100,'r','filled')
axis tight;
title('Norm Slopes')

subplot(2,4,8)
for i=1:num_sub
    scatter(squeeze(beta_all(i,1,:)),squeeze(beta_all(i,2,:)),[],'MarkerFaceColor',rand(1,3));
    hold on;
end
hold on;
for i=1:num_sub
    scatter(mean(squeeze(beta_all(i,1,:))),mean(squeeze(beta_all(i,2,:))),100,'k','filled');
    hold on;
end
mn_b_all = min(min(min(beta_all)));
mx_b_all = max(max(max(beta_all)));
hold on;
plot(linspace(mn_b_all,mx_b_all,10),linspace(mn_b_all,mx_b_all,10),'k','LineWidth',2);
hold on;
scatter(mean(mean(beta_all(:,1,:))),mean(mean(beta_all(:,2,:))),100,'r','filled')

xlabel('Beta for Foveal Stimulus ');
ylabel('Beta for Peripheral Stimulus');
axis tight;
title('Beta');
%% Precise Summary Figure

figure(num_sub+5)
ax1 = subplot(2,2,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'b');
    hold('on');
end
hold('on');
axis('tight');
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'b','LineWidth',2);
title('Foveal Stimulus Norm Weights')

ax2 = subplot(2,2,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'g');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'g','LineWidth',2);
title('Peripheral Stimulus Norm Weights')
hold('on');
axis('tight');

linkaxes([ax1,ax2, ax3],'y')

ax5 = subplot(2,2,3);
scatter(norm_slope(:,1),norm_slope(:,2),'k','filled');
for sb=1:num_sub
    mn1(sb) = mean(squeeze(norm_slope_all(sb,1,:)));
    mn2(sb) = mean(squeeze(norm_slope_all(sb,2,:)));
    v1(sb) = (std(squeeze(norm_slope_all(sb,1,:)))/sqrt(size(norm_slope_all,3)))^2;%var(squeeze(norm_slope_all(sb,1,:)));
    v2(sb) = (std(squeeze(norm_slope_all(sb,2,:)))/sqrt(size(norm_slope_all,3)))^2;%var(squeeze(norm_slope_all(sb,2,:)));
    l1(sb) = prctile(squeeze(norm_slope_all(sb,1,:)),50)-prctile(squeeze(norm_slope_all(sb,1,:)),16);
    h1(sb) = prctile(squeeze(norm_slope_all(sb,1,:)),84)-prctile(squeeze(norm_slope_all(sb,1,:)),50);
    l2(sb) = prctile(squeeze(norm_slope_all(sb,2,:)),50)-prctile(squeeze(norm_slope_all(sb,2,:)),16);
    h2(sb) = prctile(squeeze(norm_slope_all(sb,2,:)),84)-prctile(squeeze(norm_slope_all(sb,2,:)),50);
end
hold on;
eb(1) = errorbar(norm_slope(:,1),norm_slope(:,2),l1,h1, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(norm_slope(:,1),norm_slope(:,2),l2,h2, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'k', 'LineWidth', 2)
mn_s = min(min(norm_slope));
mx_s = max(max(norm_slope));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
mn_slp1 = ((1./v1)./sum(1./v1)) .* norm_slope(:,1)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* norm_slope(:,2)';
scatter(sum(mn_slp1),sum(mn_slp2),200,'r','filled');
xlabel('Norm Slopes for Foveal Stimulus ');
ylabel('Norm Slopes for Peripheral Stimulus');
hold('on');
title('Norm Slope Comparison')
axis('tight');

ax6 = subplot(2,2,4);
scatter(beta(:,1),beta(:,2),'k','filled');
for sb=1:num_sub
    mn1(sb) = mean(squeeze(beta_all(sb,1,:)));
    mn2(sb) = mean(squeeze(beta_all(sb,2,:)));
    v1(sb) = (std(squeeze(beta_all(sb,1,:)))/sqrt(size(beta_all,3)))^2;%var(squeeze(norm_slope_all(sb,1,:)));
    v2(sb) = (std(squeeze(beta_all(sb,2,:)))/sqrt(size(beta_all,3)))^2;%var(squeeze(norm_slope_all(sb,2,:)));
    l1(sb) = prctile(squeeze(beta_all(sb,1,:)),50)-prctile(squeeze(beta_all(sb,1,:)),16);
    h1(sb) = prctile(squeeze(beta_all(sb,1,:)),84)-prctile(squeeze(beta_all(sb,1,:)),50);
    l2(sb) = prctile(squeeze(beta_all(sb,2,:)),50)-prctile(squeeze(beta_all(sb,2,:)),16);
    h2(sb) = prctile(squeeze(beta_all(sb,2,:)),84)-prctile(squeeze(beta_all(sb,2,:)),50);
end
hold on;
eb(1) = errorbar(beta(:,1),beta(:,2),l1,h1, 'horizontal', 'LineStyle', 'none');
hold on;
eb(2) = errorbar(beta(:,1),beta(:,2),l2,h2, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'k', 'LineWidth', 2)
mn_b = min(min(beta));
mx_b = max(max(beta));
hold on;
plot(linspace(mn_b,mx_b,10),linspace(mn_b,mx_b,10),'k','LineWidth',2);
hold on;
mn_slp1 = ((1./v1)./sum(1./v1)) .* beta(:,1)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* beta(:,2)';
scatter(sum(mn_slp1),sum(mn_slp2),200,'r','filled');
xlabel('Beta for Foveal Stimulus ');
ylabel('Beta for Peripheral Stimulus');
hold('on');
title('Beta Comparison')
hold('on');
axis('tight');

%% CCN Abstract Figure
figure()
axis image
subaxis(2,3,1, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'b');
    hold on;
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'g');
    hold on;
end
hold('on');
xlabel('Frames','FontSize',20);
ylabel('Weights','FontSize',20);
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'-ob','LineWidth',2);
hold on;
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'-og','LineWidth',2);
hold('on');
hold on;
plot(1:num_frames,zeros(1,num_frames),'k','LineWidth',2);
legend ({['Foveal Stim'] ['Peripheral Stim']},'Box','off','Fontsize',10)
ylim([-5 5]);
xlim([1 10]);
xticks([1 : 1:num_frames]);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subaxis(2,3,2, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
for sb=1:num_sub
    mn1(sb) = mean(squeeze(thresh(sb,1,:)));
    mn2(sb) = mean(squeeze(thresh(sb,2,:)));
    v1(sb) = (std(squeeze(thresh(sb,1,:)))/sqrt(size(thresh,3)))^2;%var(thresh(sb,1,:));
    v2(sb) = (std(squeeze(thresh(sb,2,:)))/sqrt(size(thresh,3)))^2;%var(thresh(sb,2,:));
    l1(sb) = prctile(squeeze(thresh(sb,1,:)),50) - prctile(squeeze(thresh(sb,1,:)),16);
    h1(sb) = prctile(squeeze(thresh(sb,1,:)),84) - prctile(squeeze(thresh(sb,1,:)),50);
    l2(sb) = prctile(squeeze(thresh(sb,2,:)),50) - prctile(squeeze(thresh(sb,2,:)),16);
    h2(sb) = prctile(squeeze(thresh(sb,2,:)),84) - prctile(squeeze(thresh(sb,2,:)),50);
end
hold on;
eb(1) = errorbar(mn1,mn2,l1,h1, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(mn1,mn2,l2,h2, 'vertical', 'LineStyle', 'none');
set(eb(1), 'color', 'b', 'LineWidth', 2)
set(eb(2), 'color', 'g', 'LineWidth', 2)
hold on;
scatter(mn1,mn2,50,'k','filled');
hold on;
mn = min(min(min(thresh)));
mx = max(max(max(thresh)));
plot(linspace(mn,mx,10),linspace(mn,mx,10),'k','LineWidth',2);
hold on;
mn_slp1 = ((1./v1)./sum(1./v1)) .* mean(squeeze(thresh(:,1,:)),2)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* mean(squeeze(thresh(:,2,:)),2)';
scatter(sum(mn_slp1),sum(mn_slp2),200,'r','filled');
xlabel('Thresh \kappa for Fovea Stim','FontSize',20);
ylabel('Thresh \kappa for Periphery Stim','FontSize',20);
hold('on');
hold('on');
axis('tight');
get(gca,'XTickLabel');
set(gca,'fontsize',18)
get(gca,'YTickLabel');
set(gca,'fontsize',18)
hold on;
ax = gca;
ax.LineWidth = 2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subaxis(2,3,3, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
scatter(norm_slope(:,1),norm_slope(:,2),50,'k','filled');
for sb=1:num_sub
    mn1(sb) = mean(squeeze(norm_slope_all(sb,1,:)));
    mn2(sb) = mean(squeeze(norm_slope_all(sb,2,:)));
    v1(sb) = (std(squeeze(norm_slope_all(sb,1,:)))/sqrt(size(norm_slope_all,3)))^2;%var(squeeze(norm_slope_all(sb,1,:)));
    v2(sb) = (std(squeeze(norm_slope_all(sb,2,:)))/sqrt(size(norm_slope_all,3)))^2;%var(squeeze(norm_slope_all(sb,2,:)));
    l1(sb) = prctile(squeeze(norm_slope_all(sb,1,:)),50)-prctile(squeeze(norm_slope_all(sb,1,:)),16);
    h1(sb) = prctile(squeeze(norm_slope_all(sb,1,:)),84)-prctile(squeeze(norm_slope_all(sb,1,:)),50);
    l2(sb) = prctile(squeeze(norm_slope_all(sb,2,:)),50)-prctile(squeeze(norm_slope_all(sb,2,:)),16);
    h2(sb) = prctile(squeeze(norm_slope_all(sb,2,:)),84)-prctile(squeeze(norm_slope_all(sb,2,:)),50);
end
hold on;
eb(1) = errorbar(norm_slope(:,1), norm_slope(:,2), l1, h1, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(norm_slope(:,1), norm_slope(:,2), l2, h2, 'vertical', 'LineStyle', 'none');
set(eb(1), 'color', 'b', 'LineWidth', 2)
set(eb(2), 'color', 'g', 'LineWidth', 2)
mn_s = min(min([norm_slope(:,1)-l1 norm_slope(:,2)-l2]));
mx_s = max(max([norm_slope(:,1)+h1 norm_slope(:,2)+h2]));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
scatter(norm_slope(:,1),norm_slope(:,2),50,'k','filled')
mn_slp1 = ((1./v1)./sum(1./v1)) .* norm_slope(:,1)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* norm_slope(:,2)';
scatter(sum(mn_slp1),sum(mn_slp2),200,'r','filled');
xlabel('Slopes for Foveal Stimulus ','FontSize',20);
ylabel('Slopes for Peripheral Stimulus','FontSize',20);
hold('on');
hold('on');
axis('tight')
xticks([-0.75 : 0.25:0.2])
yticks([-0.75 : 0.25:0.2])
a = get(gca,'XTickLabel');
set(gca,'fontsize',18)
b = get(gca,'YTickLabel');
set(gca,'fontsize',18)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subaxis(2,3,6, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
for i=1:(num_sub)
    plot(vals,squeeze(subj_resp_pred(i,1,:)),'b','Linewidth',0.2);
    hold('on');
end
plot(vals,squeeze(mean(subj_resp_pred(:,1,:),1)),'-ob','Linewidth',4);
hold('on');
yline(0.5,'k','Linewidth',1);
xline(0.0,'k','Linewidth',1);
axis('tight')
xlabel('Signed Kappa');
ylabel('Percent chose left for Foveal Stimulus');
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subaxis(2,3,5, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
for i=1:(num_sub)
    plot(vals,squeeze(subj_resp_pred(i,2,:)),'g','Linewidth',0.2);
    hold('on');
end
plot(vals,squeeze(mean(subj_resp_pred(:,2,:),1)),'-og','Linewidth',4);
hold('on');
yline(0.5,'k','Linewidth',1);
xline(0.0,'k','Linewidth',1);
axis('tight')
xlabel('Signed Kappa');
ylabel('Percent chose left for Peripheral Stimulus');
set(gca,'fontsize',18)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
filename_save = strcat('SavedWorkspace/FoveaPeripheryEyeTracked_',date,'.mat');
save(filename_save);
