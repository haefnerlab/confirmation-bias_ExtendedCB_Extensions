clear all; close all; clc;
boots = 500;
bins = 15;
hpr_ridge = 0.0;%logspace(-1, 5, 7);
hpr_ar1 = 0.0;
hpr_curvature = logspace(-1, 5, 7);
num_frames = 10;
cases = 2; % short and long frames, hence, 2 cases
% for time-locked computations
time(1) = 1/24;
time(2) = 1/6;
dir = 'RawDataLongShortFrameDuration';
expt_type = 2;
standardize = 0;
folds = 10;
performance_bnd = 0.7;
thresh_boots = 10;
subjects = {...
    'bpgshorterframes-subject02', 'bpglongerframes-subject02';
    'bpgshorterframes-subject04', 'bpglongerframes-subject04';
    'bpgshorterframes-subject05', 'bpglongerframes-subject06';
    'bpgshorterframes-subject03', 'bpglongerframes-subject07';
    'bpgshorterframes-subject11', 'bpglongerframes-subject09';
    'bpgshorterframes-subject08', 'bpglongerframes-subject12';%colleen
    'bpgshorterframes-subject10', 'bpglongerframes-subject15';
    'bpgshorterframes-subject07', 'bpglongerframes-subject10';
    'bpgshorterframes-subject14', 'bpglongerframes-subject16';
    'bpgshorterframes-subject13', 'bpglongerframes-subject18'%eimi
    };
subjectID_short = {...
    'bpgshorterframes-subject02';...
    'bpgshorterframes-subject04';
    'bpgshorterframes-subject05';
    'bpgshorterframes-subject03';
    'bpgshorterframes-subject11';
    'bpgshorterframes-subject08';%colleen
    'bpgshorterframes-subject10';
    'bpgshorterframes-subject07';
    'bpgshorterframes-subject14';
    'bpgshorterframes-subject13'%eimi
    };
subjectID_long = {...
    'bpglongerframes-subject02';
    'bpglongerframes-subject04';
    'bpglongerframes-subject06';
    'bpglongerframes-subject07';
    'bpglongerframes-subject09';
    'bpglongerframes-subject12';%colleen
    'bpglongerframes-subject15';
    'bpglongerframes-subject10';
    'bpglongerframes-subject16';
    'bpglongerframes-subject18'%eimi
    };

[num_sub,~] = size(subjects);
disp('Starting to find best hyperparameters for SHORT FRAMES data across subjects....');
[best_hprs_short] = CustomRegression.combined_hprs_search_frame_duration(subjectID_short, expt_type, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds, dir);
disp('Starting to find best hyperparameters for LONG FRAMES data across subjects....');
[best_hprs_long] = CustomRegression.combined_hprs_search_frame_duration(subjectID_long, expt_type, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds, dir);

for i=1:(num_sub)
    tic;
    figure(i)
    for j=1:cases
        if j==1
            best_hprs = best_hprs_short;
            disp(['Starting analysis for SHORT FRAMES data of Subject ' num2str(i) ' ...']);
        elseif j==2
            best_hprs = best_hprs_long;
            disp(['Starting analysis for LONG FRAMES data of Subject ' num2str(i) ' ...']);
        end
        [params_boot,params_boot_first_half,params_boot_second_half,...
            sobl,sobl_first_half,sobl_second_half,...
            abbl_exp,abbl_exp_first_half,abbl_exp_second_half,...
            trials,bin_centers,~,~,data,frame_signals,...
            sobl_time_locked,sobl_time_locked_first_half,sobl_time_locked_second_half,...
            log_bernoulli{i,j}] = run_analysis_noise_only(subjects{i,j}, expt_type, time(j), boots, best_hprs, standardize, dir);
        %         alpha(i,j,:) = [prctile(params_boot(:, end).^2, 50) std(params_boot(:, end).^2)];
        alpha(i,j,:) = [prctile(1e-4+(1-1e-4) * sigmoid(params_boot(:,end)), 50) std(1e-4+(1-1e-4) * sigmoid(params_boot(:,end)))/sqrt(size(params_boot,1))];
        bias(i,j) = prctile(params_boot(:, end-1), 50);
        temporal_kernel(i,j,:) = prctile(params_boot(:, 1:num_frames), 50);
        norm_temporal_kernel(i,j,:) = temporal_kernel(i,j,:)/mean(temporal_kernel(i,j,:));
        lo_temporal_kernel(i,j,:) = prctile(params_boot(:, 1:num_frames), 50) - prctile(params_boot(:, 1:num_frames), 16);
        hi_temporal_kernel(i,j,:) = prctile(params_boot(:, 1:num_frames), 84) - prctile(params_boot(:, 1:num_frames), 50);
        num_trials(i,j) = trials;
        all_exp(i,j,:,:) = abbl_exp;
        beta_all(i,j,:) = abbl_exp(:,2);
        beta(i,j) = prctile(squeeze(all_exp(i,j,:,2)),50);
        all_linear(i,j,:,:) = sobl;
        norm_all_linear(i,j,:,:) = [sobl(:,1)/mean(temporal_kernel(i,j,:)) sobl(:,2)/mean(temporal_kernel(i,j,:)) sobl(:,3) sobl(:,4)];
        norm_all_linear_time_locked(i,j,:,:) = [sobl_time_locked(:,1)/mean(temporal_kernel(i,j,:)) sobl_time_locked(:,2)/mean(temporal_kernel(i,j,:)) sobl_time_locked(:,3) sobl_time_locked(:,4)];%sobl_norm;
        all_linear_time_locked(i,j,:,:) = sobl_time_locked;
        slope(i,j) = prctile(squeeze(all_linear(i,j,:,1)),50);
        slope_all(i,j,:) = sobl(:,1);
        norm_slope_all(i,j,:) = norm_all_linear(i,j,:,1);
        norm_slope_all_time_locked(i,j,:) = norm_all_linear_time_locked(i,j,:,1);
        norm_slope(i,j) = prctile(squeeze(norm_all_linear(i,j,:,1)),50);
        norm_slope_time_locked(i,j) = prctile(squeeze(norm_all_linear_time_locked(i,j,:,1)),50);
        hprs_used(i,j,:) = best_hprs;
        data_sub{i,j} = data;
        frame_signals_sub{i,j} = frame_signals;
        
        temporal_kernel_first_half(i,j,:) = prctile(params_boot_first_half(:, 1:num_frames), 50);
        norm_temporal_kernel_first_half(i,j,:) = temporal_kernel_first_half(i,j,:)/mean(temporal_kernel_first_half(i,j,:));
        lo_temporal_kernel_first_half(i,j,:) = prctile(params_boot_first_half(:, 1:num_frames), 50) - prctile(params_boot_first_half(:, 1:num_frames), 16);
        hi_temporal_kernel_first_half(i,j,:) = prctile(params_boot_first_half(:, 1:num_frames), 84) - prctile(params_boot_first_half(:, 1:num_frames), 50);
        
        all_linear_first_half(i,j,:,:) = sobl_first_half;
        norm_all_linear_first_half(i,j,:,:) = [sobl_first_half(:,1)/mean(temporal_kernel_first_half(i,j,:)) sobl_first_half(:,2)/mean(temporal_kernel_first_half(i,j,:)) sobl_first_half(:,3) sobl_first_half(:,4)];
        norm_all_linear_time_locked_first_half(i,j,:,:) = [sobl_time_locked_first_half(:,1)/mean(temporal_kernel_first_half(i,j,:)) sobl_time_locked_first_half(:,2)/mean(temporal_kernel_first_half(i,j,:)) sobl_time_locked_first_half(:,3) sobl_time_locked_first_half(:,4)];%sobl_norm;
        all_linear_time_locked_first_half(i,j,:,:) = sobl_time_locked_first_half;
        slope_first_half(i,j) = prctile(squeeze(all_linear_first_half(i,j,:,1)),50);
        slope_all_first_half(i,j,:) = sobl_first_half(:,1);
        norm_slope_all_first_half(i,j,:) = norm_all_linear_first_half(i,j,:,1);
        norm_slope_all_time_locked_first_half(i,j,:) = norm_all_linear_time_locked_first_half(i,j,:,1);
        norm_slope_first_half(i,j) = prctile(squeeze(norm_all_linear_first_half(i,j,:,1)),50);
        norm_slope_time_locked_first_half(i,j) = prctile(squeeze(norm_all_linear_time_locked_first_half(i,j,:,1)),50);
        
        
        temporal_kernel_second_half(i,j,:) = prctile(params_boot_second_half(:, 1:num_frames), 50);
        norm_temporal_kernel_second_half(i,j,:) = temporal_kernel_second_half(i,j,:)/mean(temporal_kernel_second_half(i,j,:));
        lo_temporal_kernel_second_half(i,j,:) = prctile(params_boot_second_half(:, 1:num_frames), 50) - prctile(params_boot_second_half(:, 1:num_frames), 16);
        hi_temporal_kernel_second_half(i,j,:) = prctile(params_boot_second_half(:, 1:num_frames), 84) - prctile(params_boot_second_half(:, 1:num_frames), 50);
        
        all_linear_second_half(i,j,:,:) = sobl_second_half;
        norm_all_linear_second_half(i,j,:,:) = [sobl_second_half(:,1)/mean(temporal_kernel_second_half(i,j,:)) sobl_second_half(:,2)/mean(temporal_kernel_second_half(i,j,:)) sobl_second_half(:,3) sobl_second_half(:,4)];
        norm_all_linear_time_locked_second_half(i,j,:,:) = [sobl_time_locked_second_half(:,1)/mean(temporal_kernel_second_half(i,j,:)) sobl_time_locked_second_half(:,2)/mean(temporal_kernel_second_half(i,j,:)) sobl_time_locked_second_half(:,3) sobl_time_locked_second_half(:,4)];%sobl_norm;
        all_linear_time_locked_second_half(i,j,:,:) = sobl_time_locked_second_half;
        slope_second_half(i,j) = prctile(squeeze(all_linear_second_half(i,j,:,1)),50);
        slope_all_second_half(i,j,:) = sobl_second_half(:,1);
        norm_slope_all_second_half(i,j,:) = norm_all_linear_second_half(i,j,:,1);
        norm_slope_all_time_locked_second_half(i,j,:) = norm_all_linear_time_locked_second_half(i,j,:,1);
        norm_slope_second_half(i,j) = prctile(squeeze(norm_all_linear_second_half(i,j,:,1)),50);
        norm_slope_time_locked_second_half(i,j) = prctile(squeeze(norm_all_linear_time_locked_second_half(i,j,:,1)),50);
        
        
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
        subject_pm_curve = [];
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
    sgtitle(['Top Row: Short Frames and Bottom Row: Longer Frames for Subject ' num2str(i)],'fontsize',30);
    disp(['All analysis complete for Subject ' num2str(i) ' !!!!']);
    toc;
    disp('-----------------------------------------------------------------------------------------------------');
end

%%
figure()
j=1;
for i=1:num_sub
    subplot(2,5,i)
    bins = 10;
    subject_pm_curve =[];
    uniq_vals = linspace(-0.8,0.8,bins);
    tr_kappa = data_sub{i,j}.sign_noise;
    noise_signal = data_sub{i,j}.ideal_frame_signals;
    for tt=1:(length(uniq_vals)-1)
        subj_resp(i,j,tt) = mean(data_sub{i,j}.choice(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)));
        ntrial_subj(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
    end
    vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
    errorbar(vals,squeeze(subj_resp(i,j,:)),squeeze((subj_resp(i,j,:)).*(1-subj_resp(i,j,:))./sqrt(ntrial_subj(i,j,:))),'ob','Linestyle','none','linewidth',2);
    subject_pm_curve = (1./(1+exp(-(noise_signal*squeeze(temporal_kernel(i,j,:))+bias(i,j)))))*( 1-(alpha(i,j,1)))+(alpha(i,j,1)/2);
    for tt=1:(length(uniq_vals)-1)
        subj_resp_pred(i,j,tt) = mean(subject_pm_curve(((tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)))));
        ntrial_subj_pred(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
    end
    hold on;
    err = sqrt(squeeze(subj_resp_pred(i,j,:)).*(1-squeeze(subj_resp_pred(i,j,:)))./squeeze(ntrial_subj(i,j,:)));
    errorbar(vals,squeeze(subj_resp_pred(i,j,:)),err,'-or','Linewidth',2);
    yline(0.5,'--k');
    xline(0.0,'--k');
    xlabel('Signed Kappa');
    ylabel('Percent chose left');
    xlim([-0.8 0.8])
    ylim([0.0 1.0])
end
sgtitle('Response predicted by fitted weights to real data for short duration stimuli','fontsize',30);

figure()
j=2;
for i=1:num_sub
    subplot(2,5,i)
    bins = 10;
    subject_pm_curve =[];
    uniq_vals = linspace(-0.8,0.8,bins);
    tr_kappa = data_sub{i,j}.sign_noise;
    noise_signal = data_sub{i,j}.ideal_frame_signals;
    for tt=1:(length(uniq_vals)-1)
        subj_resp(i,j,tt) = mean(data_sub{i,j}.choice(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)));
        ntrial_subj(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
    end
    vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
    errorbar(vals,squeeze(subj_resp(i,j,:)),squeeze((subj_resp(i,j,:)).*(1-subj_resp(i,j,:))./sqrt(ntrial_subj(i,j,:))),'ob','Linestyle','none','linewidth',2);
    subject_pm_curve = (1./(1+exp(-(noise_signal*squeeze(temporal_kernel(i,j,:))+bias(i,j)))))*( 1-(alpha(i,j,1)))+(alpha(i,j,1)/2);
    for tt=1:(length(uniq_vals)-1)
        subj_resp_pred(i,j,tt) = mean(subject_pm_curve(((tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)))));
        ntrial_subj_pred(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
    end
    hold on;
    err = sqrt(squeeze(subj_resp_pred(i,j,:)).*(1-squeeze(subj_resp_pred(i,j,:)))./squeeze(ntrial_subj(i,j,:)));
    errorbar(vals,squeeze(subj_resp_pred(i,j,:)),err,'-or','Linewidth',2);
    yline(0.5,'--k');
    xline(0.0,'--k');
    xlabel('Signed Kappa');
    ylabel('Percent chose left');
    xlim([-0.8 0.8])
    ylim([0.0 1.0])
end
sgtitle('Response predicted by fitted weights to real data for long duration stimuli','fontsize',30);

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
    [sorted_llo,order_llo] = sort(log_bernoulli{sub,1});
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
sgtitle('Check how good weights predict short frame behavior','fontsize',30);

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
    [sorted_llo,order_llo] = sort(log_bernoulli{sub,2});
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
sgtitle('Check how good weights predict long frame behavior','fontsize',30);
%%

figure()
for i=1:(num_sub)
    subplot(2,5,i);
    errorbar(1:num_frames,squeeze(temporal_kernel(i,1,1:num_frames)),squeeze(lo_temporal_kernel(i,1,:)),squeeze(hi_temporal_kernel(i,1,:)),'-om','LineWidth',2);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    errorbar(1:num_frames,squeeze(temporal_kernel(i,2,1:num_frames)),squeeze(lo_temporal_kernel(i,2,:)),squeeze(hi_temporal_kernel(i,2,:)),'-oc','LineWidth',2);
    legend({['short' '(' num2str(num_trials(i,1)) ' trials)'],['long' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Short Vs Long Frames for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Temporal weights for all subjects','fontsize',30);

figure()
for i=1:(num_sub)
    subplot(2,5,i);
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'-om','LineWidth',2);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'-oc','LineWidth',2);
    legend({['short' '(' num2str(num_trials(i,1)) ' trials)'],['long' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Short Vs Long Frames for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Normalized temporal weights for all subjects','fontsize',30);

%%
figure()
for i=1:(num_sub)
    subplot(2,5,i);
    plot(1:num_frames,prctile(squeeze(all_exp(i,1,:,1)),50) * exp(prctile(squeeze(all_exp(i,1,:,2)),50) * (0:num_frames-1)),'-om','LineWidth',2);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(all_exp(i,2,:,1)),50) * exp(prctile(squeeze(all_exp(i,2,:,2)),50) * (0:num_frames-1)),'-oc','LineWidth',2);
    legend({['short' '(' num2str(num_trials(i,1)) ' trials)'],['long' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Short Vs Long Frames for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Exponential temporal weights for all subjects','fontsize',30);


figure()
for i=1:(num_sub)
    subplot(2,5,i);
    plot(1:num_frames,prctile(squeeze(all_linear(i,1,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,1,:,1)),50),'-om','LineWidth',2);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(all_linear(i,2,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,2,:,1)),50),'-oc','LineWidth',2);
    legend({['short' '(' num2str(num_trials(i,1)) ' trials)'],['long' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Short Vs Long Frames for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
sgtitle('Linear temporal weights for all subjects','fontsize',30);

%%
f = figure();
set(f,'defaultLegendAutoUpdate','off');
subplot(2,2,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(temporal_kernel(i,1,1:num_frames)),'m','LineWidth',0.2);
    cs1_sub(i,:) = squeeze(temporal_kernel(i,1,1:num_frames));
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,squeeze(temporal_kernel(i,2,1:num_frames)),'c','LineWidth',0.2);
    cs2_sub(i,:) = squeeze(temporal_kernel(i,2,1:num_frames));
    legend({['short trials'],['long trials']});
    title(['Short Vs Long Frames for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-om','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-oc','LineWidth',4);
title('Temporal Weights','fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(2,2,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'m','LineWidth',0.2);
    cs1_sub(i,:) = squeeze(norm_temporal_kernel(i,1,1:num_frames));
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'c','LineWidth',0.2);
    cs2_sub(i,:) = squeeze(norm_temporal_kernel(i,2,1:num_frames));
    legend({['short trials'],['long trials']});
    title(['Short Vs Long Frames for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-om','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-oc','LineWidth',4);
title('Normalized Temporal Weights','fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(2,2,4)
for i=1:(num_sub)
    plot(1:num_frames,prctile(squeeze(all_linear(i,1,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,1,:,1)),50),'m','LineWidth',0.2);
    cs1_sub(i,:) = prctile(squeeze(all_linear(i,1,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,1,:,1)),50);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(all_linear(i,2,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,2,:,1)),50),'c','LineWidth',0.2);
    cs2_sub(i,:) = prctile(squeeze(all_linear(i,2,:,2)),50) + (0:num_frames-1) * prctile(squeeze(all_linear(i,2,:,1)),50);
    legend({['short trials'],['long trials']});
    title(['Short Vs Long Frames for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-om','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-oc','LineWidth',4);
title('Linear Weights','fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subplot(2,2,3);
for i=1:(num_sub)
    plot(1:num_frames,prctile(squeeze(all_exp(i,1,:,1)),50) * exp(prctile(squeeze(all_exp(i,1,:,2)),50) * (0:num_frames-1)),'m','LineWidth',0.2);
    cs1_sub(i,:) = prctile(squeeze(all_exp(i,1,:,1)),50) * exp(prctile(squeeze(all_exp(i,1,:,2)),50) * (0:num_frames-1));
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    plot(1:num_frames,prctile(squeeze(all_exp(i,2,:,1)),50) * exp(prctile(squeeze(all_exp(i,2,:,2)),50) * (0:num_frames-1)),'c','LineWidth',0.2);
    cs2_sub(i,:) = prctile(squeeze(all_exp(i,2,:,1)),50) * exp(prctile(squeeze(all_exp(i,2,:,2)),50) * (0:num_frames-1));
    legend({['short trials'],['long trials']});
    title(['Short Vs Long Frames for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
yline(0.0,'k','linewidth',2);
plot(1:num_frames,mean(cs1_sub,1),'-om','LineWidth',4);
plot(1:num_frames,mean(cs2_sub,1),'-oc','LineWidth',4);
title('Exponential Weights','fontsize',20); 
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

sgtitle('All types of temporal weights fit across subjects','fontsize',30);
%%
figure(num_sub+3)
ax1 = subplot(2,5,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(temporal_kernel(i,1,1:num_frames)),'m');
    hold('on');
end
hold('on');
axis('tight')
xlabel('Frames');
ylabel('Weights');
hold('on');
plot(1:num_frames,mean(squeeze((temporal_kernel(:,1,1:num_frames))),1),'m','LineWidth',2);
title('Short Frames Weights')

ax2 = subplot(2,5,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(temporal_kernel(i,2,1:num_frames)),'c');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Weights');
hold('on');
plot(1:num_frames,mean(squeeze((temporal_kernel(:,2,1:num_frames))),1),'c','LineWidth',2);
title('Long Frames Weights')
hold('on');
axis('tight')

ax3 = subplot(2,5,3);
for i=1:(num_sub)
    plot((1:2),slope(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, slope(i,1), 'm');
    hold on;
    scatter(2, slope(i,2), 'c');
    hold on;
end
hold('on');
xticklabels({'Short Frames','','Long Frames'});
xlabel('Short Frames (blue) to Long Frames (green)');
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
xlabel('Slopes for Short Frames ');
ylabel('Slopes for Long Frames');
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
xlabel('Beta for Short Frames ');
ylabel('Beta for Long Frames');
hold('on');
title('Beta Comparison')
hold('on');
axis('tight')

linkaxes([ax1,ax2],'y')

ax6 = subplot(2,5,6);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'m');
    hold('on');
end
hold('on');
axis('tight')
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'m','LineWidth',2);
title('Short Frames Norm Weights')

ax7 = subplot(2,5,7);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'c');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'c','LineWidth',2);
title('Long Frames Norm Weights')
hold('on');
axis('tight')

ax8 = subplot(2,5,8);
for i=1:(num_sub)
    plot((1:2),norm_slope(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, norm_slope(i,1), 'm');
    hold on;
    scatter(2, norm_slope(i,2), 'c');
    hold on;
end
hold('on');
xticklabels({'Short Frames','','Long Frames'});
xlabel('Short Frames (blue) to Long Frames (green)');
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
scatter(mean(mean(norm_slope_all(:,1,:))),mean(mean(norm_slope_all(:,2,:))),100,'r','filled');
xlabel('Norm Slopes for Short Frames ');
ylabel('Norm Slopes for Long Frames');
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
scatter(mean(mean(beta_all(:,1,:))),mean(mean(beta_all(:,2,:))),100,'r','filled');
xlabel('Beta for Short Frames ');
ylabel('Beta for Long Frames');
hold('on');
title('Beta Comparison')
hold('on');
axis('tight')
linkaxes([ax6,ax7],'y')

%% Summary Fig
figure(num_sub+4)
ax1 = subplot(2,4,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'m');
    hold('on');
end
hold('on');
axis('tight')
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'m','LineWidth',2);
title('Short Frames Norm Weights')

ax2 = subplot(2,4,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'c');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'c','LineWidth',2);
title('Long Frames Norm Weights')
hold('on');
axis('tight')

ax3 = subplot(2,4,3);
for i=1:(num_sub)
    plot((1:2),norm_slope(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, norm_slope(i,1), 'm');
    hold on;
    scatter(2, norm_slope(i,2), 'c');
    hold on;
end
hold('on');
xticklabels({'Short Frames','','Long Frames'});
ylabel('Slopes');
hold('on');
title('Norm Slope Comparison')
hold('on');

ax4 = subplot(2,4,4);
for i=1:(num_sub)
    plot((1:2), beta(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, beta(i,1), 'm');
    hold on;
    scatter(2, beta(i,2), 'c');
    hold on;
end
hold('on');
xticklabels({'Short Frames','','Long Frames'});
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
xlabel('Norm Slopes for Short Frames ');
ylabel('Norm Slopes for Long Frames');
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
xlabel('Beta for Short Frames ');
ylabel('Beta for Long Frames');
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
xlabel('Norm Slope for Short Frames ');
ylabel('Norm Slope for Long Frames');
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
xlabel('Beta for Short Frames ');
ylabel('Beta for Long Frames');
axis tight;
title('Beta');


%% Precise Summary Figure
figure(num_sub+5)
ax1=subplot(2,2,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'m');
    hold('on');
end
hold('on');
axis('tight')
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'m','LineWidth',2);
title('Short Frames Norm Weights')

ax2=subplot(2,2,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'c');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'c','LineWidth',2);
title('Long Frames Norm Weights')
hold('on');
axis('tight')

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
xlabel('Norm Slopes for Short Frames ');
ylabel('Norm Slopes for Long Frames');
hold('on');
title('Norm Slope Comparison')
axis('tight')

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
xlabel('Beta for Short Frames ');
ylabel('Beta for Long Frames');
hold('on');
title('Beta Comparison')
hold('on');
axis('tight')


%% CCN Abstract Figure
figure(num_sub+6)
time(1) = 1/24;
time(2) = 1/6;
axis image
subaxis(2,4,1, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'m');
    hold on;
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'c');
    hold on;
end
hold('on');
xlabel('Frames','FontSize',20);
ylabel('Weights','FontSize',20);
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'-om','LineWidth',2);
hold on;
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'-oc','LineWidth',2);
hold('on');
hold on;
plot(1:num_frames,zeros(1,num_frames),'k','LineWidth',2);
legend ({['Short frame'] ['Long frame']},'Box','off','Fontsize',10);
axis('tight');
xticks([1 : 1:num_frames])
a = get(gca,'XTickLabel');
set(gca,'fontsize',18)
b = get(gca,'YTickLabel');
set(gca,'fontsize',18)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subaxis(2,4,2, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
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
set(eb(1), 'color', 'm', 'LineWidth', 2)
set(eb(2), 'color', 'c', 'LineWidth', 2)
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
xlabel('Thresh \kappa for Short Frame','FontSize',20);
ylabel('Thresh \kappa for Long Frame','FontSize',20);
hold('on');
hold('on');
axis('tight');
get(gca,'XTickLabel');
set(gca,'fontsize',18)
get(gca,'YTickLabel');
set(gca,'fontsize',18)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subaxis(2,4,3, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
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
set(eb(1), 'color', 'm', 'LineWidth', 2)
set(eb(2), 'color', 'c', 'LineWidth', 2)
mn_s = min(min([norm_slope(:,1)-l1 norm_slope(:,2)-l2]));
mx_s = max(max([norm_slope(:,1)+h1 norm_slope(:,2)+h2]));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
scatter(norm_slope(:,1),norm_slope(:,2),50,'k','filled')
mn_slp1 = ((1./v1)./sum(1./v1)) .* norm_slope(:,1)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* norm_slope(:,2)';
scatter(sum(mn_slp1),sum(mn_slp2),200,'r','filled');
xlabel('Slopes for Short Frame ','FontSize',20);
ylabel('Slopes for Long Frame','FontSize',20);
hold('on');
hold('on');
axis('tight');
get(gca,'XTickLabel');
set(gca,'fontsize',18)
get(gca,'YTickLabel');
set(gca,'fontsize',18)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subaxis(2,4,4, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
norm_slope_time_locked(:,1) = norm_slope(:,1)/time(1);
norm_slope_time_locked(:,2) = norm_slope(:,2)/time(2);
norm_slope_all_time_locked(:,1,:) = norm_slope_all(:,1,:)/time(1);
norm_slope_all_time_locked(:,2,:) = norm_slope_all(:,2,:)/time(2);
scatter(norm_slope_time_locked(:,1),norm_slope_time_locked(:,2),50,'k','filled');
for sb=1:num_sub
    mn1(sb) = mean(squeeze(norm_slope_all_time_locked(sb,1,:)));
    mn2(sb) = mean(squeeze(norm_slope_all_time_locked(sb,2,:)));
    v1(sb) = (std(squeeze(norm_slope_all_time_locked(sb,1,:)))/sqrt(size(norm_slope_all_time_locked,3)))^2;%var(squeeze(norm_slope_all_time_locked(sb,1,:)));
    v2(sb) = (std(squeeze(norm_slope_all_time_locked(sb,2,:)))/sqrt(size(norm_slope_all_time_locked,3)))^2;%var(squeeze(norm_slope_all_time_locked(sb,2,:)));
    l1(sb) = prctile(squeeze(norm_slope_all_time_locked(sb,1,:)),50)-prctile(squeeze(norm_slope_all_time_locked(sb,1,:)),16);
    h1(sb) = prctile(squeeze(norm_slope_all_time_locked(sb,1,:)),84)-prctile(squeeze(norm_slope_all_time_locked(sb,1,:)),50);
    l2(sb) = prctile(squeeze(norm_slope_all_time_locked(sb,2,:)),50)-prctile(squeeze(norm_slope_all_time_locked(sb,2,:)),16);
    h2(sb) = prctile(squeeze(norm_slope_all_time_locked(sb,2,:)),84)-prctile(squeeze(norm_slope_all_time_locked(sb,2,:)),50);
end
hold on;
eb(1) = errorbar(norm_slope_time_locked(:,1),norm_slope_time_locked(:,2), l1, h1, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(norm_slope_time_locked(:,1),norm_slope_time_locked(:,2), l2, h2, 'vertical', 'LineStyle', 'none');
set(eb(1), 'color', 'm', 'LineWidth', 2)
set(eb(2), 'color', 'c', 'LineWidth', 2)
mn_s = min(min([norm_slope_time_locked(:,1)-l1 norm_slope_time_locked(:,2)-l2]));
mx_s = max(max([norm_slope_time_locked(:,1)+h1 norm_slope_time_locked(:,2)+h2]));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
scatter(norm_slope_time_locked(:,1),norm_slope_time_locked(:,2),50,'k','filled')
mn_slp1 = ((1./v1)./sum(1./v1)) .* norm_slope_time_locked(:,1)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* norm_slope_time_locked(:,2)';
scatter(sum(mn_slp1),sum(mn_slp2),200,'r','filled');
xlabel('Short time-locked slope ','FontSize',20);
ylabel('Long time-locked slope','FontSize',20);
hold('on');
axis('tight');
get(gca,'XTickLabel');
set(gca,'fontsize',18)
get(gca,'YTickLabel');
set(gca,'fontsize',18)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subaxis(2,4,5, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
time1 = 1/24;
time2 = 1/6;
for i=1:(num_sub)
    plot((0:num_frames-1)*time1,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'m');
    hold('on');
    plot((0:num_frames-1)*time2,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'c');
    hold('on');
end
hold('on');
xlabel('Time elapsed in msec');
ylabel('Weights');
hold('on');
plot((0:num_frames-1)*time1,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'-om','LineWidth',2);
hold on;
plot((0:num_frames-1)*time2,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'-oc','LineWidth',2);
hold on;
plot((0:num_frames-1)*time2,0.0*linspace(0,(num_frames-1)*time2,num_frames),'k','LineWidth',2);
axis('tight');
xticks([0 : 0.25:(num_frames-1)*time2])
a = get(gca,'XTickLabel');
set(gca,'fontsize',18)
b = get(gca,'YTickLabel');
set(gca,'fontsize',18);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
legend ({['Short time-locked'] ['Long time-locked']},'Box','off','Fontsize',10)

subaxis(2,4,6, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
for i=1:(num_sub)
    plot(vals,squeeze(subj_resp_pred(i,1,:)),'m','Linewidth',0.5);
    hold('on');
end
plot(vals,squeeze(mean(subj_resp_pred(:,1,:),1)),'-om','Linewidth',2);
hold('on');
yline(0.5,'k','Linewidth',1);
xline(0.0,'k','Linewidth',1);
axis('tight');
xlabel('Signed Kappa');
ylabel('Percent chose left for short stimulus');
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subaxis(2,4,7, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
for i=1:(num_sub)
    plot(vals,squeeze(subj_resp_pred(i,2,:)),'c','Linewidth',0.5);
    hold('on');
end
plot(vals,squeeze(mean(subj_resp_pred(:,2,:),1)),'-oc','Linewidth',2);
hold('on');
yline(0.5,'k','Linewidth',1);
xline(0.0,'k','Linewidth',1);
axis('tight');
xlabel('Signed Kappa');
ylabel('Percent chose left for long stimulus');
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

%% first half second half comparison for short frame trials
figure(num_sub+7)
time(1) = 1/24;
time(2) = 1/6;
axis image
subaxis(1,2,1, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel_first_half(i,1,1:num_frames)),'r');
    hold on;
    plot(1:num_frames,squeeze(norm_temporal_kernel_first_half(i,1,1:num_frames)),'b');
    hold on;
end
hold('on');
xlabel('Frames','FontSize',20);
ylabel('Weights','FontSize',20);
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel_first_half(:,1,1:num_frames))),1),'-or','LineWidth',2);
hold on;
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel_second_half(:,1,1:num_frames))),1),'-ob','LineWidth',2);
hold('on');
hold on;
plot(1:num_frames,zeros(1,num_frames),'k','LineWidth',2);
legend ({['Short frame first half'] ['Short frame second frame']},'Box','off','Fontsize',10);
axis('tight');
xticks([1 : 1:num_frames])
a = get(gca,'XTickLabel');
set(gca,'fontsize',18)
b = get(gca,'YTickLabel');
set(gca,'fontsize',18)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subaxis(1,2,2, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
scatter(norm_slope_first_half(:,1),norm_slope_second_half(:,1),50,'k','filled');
for sb=1:num_sub
    mn1(sb) = mean(squeeze(norm_slope_all_first_half(sb,1,:)));
    mn2(sb) = mean(squeeze(norm_slope_all_second_half(sb,1,:)));
    v1(sb) = (std(squeeze(norm_slope_all_first_half(sb,1,:)))/sqrt(size(norm_slope_all_first_half,3)))^2;%var(squeeze(norm_slope_all(sb,1,:)));
    v2(sb) = (std(squeeze(norm_slope_all_second_half(sb,1,:)))/sqrt(size(norm_slope_all_second_half,3)))^2;%var(squeeze(norm_slope_all(sb,2,:)));
    l1(sb) = prctile(squeeze(norm_slope_all_first_half(sb,1,:)),50)-prctile(squeeze(norm_slope_all_first_half(sb,1,:)),16);
    h1(sb) = prctile(squeeze(norm_slope_all_first_half(sb,1,:)),84)-prctile(squeeze(norm_slope_all_first_half(sb,1,:)),50);
    l2(sb) = prctile(squeeze(norm_slope_all_second_half(sb,1,:)),50)-prctile(squeeze(norm_slope_all_second_half(sb,1,:)),16);
    h2(sb) = prctile(squeeze(norm_slope_all_second_half(sb,1,:)),84)-prctile(squeeze(norm_slope_all_second_half(sb,1,:)),50);
end
hold on;
eb(1) = errorbar(norm_slope_first_half(:,1), norm_slope_second_half(:,1), l1, h1, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(norm_slope_first_half(:,1), norm_slope_second_half(:,1), l2, h2, 'vertical', 'LineStyle', 'none');
set(eb(1), 'color', 'r', 'LineWidth', 2)
set(eb(2), 'color', 'b', 'LineWidth', 2)
mn_s = min(min([norm_slope_first_half(:,1)-l1 norm_slope_second_half(:,1)-l2]));
mx_s = max(max([norm_slope_first_half(:,1)+h1 norm_slope_second_half(:,1)+h2]));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
scatter(norm_slope_first_half(:,1),norm_slope_second_half(:,1),50,'k','filled')
mn_slp1 = ((1./v1)./sum(1./v1)) .* norm_slope_first_half(:,1)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* norm_slope_second_half(:,1)';
scatter(sum(mn_slp1),sum(mn_slp2),200,'g','filled');
xlabel('Slopes for Short Frame First Half','FontSize',20);
ylabel('Slopes for Short Frame Second Half','FontSize',20);
hold('on');
hold('on');
axis('tight');
get(gca,'XTickLabel');
set(gca,'fontsize',18)
get(gca,'YTickLabel');
set(gca,'fontsize',18)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;



%% first and second half comparison for long frame trials
figure(num_sub+8)
time(1) = 1/24;
time(2) = 1/6;
axis image
subaxis(1,2,1, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel_first_half(i,2,1:num_frames)),'r');
    hold on;
    plot(1:num_frames,squeeze(norm_temporal_kernel_first_half(i,2,1:num_frames)),'b');
    hold on;
end
hold('on');
xlabel('Frames','FontSize',20);
ylabel('Weights','FontSize',20);
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel_first_half(:,2,1:num_frames))),1),'-or','LineWidth',2);
hold on;
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel_second_half(:,2,1:num_frames))),1),'-ob','LineWidth',2);
hold('on');
hold on;
plot(1:num_frames,zeros(1,num_frames),'k','LineWidth',2);
legend ({['Long frame first half'] ['Long frame second frame']},'Box','off','Fontsize',10);
axis('tight');
xticks([1 : 1:num_frames])
a = get(gca,'XTickLabel');
set(gca,'fontsize',18)
b = get(gca,'YTickLabel');
set(gca,'fontsize',18)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

subaxis(1,2,2, 'Spacing', 0.025, 'Padding', 0.025, 'Margin', 0.075);
scatter(norm_slope_first_half(:,2),norm_slope_second_half(:,2),50,'k','filled');
for sb=1:num_sub
    mn1(sb) = mean(squeeze(norm_slope_all_first_half(sb,2,:)));
    mn2(sb) = mean(squeeze(norm_slope_all_second_half(sb,2,:)));
    v1(sb) = (std(squeeze(norm_slope_all_first_half(sb,2,:)))/sqrt(size(norm_slope_all_first_half,3)))^2;%var(squeeze(norm_slope_all(sb,1,:)));
    v2(sb) = (std(squeeze(norm_slope_all_second_half(sb,2,:)))/sqrt(size(norm_slope_all_second_half,3)))^2;%var(squeeze(norm_slope_all(sb,2,:)));
    l1(sb) = prctile(squeeze(norm_slope_all_first_half(sb,2,:)),50)-prctile(squeeze(norm_slope_all_first_half(sb,2,:)),16);
    h1(sb) = prctile(squeeze(norm_slope_all_first_half(sb,2,:)),84)-prctile(squeeze(norm_slope_all_first_half(sb,2,:)),50);
    l2(sb) = prctile(squeeze(norm_slope_all_second_half(sb,2,:)),50)-prctile(squeeze(norm_slope_all_second_half(sb,2,:)),16);
    h2(sb) = prctile(squeeze(norm_slope_all_second_half(sb,2,:)),84)-prctile(squeeze(norm_slope_all_second_half(sb,2,:)),50);
end
hold on;
eb(1) = errorbar(norm_slope_first_half(:,2), norm_slope_second_half(:,2), l1, h1, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(norm_slope_first_half(:,2), norm_slope_second_half(:,2), l2, h2, 'vertical', 'LineStyle', 'none');
set(eb(1), 'color', 'r', 'LineWidth', 2)
set(eb(2), 'color', 'b', 'LineWidth', 2)
mn_s = min(min([norm_slope_first_half(:,2)-l1 norm_slope_second_half(:,2)-l2]));
mx_s = max(max([norm_slope_first_half(:,2)+h1 norm_slope_second_half(:,2)+h2]));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
scatter(norm_slope_first_half(:,2),norm_slope_second_half(:,2),50,'k','filled')
mn_slp1 = ((1./v1)./sum(1./v1)) .* norm_slope_first_half(:,2)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* norm_slope_second_half(:,2)';
scatter(sum(mn_slp1),sum(mn_slp2),200,'g','filled');
xlabel('Slopes for Long Frame First Half','FontSize',20);
ylabel('Slopes for Long Frame Second Half','FontSize',20);
hold('on');
hold('on');
axis('tight');
get(gca,'XTickLabel');
set(gca,'fontsize',18)
get(gca,'YTickLabel');
set(gca,'fontsize',18)
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
filename_save = strcat('SavedWorkspace/ShortLongFrames_',date,'.mat');
save(filename_save);
