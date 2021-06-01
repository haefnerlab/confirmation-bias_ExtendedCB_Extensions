clear all; close all; clc;
boots = 500;% number of bootstraps to get PK

hpr_ridge = logspace(-1, 5, 7);
hpr_ar1 = 0.0;
hpr_curvature = logspace(-1, 5, 7);
expt_type_noise = 2;
expt_type_ratio = 1;
standardize = 0; % z-score data or not
folds = 10; % how many folds of cross validation
num_frames = 10; % number of stimulus frames
cases = 2; % ratio case and noise case, so 2 cases. phase 2 for noise and 1 ratio

dir = 'RawDataLongDuration_NoiseRatio';
subjects = {...
    'bpglongerframes-subject12';
    'bpglongerframes-subject18'
    };

disp('Preloading big stimulus statistics for faster post computations....');
[num_sub,~] = size(subjects);

disp('Starting to find best hyperparameters for NOISE data across subjects....');
[best_hprs_noise] = CustomRegression.combined_hprs_search(subjects, expt_type_noise, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds, dir);
disp('Starting to find best hyperparameters for RATIO data across subjects....');
[best_hprs_ratio] = CustomRegression.combined_hprs_search(subjects, expt_type_ratio, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds, dir);

for i=1:(num_sub)
    figure()
    for j=1:cases
        phase = 3-j;
        if phase==expt_type_noise
            best_hprs = best_hprs_noise;
            disp(['Starting analysis for NOISE data of Subject ' num2str(i) ' ...']);
        elseif phase==expt_type_ratio
            best_hprs = best_hprs_ratio;
            disp(['Starting analysis for RATIO data of Subject ' num2str(i) ' ...']);
        end
        [params_boot,sobl,abbl_exp,trials,bin_centers,means,stderrs,data,log_bernoulli{i,j}] = run_analysis_both(subjects{i}, phase, boots, best_hprs, dir);
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
        slope(i,j) = prctile(squeeze(all_linear(i,j,:,1)),50);
        slope_all(i,j,:) = sobl(:,1);
        norm_slope_all(i,j,:) = norm_all_linear(i,j,:,1);
        norm_slope(i,j) = prctile(squeeze(norm_all_linear(i,j,:,1)),50);
        hprs_used(i,j,:) = best_hprs;
        data_sub{i,j} = data;
        
        subplot(cases,3,3+3*(j-1))
        errorbar(1:num_frames,squeeze(temporal_kernel(i,j,1:num_frames)),squeeze(lo_temporal_kernel(i,j,:)),squeeze(hi_temporal_kernel(i,j,:)),'k','LineWidth',1)
        xlabel('Frames');
        ylabel('Weights');
        axis('tight');
        
        if phase==2
            
            subplot(cases,3,1+3*(j-1))
            plot((1:length(data.noise)), data.noise);
            hold on;
            plot(linspace(1,length(data.noise),100),ones(1,100),'r','LineWidth',1);
            xlabel('Trials');
            ylabel('Noise Level');
            axis('tight');
            
            subplot(cases,3,2+3*(j-1))
            bins = 10;
            subject_pm_curve =[];
            uniq_vals = linspace(-0.8,0.8,10);
            tr_kappa = data.sign_noise;
            noise_signal = data.ideal_frame_signals;
            for tt=1:(length(uniq_vals)-1)
                subj_resp(i,j,tt) = mean(data.choice(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)));
                ntrial_subj(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
            end
            noise_vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
            errorbar(noise_vals,squeeze(subj_resp(i,j,:)),squeeze((subj_resp(i,j,:)).*(1-subj_resp(i,j,:))./sqrt(ntrial_subj(i,j,:))),'o');
            subject_pm_curve = (1./(1+exp(-(noise_signal*squeeze(temporal_kernel(i,j,:))+bias(i,j)))))*( 1-(alpha(i,j,1)))+(alpha(i,j,1)/2);
            for tt=1:(length(uniq_vals)-1)
                subj_resp_pred(i,j,tt) = mean(subject_pm_curve(((tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1)))));
                ntrial_subj_pred(i,j,tt) = sum(tr_kappa>uniq_vals(tt)&tr_kappa<=uniq_vals(tt+1));
            end
            hold on;
            plot(noise_vals,squeeze(subj_resp_pred(i,j,:)),'Linewidth',2);
            yline(0.5,'--k');
            xline(0.0,'--k');
            xlabel('Signed Kappa');
            ylabel('Percent chose left');
            xlim([-0.8 0.8])
            ylim([0.0 1.0])
            
        else
            
            subplot(cases,3,1+3*(j-1))
            plot((1:length(data.ratio)), data.ratio);
            hold on;
            plot(linspace(1,length(data.ratio),100),ones(1,100),'r','LineWidth',1);
            xlabel('Trials');
            ylabel('Ratio Level');
            axis('tight');
            
            subplot(cases,3,2+3*(j-1))
            bins = 10;
            subject_pm_curve =[];
            uniq_vals = linspace(0,1,10);
            tr_ratio = data.true_ratio;
            noise_signal = data.ideal_frame_signals;
            for tt=1:(length(uniq_vals)-1)
                subj_resp(i,j,tt) = mean(data.choice(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)));
                ntrial_subj(i,j,tt) = sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
            end
            ratio_vals = uniq_vals(1:end-1) + (uniq_vals(2) - uniq_vals(1))/2;
            errorbar(ratio_vals,squeeze(subj_resp(i,j,:)),squeeze((subj_resp(i,j,:)).*(1-subj_resp(i,j,:))./sqrt(ntrial_subj(i,j,:))),'o');
            subject_pm_curve = (1./(1+exp(-(noise_signal*squeeze(temporal_kernel(i,j,:))+bias(i,j)))))*( 1-(alpha(i,j,1)))+(alpha(i,j,1)/2);
            for tt=1:(length(uniq_vals)-1)
                subj_resp_pred(i,j,tt) = mean(subject_pm_curve(((tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1)))));
                ntrial_subj_pred(i,j,tt) = sum(tr_ratio>uniq_vals(tt)&tr_ratio<=uniq_vals(tt+1));
            end
            hold on;
            plot(ratio_vals,squeeze(subj_resp_pred(i,j,:)),'Linewidth',2);
            yline(0.5,'--k');
            xline(0.5,'--k');
            xlabel('Ratio');
            ylabel('Percent chose left');
            xlim([0.0 1.0])
            ylim([0.0 1.0])
            title(['Fitted Psychometric Curve'])
            
        end
    end
    sgtitle(['Top Row: Noise and Bottom Row: Longer Frames for Subject ' num2str(i)])
end

%%

figure()
for i=1:(num_sub)
    subplot(2,5,i);
    errorbar(1:num_frames,squeeze(temporal_kernel(i,1,1:num_frames)),squeeze(lo_temporal_kernel(i,1,:)),squeeze(hi_temporal_kernel(i,1,:)),'m','LineWidth',2);
    xlabel('Frames');
    ylabel('Weights');
    hold('on');
    errorbar(1:num_frames,squeeze(temporal_kernel(i,2,1:num_frames)),squeeze(lo_temporal_kernel(i,2,:)),squeeze(hi_temporal_kernel(i,2,:)),'c','LineWidth',2);
    legend({['short' '(' num2str(num_trials(i,1)) ' trials)'],['long' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Noise Vs Ratio for Ratio for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end
figure()
for i=1:(num_sub)
    subplot(2,5,i);
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'m','LineWidth',2);
    xlabel('Frames');
    ylabel('Norm Weights');
    hold('on');
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'c','LineWidth',2);
    legend({['short' '(' num2str(num_trials(i,1)) ' trials)'],['long' '(' num2str(num_trials(i,2)) ' trials)']});
    title(['Noise Vs Ratio for Ratio for Subject ' num2str(i)])
    hold('on');
    axis('tight');
end


%%
figure()
ax1=subplot(2,5,1);
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
title('Noise Weights')

ax2=subplot(2,5,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(temporal_kernel(i,2,1:num_frames)),'c');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Weights');
hold('on');
plot(1:num_frames,mean(squeeze((temporal_kernel(:,2,1:num_frames))),1),'c','LineWidth',2);
title('Ratio Weights')
hold('on');
axis('tight')

ax3=subplot(2,5,3);
for i=1:(num_sub)
    plot((1:2),slope(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, slope(i,1), 'm');
    hold on;
    scatter(2, slope(i,2), 'c');
    hold on;
end
hold('on');
xticklabels({'Noise','','Ratio'});
xlabel('Noise (blue) to Ratio (green)');
ylabel('Slopes');
hold('on');
title('Slope Comparison')
hold('on');

ax4=subplot(2,5,4);
scatter(slope(:,1),slope(:,2),'r');
mn_s = min(min(slope));
mx_s = max(max(slope));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
scatter(mean(mean(slope_all(:,1,:))),mean(mean(slope_all(:,2,:))),100,'r','filled');
xlabel('Slopes for Noise ');
ylabel('Slopes for Ratio');
hold('on');
title('Slope Comparison')
hold('on');
axis('tight')

ax5=subplot(2,5,5);
scatter(beta(:,1),beta(:,2),'r');
mn_b = min(min(beta));
mx_b = max(max(beta));
hold on;
plot(linspace(mn_b,mx_b,10),linspace(mn_b,mx_b,10),'k','LineWidth',2);
hold on;
scatter(mean(mean(beta_all(:,1,:))),mean(mean(beta_all(:,2,:))),100,'r','filled');
xlabel('Beta for Noise ');
ylabel('Beta for Ratio');
hold('on');
title('Beta Comparison')
hold('on');
axis('tight')
linkaxes([ax1,ax2],'y')

ax6=subplot(2,5,6);
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
title('Noise Norm Weights')

ax7=subplot(2,5,7);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'c');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'c','LineWidth',2);
title('Ratio Norm Weights')
hold('on');
axis('tight')

ax8=subplot(2,5,8);
for i=1:(num_sub)
    plot((1:2),norm_slope(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, norm_slope(i,1), 'm');
    hold on;
    scatter(2, norm_slope(i,2), 'c');
    hold on;
end
hold('on');
xticklabels({'Noise','','Ratio'});
xlabel('Noise (blue) to Ratio (green)');
ylabel('Norm Slopes');
hold('on');
title('Norm Slope Comparison')
hold('on');

ax9=subplot(2,5,9);
scatter(norm_slope(:,1),norm_slope(:,2),'r');
mn_s = min(min(norm_slope));
mx_s = max(max(norm_slope));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
scatter(mean(norm_slope(:,1)),mean(norm_slope(:,2)),100,'r','filled');
xlabel('Norm Slopes for Noise ');
ylabel('Norm Slopes for Ratio');
hold('on');
title('Norm Slope Comparison')
hold('on');
axis('tight')

ax10=subplot(2,5,10);
scatter(beta(:,1),beta(:,2),'r');
mn_b = min(min(beta));
mx_b = max(max(beta));
hold on;
plot(linspace(mn_b,mx_b,10),linspace(mn_b,mx_b,10),'k','LineWidth',2);
hold on;
scatter(mean(beta(:,1)),mean(beta(:,2)),100,'r','filled');
xlabel('Beta for Noise ');
ylabel('Beta for Ratio');
hold('on');
title('Beta Comparison')
hold('on');
axis('tight')
linkaxes([ax6,ax7],'y')

%% Summary Fig
ax1=subplot(2,4,1);
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
title('Noise Norm Weights')

ax2=subplot(2,4,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'c');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'c','LineWidth',2);
title('Ratio Norm Weights')
hold('on');
axis('tight')

ax3=subplot(2,4,3);
for i=1:(num_sub)
    plot((1:2),norm_slope(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, norm_slope(i,1), 'm');
    hold on;
    scatter(2, norm_slope(i,2), 'c');
    hold on;
end
hold('on');
xticklabels({'Noise','','Ratio'});
ylabel('Slopes');
hold('on');
title('Norm Slope Comparison')
hold('on');

ax4=subplot(2,4,4);
for i=1:(num_sub)
    plot((1:2),beta(i,:),'-ok','LineWidth',2);
    hold on;
    scatter(1, beta(i,1), 'm');
    hold on;
    scatter(2, beta(i,2), 'c');
    hold on;
end
hold('on');
xticklabels({'Noise','','Ratio'});
ylabel('Beta');
hold('on');
title('Beta Comparison')
hold('on');

ax5=subplot(2,4,5);
scatter(norm_slope(:,1),norm_slope(:,2),'k');
mn_s = min(min(norm_slope));
mx_s = max(max(norm_slope));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
scatter(mean(norm_slope(:,1)),mean(norm_slope(:,2)),100,'k','filled');
xlabel('Norm Slopes for Noise ');
ylabel('Norm Slopes for Ratio');
hold('on');
title('Norm Slope Comparison')
hold('on');
axis('tight')

ax6=subplot(2,4,6);
scatter(beta(:,1),beta(:,2),'k');
mn_b = min(min(beta));
mx_b = max(max(beta));
hold on;
plot(linspace(mn_b,mx_b,10),linspace(mn_b,mx_b,10),'k','LineWidth',2);
hold on;
scatter(mean(beta(:,1)),mean(beta(:,2)),100,'k','filled');
xlabel('Beta for Noise ');
ylabel('Beta for Ratio');
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

xlabel('Norm Slope for Noise ');
ylabel('Norm Slope for Ratio');
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

xlabel('Beta for Noise ');
ylabel('Beta for Ratio');
axis tight;
title('Beta');
%% Precise Summary Figure

figure(1000)
ax1=subplot(2,3,1);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,1,1:num_frames)),'m');
    hold('on');
end
hold('on');
axis('tight')
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,1,1:num_frames))),1),'-om','LineWidth',2);
title('Noise Norm Weights')

ax2=subplot(2,3,2);
for i=1:(num_sub)
    plot(1:num_frames,squeeze(norm_temporal_kernel(i,2,1:num_frames)),'c');
    hold('on');
end
hold('on');
xlabel('Frames');
ylabel('Norm Weights');
hold('on');
plot(1:num_frames,mean(squeeze((norm_temporal_kernel(:,2,1:num_frames))),1),'-oc','LineWidth',2);
title('Ratio Norm Weights')
hold('on');
axis('tight')

linkaxes([ax1,ax2, ax3],'y')

ax5=subplot(2,3,3);
scatter(norm_slope(:,1),norm_slope(:,2),'k','filled');
err1 = std(squeeze(norm_slope_all(:,1,:))',1)./sqrt(boots);
err2 = std(squeeze(norm_slope_all(:,2,:))',1)./sqrt(boots);
v1 = var(squeeze(norm_slope(:,1))',1);
v2 = var(squeeze(norm_slope(:,2))',1);
hold on;
eb(1) = errorbar(norm_slope(:,1),norm_slope(:,2),err1, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(norm_slope(:,1),norm_slope(:,2),err2, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'k', 'LineWidth', 2)
mn_s = min(min(norm_slope));
mx_s = max(max(norm_slope));
hold on;
plot(linspace(mn_s,mx_s,10),linspace(mn_s,mx_s,10),'k','LineWidth',2);
hold on;
mn_slp1 = ((1./v1)./sum(1./v1)) .* norm_slope(:,1)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* norm_slope(:,2)';
scatter(sum(mn_slp1),sum(mn_slp2),200,'r','filled');
xlabel('Norm Slopes for Noise ');
ylabel('Norm Slopes for Ratio');
hold('on');
title('Norm Slope Comparison')
axis('tight')

ax6=subplot(2,3,4);
scatter(beta(:,1),beta(:,2),'k','filled');
err1 = std(squeeze(beta_all(:,1,:))',1)./sqrt(boots);
err2 = std(squeeze(beta_all(:,2,:))',1)./sqrt(boots);
v1 = var(squeeze(beta(:,1))',1);
v2 = var(squeeze(beta(:,2))',1);
hold on;
eb(1) = errorbar(beta(:,1),beta(:,2),err1, 'horizontal', 'LineStyle', 'none');
hold on;
eb(2) = errorbar(beta(:,1),beta(:,2),err2, 'vertical', 'LineStyle', 'none');
set(eb, 'color', 'k', 'LineWidth', 2)
mn_b = min(min(beta));
mx_b = max(max(beta));
hold on;
plot(linspace(mn_b,mx_b,10),linspace(mn_b,mx_b,10),'k','LineWidth',2);
hold on;
ss = size(beta(:,1),1);
mn_slp1 = ((1./v1)./sum(1./v1)) .* norm_slope(:,1)';
mn_slp2 = ((1./v2)./sum(1./v2)) .* norm_slope(:,2)';
scatter(sum(mn_slp1),sum(mn_slp2),200,'r','filled');
xlabel('Beta for Noise ');
ylabel('Beta for Ratio');
hold('on');
title('Beta Comparison')
hold('on');
axis('tight')

subplot(2,3,5)
for i=1:(num_sub)
    plot(noise_vals,squeeze(subj_resp_pred(i,1,:)),'m','Linewidth',0.5);
    hold('on');
end
plot(noise_vals,squeeze(mean(subj_resp_pred(:,1,:),1)),'-om','Linewidth',2);
hold('on');
yline(0.5,'k','Linewidth',1);
xline(0.0,'k','Linewidth',1);
axis('tight')
xlabel('Signed Kappa');
ylabel('Percent chose left for noise stimulus');

subplot(2,3,6)
for i=1:(num_sub)
    plot(ratio_vals,squeeze(subj_resp_pred(i,2,:)),'c','Linewidth',0.5);
    hold('on');
end
plot(ratio_vals,squeeze(mean(subj_resp_pred(:,2,:),1)),'-oc','Linewidth',2);
hold('on');
yline(0.5,'k','Linewidth',1);
xline(0.5,'k','Linewidth',1);
axis('tight')
xlabel('Ratio');
ylabel('Percent chose left for ratio stimulus');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename_save = strcat('SavedWorkspace/ShortLongFrames_NoiseRatio_',date,'.mat');
save(filename_save);
