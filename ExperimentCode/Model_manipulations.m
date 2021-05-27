close all; clear all; clc;
ps = 0.51:0.02:0.99;
pts = 35;
THRESHOLD = 0.7;
% updates = [5 10 20];
% fb_set= [0.2 1.0 2.0];
updates = [5 10 15];
fb_set= [0.1 0.5 1.0];
samples_set = [5 10 20];
models = {'is';'is_parallel';'is_race'};
models_txt1 = {'is';'is parallel';'is race';'samples and fb half'};
models_txt2 = {'is';'is parallel with cognitive load';'is race with cognitive load';'samples and fb half'};
num_frames = 10;
beta_range = [-.32 .1]; % min and max beta expected (to get maximum use of colorbar range)

%% Slopes with change in updates

for i=1:length(updates)
    disp(i);
    params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
    params.updates = floor(updates(i));
    sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, pts);
    [weights(i,:), errors(i,:), ~,~] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts(2,:));
    results = Model.runVectorized(params);
    [~, answers] = Model.genDataWithParams(params);
    perf(i) = mean(results.choices == answers);
    errors(i,:) = errors(i,:) ./ mean(weights(i,:));
    weights(i,:) = weights(i,:) ./ mean(weights(i,:));
end
close all;
figure();
for i=1:length(updates)
    txt = ['updates = ',num2str(updates(i))];
    errorbar(1:num_frames, weights(i,:),errors(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number');
    ylabel('Weight');
    hold on;
end

%% Fixed updates in clock time

for i=1:length(updates)
    params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
    params.constant_update = updates(ceil(length(updates)/2));
    params.updates = updates(i);
    sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, pts);
    [weights_clock(i,:), errors_clock(i,:), ~, ~] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts(2,:));
    results = Model.runVectorized(params);
    [~, answers] = Model.genDataWithParams(params);
    perf_clock(i) = mean(results.choices == answers);
    errors_clock(i,:) = errors_clock(i,:) ./ mean(weights_clock(i,1:end-1));
    weights_clock(i,:) = weights_clock(i,:) ./ mean(weights_clock(i,1:end-1));
end
close all;
figure();
for i=1:length(updates)
    txt = ['updates = ',num2str(updates(i))];
    errorbar(1:num_frames, weights_clock(i,:),errors_clock(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number');
    ylabel('Weight');
    hold on;
end

%% Change of slopes with change in fb strength

for i=1:length(fb_set)
    params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
    params.fb_strength = fb_set(i);
    sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, pts);
    [weights_fb(i,:), errors_fb(i,:), ~, ~] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts(2,:));
    results = Model.runVectorized(params);
    [~, answers] = Model.genDataWithParams(params);
    perf_fb(i) = mean(results.choices == answers);
    errors_fb(i,:) = errors_fb(i,:) ./ mean(weights_fb(i,1:end-1));
    weights_fb(i,:) = weights_fb(i,:) ./ mean(weights_fb(i,1:end-1));
end
close all;
figure();
for i=1:length(fb_set)
    txt = ['updates = ',num2str(fb_set(i))];
    errorbar(1:num_frames, weights_fb(i,:),errors_fb(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number');
    ylabel('Weight');
    hold on;
end

%% Change of slopes with change in number of samples

for i=1:length(samples_set)
    disp(i);
    params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', samples_set(i));
    sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, pts);
    [weights_smp(i,:), errors_smp(i,:), ~,~] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts(2,:));
    results = Model.runVectorized(params);
    [~, answers] = Model.genDataWithParams(params);
    perf_smp(i) = mean(results.choices == answers);
    errors_smp(i,:) = errors_smp(i,:) ./ mean(weights_smp(i,:));
    weights_smp(i,:) = weights_smp(i,:) ./ mean(weights_smp(i,:));
end
close all;
figure();
for i=1:length(samples_set)
    txt = ['samples = ',num2str(samples_set(i))];
    errorbar(1:num_frames, weights_smp(i,:),errors_smp(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number');
    ylabel('Weight');
    hold on;
end

%% Change of slopes for different models

for i=1:3
    disp(i);
    params = Model.newModelParams('model', models{i}, 'var_x', 0.1, 'gamma', 0.0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
    sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, pts);
    [weights_model(i,:), errors_model(i,:), ~,~] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts(2,:));
    results = Model.runVectorized(params);
    [~, answers] = Model.genDataWithParams(params);
    perf_model(i) = mean(results.choices == answers);
    errors_model(i,:) = errors_model(i,:) ./ mean(weights_model(i,:));
    weights_model(i,:) = weights_model(i,:) ./ mean(weights_model(i,:));
end
params = Model.newModelParams('model', models{1}, 'var_x', 0.1, 'gamma', 0.0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
params.samples = round(params.samples/2);
params.fb_strength = double(params.fb_strength/2);
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, pts);
[weights_model(4,:), errors_model(4,:), ~,~] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts(2,:));
results = Model.runVectorized(params);
[~, answers] = Model.genDataWithParams(params);
perf_model(4) = mean(results.choices == answers);
errors_model(4,:) = errors_model(4,:) ./ mean(weights_model(4,:));
weights_model(4,:) = weights_model(4,:) ./ mean(weights_model(4,:));
close all;
%%
figure();
for i=1:4
    txt = ['model = ',(models_txt1{i})];
    errorbar(1:num_frames, weights_model(i,:),errors_model(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number');
    ylabel('Weight');
    hold on;
end


for i=1:3
    disp(i);
    params = Model.newModelParams('model', models{i}, 'var_x', 0.1, 'gamma', 0.0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
    if i>=2
        params.cognitive_load = 1;
    end
    sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, pts);
    [weights_model_cog(i,:), errors_model_cog(i,:), ~,~] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts(2,:));
    results = Model.runVectorized(params);
    [~, answers] = Model.genDataWithParams(params);
    perf_model_cog(i) = mean(results.choices == answers);
    errors_model_cog(i,:) = errors_model_cog(i,:) ./ mean(weights_model_cog(i,:));
    weights_model_cog(i,:) = weights_model_cog(i,:) ./ mean(weights_model_cog(i,:));
end
params = Model.newModelParams('model', models{1}, 'var_x', 0.1, 'gamma', 0.0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
params.samples = round(params.samples/2);
params.fb_strength = double(params.fb_strength/2);
sens_cat_pts = Model.getThresholdPoints(ps, params, THRESHOLD, pts);
[weights_model_cog(4,:), errors_model_cog(4,:), ~,~] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts(2,:));
results = Model.runVectorized(params);
[~, answers] = Model.genDataWithParams(params);
perf_model_cog(4) = mean(results.choices == answers);
errors_model_cog(4,:) = errors_model_cog(4,:) ./ mean(weights_model_cog(4,:));
weights_model_cog(4,:) = weights_model_cog(4,:) ./ mean(weights_model_cog(4,:));
close all;

figure();
for i=1:4
    txt = ['model = ',(models_txt2{i})];
    errorbar(1:num_frames, weights_model_cog(i,:),errors_model_cog(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number');
    ylabel('Weight');
    hold on;
end
%% Combined Figure of Model Predictions for various manipulations

close all;
figure();
subplot(2,3,1);
for i=1:length(updates)
    txt = ['updates n_U = \tau = ',num2str(updates(i))];
    legend_txt{i} = [txt];
    errorbar(1:num_frames, weights(i,:),errors(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number','FontSize',20);
    ylabel('Weight','FontSize',20);
    hold on;
    axis('tight')
    
end
hold on;
ylim([-0.1 5])
xticks([1 : 1:num_frames])
plot(1:num_frames,zeros(1,num_frames),'k','LineWidth',2);
legend(legend_txt,'FontSize',15);
a = get(gca,'XTickLabel');
set(gca,'fontsize',18)
b = get(gca,'YTickLabel');
set(gca,'fontsize',18)
title('Unchanged slope for change of updates','fontsize',15)

subplot(2,3,2);
for i=1:length(updates)
    txt = ['updates n_U = ',num2str(updates(i)) ' and \tau = ', num2str(updates(ceil(length(updates)/2)))];
    legend_txt{i} = [txt];
    errorbar(1:num_frames, weights_clock(i,:),errors_clock(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number','FontSize',20);
    ylabel('Weight','FontSize',20);
    hold on;
    axis('tight')
end
hold on;
ylim([-0.1 5])
xticks([1 : 1:num_frames])
plot(1:num_frames,zeros(1,num_frames),'k','LineWidth',2);
legend(legend_txt,'FontSize',15);
a = get(gca,'XTickLabel');
set(gca,'fontsize',18)
b = get(gca,'YTickLabel');
set(gca,'fontsize',18)
title(['Fixed amortization rate of ' num2str(updates(ceil(length(updates)/2))) ' causes slope change'],'fontsize',15)

subplot(2,3,3);
for i=1:length(fb_set)
    txt = ['feedback strength = ' num2str(fb_set(i))];
    legend_txt{i} = [txt];
    errorbar(1:num_frames, weights_fb(i,:),errors_fb(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number','FontSize',20);
    ylabel('Weight','FontSize',20);
    hold on;
    axis('tight')
end
hold on;
ylim([-0.1 5])
xticks([1 : 1:num_frames])
plot(1:num_frames,zeros(1,num_frames),'k','LineWidth',2);
legend(legend_txt,'FontSize',15);
a = get(gca,'XTickLabel');
set(gca,'fontsize',18)
b = get(gca,'YTickLabel');
set(gca,'fontsize',18)
title('Slope steeper for stronger feedback','fontsize',15)

subplot(2,3,4);
for i=1:length(samples_set)
    txt = ['samples = ' num2str(samples_set(i))];
    legend_txt{i} = [txt];
    errorbar(1:num_frames, weights_smp(i,:),errors_smp(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number','FontSize',20);
    ylabel('Weight','FontSize',20);
    hold on;
    axis('tight')
end
hold on;
ylim([-0.1 5])
xticks([1 : 1:num_frames])
plot(1:num_frames,zeros(1,num_frames),'k','LineWidth',2);
legend(legend_txt,'FontSize',15);
a = get(gca,'XTickLabel');
set(gca,'fontsize',18)
b = get(gca,'YTickLabel');
set(gca,'fontsize',18)
title('Slope steeper for lesser samples','fontsize',15)

subplot(2,3,5);
% figure();
for i=1:4
    txt = ['model = ' (models_txt1{i})];
    legend_txt{i} = [txt];
    errorbar(1:num_frames, weights_model(i,:),errors_model(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number','FontSize',20);
    ylabel('Weight','FontSize',20);
    hold on;
    axis('tight')
end
hold on;
ylim([-0.1 10])
xticks([1 : 1:num_frames])
plot(1:num_frames,zeros(1,num_frames),'k','LineWidth',2);
legend(legend_txt,'FontSize',15);
a = get(gca,'XTickLabel');
set(gca,'fontsize',18)
b = get(gca,'YTickLabel');
set(gca,'fontsize',18)
title('Slope changes for different models','fontsize',15)

subplot(2,3,6);
% figure();
for i=1:4
    txt = ['model = ' (models_txt2{i})];
    legend_txt{i} = [txt];
    errorbar(1:num_frames, weights_model_cog(i,:),errors_model_cog(i,:),'DisplayName',txt,'LineWidth',2);
    xlabel('Frame Number','FontSize',20);
    ylabel('Weight','FontSize',20);
    hold on;
    axis('tight')
end
hold on;
ylim([-0.1 6.5])
xticks([1 : 1:num_frames])
plot(1:num_frames,zeros(1,num_frames),'k','LineWidth',2);
legend(legend_txt,'FontSize',15);
a = get(gca,'XTickLabel');
set(gca,'fontsize',18)
b = get(gca,'YTickLabel');
set(gca,'fontsize',18)
title('Slope changes for different models','fontsize',15)