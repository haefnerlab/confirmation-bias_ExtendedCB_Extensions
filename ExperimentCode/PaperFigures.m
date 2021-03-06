%% Setup

clear;
RATIO_PHASE = 1; % aka HSLC
NOISE_PHASE = 2; % aka LSHC
THRESHOLD = 0.7;
DATADIR = fullfile(pwd, '..', 'RawData');
MEMODIR = fullfile(pwd, '..', 'Precomputed');

DARK_RED = [204 0 0] / 255;
DARK_BLUE = [32 74 135] / 255;
% fade 50% towards white
LIGHT_RED = DARK_RED * .5 + .5;
LIGHT_BLUE = DARK_BLUE * .5 + .5;

ratioSubjects = arrayfun(@(i) sprintf('bpgFinaltest-subject%02d', i), setdiff(1:15, [1 5 12]), 'UniformOutput', false);
noiseSubjects = arrayfun(@(i) sprintf('bpgFinaltest-subject%02d', i), setdiff(1:15, [1 5]), 'UniformOutput', false);

% informed is the opposite of naive
informedSubjects = arrayfun(@(i) sprintf('bpgFinaltest-subject%02d', i), [7 8 9], 'UniformOutput', false);

naiveRatioSubjects = setdiff(ratioSubjects, informedSubjects);
naiveNoiseSubjects = setdiff(noiseSubjects, informedSubjects);
naiveBothSubjects = intersect(naiveRatioSubjects, naiveNoiseSubjects);

% Re-order so that informed subjects are at the end
ratioSubjects = [naiveRatioSubjects informedSubjects];
noiseSubjects = [naiveNoiseSubjects informedSubjects];
bothSubjects = [intersect(naiveRatioSubjects, naiveNoiseSubjects) informedSubjects];

% Compute population-level psycho metrics
N = length(bothSubjects);

for n=N:-1:1
    hyphenSplit = strsplit(bothSubjects{n}, '-');
    shortnames{n} = hyphenSplit{2};
end

%% Figure 1

% Nothing to do - conceptual figures only created in vector graphics editor

%% Figure 2 - data analysis

sem = @(vals) std(vals)/sqrt(N);
for iSubject=length(bothSubjects):-1:1
    warning off;
    [~,~,psycho_fit] = GaborAnalysis.getThresholdWindow(bothSubjects{iSubject}, NOISE_PHASE, 0.5, THRESHOLD, DATADIR);
    warning on;
    % See @plotPsych, e.g.
    pc_zero_signal(iSubject) = (1-psycho_fit.Fit(3)-psycho_fit.Fit(4))*psycho_fit.options.sigmoidHandle(0,psycho_fit.Fit(1),psycho_fit.Fit(2))+psycho_fit.Fit(4);
    threshold(iSubject) = getThreshold(psycho_fit, 0.75);
    
    [~,~,psycho_fit] = GaborAnalysis.getThresholdWindow(bothSubjects{iSubject}, RATIO_PHASE, 0.3, 0.7, DATADIR);
    slope = getSlopePC(psycho_fit, 0.5);
    pc_60_40_adjusted(iSubject) = 0.5 + .1 * slope;
end
fprintf('Noise condition mean threshold = %f +/- %f, stdev = %f\n', mean(threshold), sem(threshold), std(threshold));
fprintf('Noise condition mean %%correct at zero signal = %f +/- %f, stdev = %f\n', mean(100*pc_zero_signal), sem(100*pc_zero_signal), std(100*pc_zero_signal));
fprintf('Ratio condition mean %%correct at 6:4 = %f +/- %f, stdev = %f\n', mean(100*pc_60_40_adjusted), sem(100*pc_60_40_adjusted), std(100*pc_60_40_adjusted));

GaborAnalysis.DeltaSlopeStatistics(bothSubjects, [1 2], 'exponential');
GaborAnalysis.DeltaSlopeStatistics(bothSubjects, [1 2], 'linear');
GaborAnalysis.DeltaPK(bothSubjects, [1 2], false);

%% Psychometric curves

figure;

% Noise condition
subplot(1,2,1);
hold on;
noises = linspace(-0.8, 0.8, 101);
avg_pm_curve = zeros(size(noises));
for iSubject=1:N
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, NOISE_PHASE);
    [pm_fit, ~, ~, ~] = LoadOrRun(@GaborPsychometric, {SubjectData, -2}, ...
        fullfile(MEMODIR, ['PM-noise-signed-' bothSubjects{iSubject} '.mat']));
    % Next line copied from @plotPsych in the psignifit toolbox
    subject_pm_curve = (1-pm_fit.Fit(3)-pm_fit.Fit(4))*arrayfun(@(x) pm_fit.options.sigmoidHandle(x,pm_fit.Fit(1),pm_fit.Fit(2)), noises)+pm_fit.Fit(4);
    avg_pm_curve = avg_pm_curve + (subject_pm_curve - avg_pm_curve) / iSubject;
    plot(noises, subject_pm_curve, 'Color', LIGHT_RED, 'LineWidth', 1);
end
plot(noises, avg_pm_curve, 'Color', DARK_RED, 'LineWidth', 2);
xlabel('noise level (\kappa)');
ylabel('percent chose left');

% Ratio condition
subplot(1,2,2);
hold on;
ratios = linspace(0, 1, 101);
avg_pm_curve = zeros(size(ratios));
for iSubject=1:N
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, RATIO_PHASE);
    [pm_fit, ~, ~, ~] = LoadOrRun(@GaborPsychometric, {SubjectData, RATIO_PHASE}, ...
        fullfile(MEMODIR, ['PM-true_ratio-signed-' bothSubjects{iSubject} '.mat']));
    % Next line copied from @plotPsych in the psignifit toolbox
    subject_pm_curve = (1-pm_fit.Fit(3)-pm_fit.Fit(4))*arrayfun(@(x) pm_fit.options.sigmoidHandle(x,pm_fit.Fit(1),pm_fit.Fit(2)), ratios)+pm_fit.Fit(4);
    avg_pm_curve = avg_pm_curve + (subject_pm_curve - avg_pm_curve) / iSubject;
    plot(ratios, subject_pm_curve, 'Color', LIGHT_BLUE, 'LineWidth', 1);
end
plot(ratios, avg_pm_curve, 'Color', DARK_BLUE, 'LineWidth', 2);
xlabel('frame ratio left:right');
ylabel('percent chose left');

%% Supplemental PK cross-validation figure

nFold = 20;
pkLLNoise = zeros(4, N, nFold);
pkLLRatio = zeros(4, N, nFold);
for iSubject=N:-1:1
    % Noise experiment PK cross-validation
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, NOISE_PHASE);
    [floor, thresh] = GaborAnalysis.getThresholdWindow(bothSubjects{iSubject}, NOISE_PHASE, 0.5, THRESHOLD);
    SubjectDataThresh = GaborThresholdTrials(SubjectData, NOISE_PHASE, thresh, floor);
    sigs = SubjectDataThresh.ideal_frame_signals;
    resps = SubjectDataThresh.choice == +1;
    memo_name = ['PK-xValidCompare-noise-' bothSubjects{iSubject} '-' num2str(thresh) '-' num2str(floor) '.mat'];
    [pkLLNoise(1, iSubject, :), pkLLNoise(2, iSubject, :), pkLLNoise(3, iSubject, :), pkLLNoise(4, iSubject, :), ~] = ...
        LoadOrRun(@CustomRegression.xValidatePKModels, {sigs, resps, nFold}, fullfile(MEMODIR, memo_name));
    
    % Ratio experiment PK cross-validation
    SubjectData = LoadAllSubjectData(bothSubjects{iSubject}, RATIO_PHASE);
    SubjectDataThresh = GaborThresholdTrials(SubjectData, RATIO_PHASE, .6, .4);
    sigs = SubjectDataThresh.ideal_frame_signals;
    resps = SubjectDataThresh.choice == +1;
    disp(bothSubjects{iSubject});
    disp(size(sigs));
    disp(size(resps));
    memo_name = ['PK-xValidCompare-true_ratio-' bothSubjects{iSubject} '-0.6-0.4.mat'];
    [pkLLRatio(1, iSubject, :), pkLLRatio(2, iSubject, :), pkLLRatio(3, iSubject, :), pkLLRatio(4, iSubject, :), ~] = ...
        LoadOrRun(@CustomRegression.xValidatePKModels, {sigs, resps, nFold}, fullfile(MEMODIR, memo_name));
end

% Note: CustomRegression.xValidatePKModels gives paired results. That is, each of the 4 models is
% tested on the exact same data fold. The relevant statistic is therefore the *average of
% differences in LL* rather than the *difference in averages*.

% Compute differences w.r.t. linear model as the baseline
baselineModel = 4;
pkLLNoiseDiffs = pkLLNoise - pkLLNoise(baselineModel, :, :);
pkLLRatioDiffs = pkLLRatio - pkLLRatio(baselineModel, :, :);

% Sum across subjects for combined result
pkLLNoiseDiffsAll = squeeze(sum(pkLLNoiseDiffs, 2));
pkLLRatioDiffsAll = squeeze(sum(pkLLNoiseDiffs, 2));

% Sort across 'folds' so it's easy to get confidence intervals
pkLLNoiseDiffs = sort(pkLLNoiseDiffs, 3);
pkLLRatioDiffs = sort(pkLLRatioDiffs, 3);
pkLLNoiseDiffsAll = sort(pkLLNoiseDiffsAll, 2);
pkLLRatioDiffsAll = sort(pkLLRatioDiffsAll, 2);

% Get mean and 50% confidence intervals across folds for each statistic
meanDiffLLNoise = mean(pkLLNoiseDiffs, 3);
loLLNoise = pkLLNoiseDiffs(:, :, round(.25 * nFold));
hiLLNoise = pkLLNoiseDiffs(:, :, round(.75 * nFold));

meanDiffLLRatio = mean(pkLLRatioDiffs, 3);
loLLRatio = pkLLRatioDiffs(:, :, round(.25 * nFold));
hiLLRatio = pkLLRatioDiffs(:, :, round(.75 * nFold));

meanDiffLLNoiseAll = mean(pkLLNoiseDiffsAll, 2);
loLLNoiseAll = pkLLNoiseDiffsAll(:, round(.25 * nFold));
hiLLNoiseAll = pkLLNoiseDiffsAll(:, round(.75 * nFold));

meanDiffLLRatioAll = mean(pkLLRatioDiffsAll, 2);
loLLRatioAll = pkLLRatioDiffsAll(:, round(.25 * nFold));
hiLLRatioAll = pkLLRatioDiffsAll(:, round(.75 * nFold));

figure;
subplot(2,1,1);
hold on;
p1 = bar([meanDiffLLNoise'; meanDiffLLNoiseAll']);
drawnow;
for iSubject=1:N
    for iType=1:4
        c = p1(iType).FaceColor;
        x = p1(iType).XData(iSubject) + p1(iType).XOffset;
        m = meanDiffLLNoise(iType, iSubject);
        errorbar(x, m, m-loLLNoise(iType, iSubject), hiLLNoise(iType, iSubject)-m, 'Color', 'k');
    end
end
for iType=1:4
    c = p1(iType).FaceColor;
    x = N + 1 + p1(iType).XOffset;
    m = meanDiffLLNoiseAll(iType);
    errorbar(x, m, m-loLLNoiseAll(iType), hiLLNoiseAll(iType)-m, 'Color', 'k');
end
legend({'Unregularized LR', 'Regularized LR', 'Exponential', 'Linear'}, 'location', 'best');
ylabel('Relative log-likelihood');
title('Noise Condition');
set(gca, 'XTick', 1:N+1, 'XTickLabel', [shortnames, {'Combined'}]);
xtickangle(45);

subplot(2,1,2);
hold on;
p1 = bar([meanDiffLLRatio'; meanDiffLLRatioAll']);
drawnow;
for iSubject=1:N
    for iType=1:4
        c = p1(iType).FaceColor;
        x = p1(iType).XData(iSubject) + p1(iType).XOffset;
        m = meanDiffLLRatio(iType, iSubject);
        errorbar(x, m, m-loLLRatio(iType, iSubject), hiLLRatio(iType, iSubject)-m, 'Color', 'k');
    end
end
for iType=1:4
    c = p1(iType).FaceColor;
    x = N + 1 + p1(iType).XOffset;
    m = meanDiffLLRatioAll(iType);
    errorbar(x, m, m-loLLRatioAll(iType), hiLLRatioAll(iType)-m, 'Color', 'k');
end
ylabel('Relative log-likelihood');
title('Ratio Condition');
set(gca, 'XTick', 1:N+1, 'XTickLabel', [shortnames, {'Combined'}]);
xtickangle(45);

%% Figure 3 - model results

fig3 = figure;

% Range of values for category and sensory information
ps = 0.51:0.02:0.99;

% Ideal observer performance
params = Model.newModelParams('model', 'ideal', 'trials', 10000);
[~, ideal_fig] = Model.plotCategorySensorySpace(ps, ps, params);
figureToPanel(ideal_fig, fig3, 2, 4, 5, parula);

% Sampling model with gamma = 0.1
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.1, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-.4 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.99 0.61; 0.85 0.65; 0.73 0.73; 0.67 0.83; 0.65 0.95];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, fig3, 2, 4, 2, parula);
figureToPanel(pk_fig, fig3, 2, 4, 3);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, fig3, 2, 4, 4, cmap);

% Variational model with gamma = 0.1
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0.1, 'noise', 0, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
beta_range = [-.4 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.93 0.61; 0.79 0.67; 0.69 0.75; 0.65 0.85; 0.63 0.97];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, fig3, 2, 4, 6, parula);
figureToPanel(pk_fig, fig3, 2, 4, 7);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, fig3, 2, 4, 8, cmap);

%% Supplemental figure model results (showing effect of gamma parameter)

figModelSupp = figure;

% Range of values for category and sensory information
ps = 0.51:0.02:0.99;

% --- SAMPLING ---

% Replicate sampling model results above (gamma = 0)
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-0.5 0]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [.99 .59; .83 .65; .69 .79; .67 .91];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 1, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 2);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 3, cmap);

% Same with gamma = 0.1
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.1, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-.43 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.99 0.61; 0.85 0.65; 0.73 0.73; 0.67 0.83; 0.65 0.95];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 4, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 5);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 6, cmap);

% Same with gamma = 0.2
params = Model.newModelParams('model', 'is', 'var_x', 0.1, 'gamma', 0.2, 'noise', 0, 'trials', 10000, 'updates', 5, 'samples', 5);
beta_range = [-.4 .23]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.99 0.61; 0.87 0.65; 0.75 0.73; 0.67 0.83; 0.65 0.95];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 7, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 8);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 9, cmap);

% --- VB-CZX ---

% Replicate vb model results above (gamma = 0)
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0, 'noise', 0, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
beta_range = [-0.5 0]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.99 0.59; 0.81 0.65; 0.71 0.73; 0.67 0.83; 0.65 0.97];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 10, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 11);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 12, cmap);

% Same with gamma = 0.1
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0.1, 'noise', 0, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
beta_range = [-.4 .1]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.93 0.61; 0.79 0.67; 0.69 0.75; 0.65 0.85; 0.63 0.97];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 13, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 14);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 15, cmap);

% Same with gamma = 0.2
params = Model.newModelParams('model', 'vb-czx', 'var_x', 0.1, 'gamma', 0.2, 'noise', 0, 'trials', 10000, 'updates', 5, 'step_size', 0.05);
beta_range = [-.26 .2]; % min and max beta expected (to get maximum use of colorbar range)
sens_cat_pts = [0.99 0.59; 0.83 0.65; 0.71 0.73; 0.63 0.87; 0.61 0.99];
[cs_fig, pk_fig] = Model.plotCSPK(ps, ps, params, [0 0 0], 'beta', beta_range, sens_cat_pts);
figureToPanel(cs_fig, figModelSupp, 6, 3, 16, parula);
figureToPanel(pk_fig, figModelSupp, 6, 3, 17);
[~,~,~,~,fig,cmap] = Model.plotCSSlopes(ps, ps, params, beta_range, THRESHOLD, sens_cat_pts);
figureToPanel(fig, figModelSupp, 6, 3, 18, cmap);

%% Helper function for figure layout

function figureToPanel(figSource, figDest, subM, subN, subI, cmap)
margin = 0.02;
widths = (1 - (subN + 1) * margin) / subN;
heights = (1 - (subM + 1) * margin) / subM;
[col, row] = ind2sub([subN, subM], subI);
left = margin + (col-1) * (widths + margin);
top = margin + (row-1) * (heights + margin);
bottom = 1 - top - heights;
figure(figSource);
ax = gca;
ax_copy = copyobj(ax, figDest);
ax_copy.Position = [left bottom widths heights];
ax_copy.Parent = figDest;
close(figSource);
if exist('cmap', 'var')
    colormap(ax_copy, cmap);
    colorbar('peer', ax_copy);
end
end