function [params_boot,params_boot_first_half,params_boot_second_half,...
    sobl,sobl_first_half,sobl_second_half,...
    abbl,abbl_first_half,abbl_second_half,...
    trials,bin_centers,means,stderrs,data,frame_signals,...
    sobl_time_locked,sobl_time_locked_first_half,sobl_time_locked_second_half,...
    log_bernoulli] = run_analysis_noise_only(subjectID, expt_type, time, boot_n, best_hprs, standardize, dir)

% initialize
datadir = fullfile(pwd, dir);

data = LoadAllSubjectData(subjectID,expt_type,datadir);
disp('Data for the subject loaded!');
if expt_type==1
    exp = '_Ratio';
elseif expt_type==2
    exp = '_Noise';
end
signal_raw = [];
choice_raw = [];
noises = [];
signal_raw_first_half = [];
choice_raw_first_half = [];
noises_first_half = [];
signal_raw_second_half = [];
choice_raw_second_half = [];
noises_second_half = [];
filename = [datadir '/' subjectID exp];
frame_signal_filename = [filename '.mat'];
if exist(frame_signal_filename)
    temp = load(filename);
    frame_signals = temp.frame_signals;
    if abs(size(frame_signals,1)-length(data.choice))>0
        [frame_signals] = ComputeFrameSignals(data, 0);
        save(filename, 'frame_signals');
    end
else
    [frame_signals] = ComputeFrameSignals(data,0);
    save(filename, 'frame_signals');
end
mark_half = floor(length(data.choice)/2);
data.ideal_frame_signals = frame_signals;
for k = 1:length(data.choice)
    noises = [noises data.noise(k)];
    signal_raw = [signal_raw; frame_signals(k, :)];
    choice_raw = [choice_raw data.choice(k)];
    if k<=mark_half
        noises_first_half = [noises_first_half data.noise(k)];
        signal_raw_first_half = [signal_raw_first_half; frame_signals(k, :)];
        choice_raw_first_half = [choice_raw_first_half data.choice(k)];
    else
        noises_second_half = [noises_second_half data.noise(k)];
        signal_raw_second_half = [signal_raw_second_half; frame_signals(k, :)];
        choice_raw_second_half = [choice_raw_second_half data.choice(k)];
    end
end
trials = size(choice_raw, 2);
disp('Preprocessing of data completed!');
for j = 1:boot_n
    if j==1 || mod(j,100)==0
        disp(['Completed ' num2str(j) '/' num2str(boot_n) ' iterations...']);
    end
    [signal, choice] = bootstrap(signal_raw, choice_raw, trials);
    [signal_first_half, choice_first_half] = bootstrap(signal_raw_first_half, choice_raw_first_half, mark_half);
    [signal_second_half, choice_second_half] = bootstrap(signal_raw_second_half, choice_raw_second_half, (trials-mark_half));
    
    [sobl(j,:), ~] = CustomRegression.LinearPK_with_lapse(signal, choice, standardize);
    [sobl_time_locked(j,:), ~] = CustomRegression.LinearPK_with_lapse_time_locked(signal, choice, time, standardize);
    [params_boot(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapseSabya(signal, choice, best_hprs(1), best_hprs(2), best_hprs(3), standardize);
    [abbl(j,:), ~, ~] = CustomRegression.ExponentialPK_with_lapse(signal, choice, standardize);

    [sobl_first_half(j,:), ~] = CustomRegression.LinearPK_with_lapse(signal_first_half, choice_first_half, standardize);
    [sobl_time_locked_first_half(j,:), ~] = CustomRegression.LinearPK_with_lapse_time_locked(signal_first_half, choice_first_half, time, standardize);
    [params_boot_first_half(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapseSabya(signal_first_half, choice_first_half, best_hprs(1), best_hprs(2), best_hprs(3), standardize);
    [abbl_first_half(j,:), ~, ~] = CustomRegression.ExponentialPK_with_lapse(signal_first_half, choice_first_half, standardize);

    [sobl_second_half(j,:), ~] = CustomRegression.LinearPK_with_lapse(signal_second_half, choice_second_half, standardize);
    [sobl_time_locked_second_half(j,:), ~] = CustomRegression.LinearPK_with_lapse_time_locked(signal_second_half, choice_second_half, time, standardize);
    [params_boot_second_half(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapseSabya(signal_second_half, choice_second_half, best_hprs(1), best_hprs(2), best_hprs(3), standardize);
    [abbl_second_half(j,:), ~, ~] = CustomRegression.ExponentialPK_with_lapse(signal_second_half, choice_second_half, standardize);
end
all_frames = size(signal_raw,2);
temporal_kernel = prctile(params_boot(:, 1:all_frames), 50);
bias =  prctile(params_boot(:, end-1), 50);
disp('Getting log odds...');
[log_bernoulli] = compute_log_odds(signal_raw, temporal_kernel, bias);

    function [logits] = compute_log_odds(data,weights,bias_computed)
        logits = data * weights(:) + bias_computed;
    end
    function [signals, choices] = bootstrap(signals_raw, choices_raw,trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signals = [];
        choices = [];
        for i = 1:trials
            trial_num = sample_nums(i);
            signals = [signals; signals_raw(trial_num, :)];
            choices = [choices choices_raw(trial_num)];
        end
    end
bin_edges = linspace(min(mean(signal_raw,2)), max(mean(signal_raw,2)), 11);
bin_halfwidth = (bin_edges(2) - bin_edges(1)) / 2;
bin_centers = bin_edges(1:end-1) + bin_halfwidth;
means = zeros(size(bin_centers));
stderrs = zeros(size(bin_centers));
for b=1:length(bin_centers)
    % Select all points for which bin i is closest.
    bin_dists = abs(mean(signal_raw,2) - bin_centers(b));
    indices = bin_dists <= bin_halfwidth;
    means(b) = mean(data.accuracy(indices));
    stderrs(b) = std(data.accuracy(indices)) / sqrt(sum(indices));
end

end