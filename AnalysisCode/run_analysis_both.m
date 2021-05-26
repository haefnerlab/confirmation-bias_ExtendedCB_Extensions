function [params_boot,sobl,abbl,trials,bin_centers,means,stderrs,data,log_bernoulli] = run_analysis_both(subjectID, expt_type, boot_n, best_hprs, dir)

% initialize
datadir = fullfile(pwd, dir);

data = LoadAllSubjectData(subjectID,expt_type,datadir);
disp('Data loaded for the subject!');
if expt_type==1
    exp = '_Ratio';
elseif expt_type==2
    exp = '_Noise';
end
signal_raw = [];
choice_raw = [];

filename = [datadir '/' subjectID exp];
frame_signal_filename = [filename '.mat'];
if exist(frame_signal_filename)
    temp = load(filename);
    frame_signals = temp.frame_signals;
    if abs(size(frame_signals,1)-length(data.choice))>0
        [frame_signals, ~, ~] = ComputeFrameSignals(data, 0);
        save(filename, 'frame_signals');
    end
else
    [frame_signals, ~, ~] = ComputeFrameSignals(data, 0);
    save(filename, 'frame_signals');
end

if (expt_type == 2)
    for k = 1:length(data.choice)
        signal_raw = [signal_raw; frame_signals(k, :)];
        choice_raw = [choice_raw data.choice(k)];
    end
end
if (expt_type == 1)
    for k = 1:length(data.choice)
        signal_raw = [signal_raw; frame_signals(k, :)];
        choice_raw = [choice_raw data.choice(k)];
    end
end
trials = size(choice_raw, 2);
disp('Starting parameters/kernels estimation ...');
for j = 1:boot_n
    if (j==1 || mod(j,100)==0)
        disp(['Computing ' num2str(j) '/' num2str(boot_n) ' bootstraps ...']);
    end
    [signal, choice] = bootstrap(signal_raw, choice_raw,trials);
    [sobl(j,:), ~] = CustomRegression.LinearPK_with_lapse(signal, choice, 0);
    [params_boot(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapseSabya(signal, choice, best_hprs(1), 0, best_hprs(3), 0);%,hprs, 0, hprs, 1);
    [abbl(j,:), ~, ~] = CustomRegression.ExponentialPK_with_lapse(signal, choice, 0);
end
all_frames = size(signal_raw,2);
temporal_kernel = prctile(params_boot(:, 1:all_frames), 50);
bias =  prctile(params_boot(:, end-1), 50);
disp('Getting log odds...');
[log_bernoulli] = compute_log_odds(signal_raw, temporal_kernel, bias);

    function [logits] = compute_log_odds(data,weights,bias_computed)
        logits = data * weights(:) + bias_computed;
    end

    function [signal, choice] = bootstrap(signal_raw, choice_raw,trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signal = [];
        choice = [];
        for i = 1:trials
            trial_num = sample_nums(i);
            signal = [signal; signal_raw(trial_num, :)];
            choice = [choice choice_raw(trial_num)];
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