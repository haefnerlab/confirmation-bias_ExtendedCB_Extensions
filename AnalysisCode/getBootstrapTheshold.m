function thresh = getBootstrapTheshold(data,frame_signals,performance,boots)
signal_raw = [];
choice_raw = [];
sign_noise_raw = [];
noise_raw = [];
accuracy_raw = [];
for k = 1:size(data.ideal_frame_signals, 1)
    signal_raw = [signal_raw; frame_signals(k, :)];
    choice_raw = [choice_raw data.choice(k)];
    sign_noise_raw = [sign_noise_raw data.sign_noise(k)];
    noise_raw = [noise_raw data.noise(k)];
    accuracy_raw = [accuracy_raw data.accuracy(k)];
end
trials = size(choice_raw,2); 
for i=1:boots
    [signal, choice, sign_noise, noise, accuracy] = bootstrap(signal_raw, choice_raw, sign_noise_raw, noise_raw, accuracy_raw, trials);
    data_temp.choice = choice;
    data_temp.sign_noise = sign_noise;
    data_temp.ideal_frame_signals = signal;
    data_temp.noise = noise;
    data_temp.accuracy = accuracy;
    [pm_fit,~,~,~] = GaborPsychometric(data_temp, 2);
    thresh(i) = getThreshold(pm_fit,performance); 
end
    function [signal, choice, sign_noise,noise,accuracy] = bootstrap(signal_raw, choice_raw, sign_noise_raw, noise_raw, accuracy_raw, trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signal = [];
        choice = [];
        sign_noise = [];
        noise = [];
        accuracy = [];
        for j = 1:trials
            trial_num = sample_nums(j);
            signal = [signal; signal_raw(trial_num, :)];
            choice = [choice choice_raw(trial_num)];
            sign_noise = [sign_noise sign_noise_raw(trial_num)];
            noise = [noise noise_raw(trial_num)];
            accuracy = [accuracy accuracy_raw(trial_num)];
        end
    end
end