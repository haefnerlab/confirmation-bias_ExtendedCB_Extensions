function [best_hprs] = combined_hprs_search_peripheral_stimulusCase2(subjectID, image_array_pre, ratios, contrasts, expt_type, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds, dir)
% initialize
datadir = fullfile(pwd, dir);
[num_sub,~] = size(subjectID);
trial_types = 5; % this is fixed from PreGen
for sub =1:num_sub
    data = LoadAllSubjectData(subjectID{sub}, expt_type, datadir);
    if expt_type==1
        exp = '_Ratio';
    elseif expt_type==2
        exp = '_Noise';
    end
    signal_raw = [];
    choice_raw = [];
    noises = [];
    
    filename = [datadir '/' subjectID{sub} exp];
    frame_signal_filename = [filename '.mat'];
    if exist(frame_signal_filename)
        temp = load(filename);
        frame_signals = temp.frame_signals;
        if abs(size(frame_signals,1)-length(data.choice))>0
            [frame_signals] = ComputeFrameSignals_bigstimCase2(data, image_array_pre, ratios, contrasts, trial_types);
            save(filename, 'frame_signals');
        end
    else
        [frame_signals] = ComputeFrameSignals_bigstimCase2(data, image_array_pre, ratios, contrasts, trial_types);
        save(filename, 'frame_signals');
    end
    for k = 1:length(data.choice)
        noises = [noises data.noise(k)];
        signal_raw = [signal_raw; frame_signals(k, :)];
        choice_raw = [choice_raw data.choice(k)];
    end
    
    trials = size(choice_raw, 2);
    disp(['Computing likelihoods for Subject ' num2str(sub)]);
    [~, log_likelihoods(sub,:,:,:,:)] = CustomRegression.xValidatePKwithlapseSabya(signal_raw, choice_raw, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds);
    log_likelihood_summed(sub,:,:,:) = mean(log_likelihoods(sub,:,:,:,:),5);
    disp(['Done computing for ' num2str(sub) '/' num2str(num_sub) ' subjects...']);
end
sz = size(log_likelihoods);
avg_ll =  mean(log_likelihood_summed,1);
[~, imax] = max(avg_ll(:));
[iRidge, iAR1, iCurve] = ind2sub(sz(2:4), imax);
% Err on the side of less regularization by choosing smoothing that is one order of magnitude less than the best.
% iRidge = max(iRidge-1, 1);
% iAR1 = max(iAR1-1, 1);
% iCurve = max(iCurve-1, 1);
best_hprs = [hpr_ridge(iRidge), hpr_ar1(iAR1), hpr_curvature(iCurve)];
disp('Searching of hyperparameters complete!!');
disp (['Best hyperparameters across subjects for peripheral stimulus is: ' num2str(best_hprs)]);
end