function [correct_signals,faulty_signals] = CreateSignalsGaborStimulus_pre_gen(GaborData,image_array_pre,ratio,correct_answer,contrast,kappa,trials)
%GABORSTIMULUS(GaborData, trial) create (or recreate) stimulus frames based
%on parameters in GaborData and the seed, contrast, ratio, and noise on the
%given trial. If 'GaborData.iid(trial)' is true, each frame's category is
%drawn iid based on the 'ratio' parameter. Otherwise, exactly
%round(ratio*num_images) frames will match the 'true' category.
%
%This function makes no modifications to GaborData.

% Set RNG state to recreate stimulus for this trail.
% rng(GaborData.seed(trial), 'twister');

% image_array_pre = zeros(trials,size(correct_answer,1),size(ratio,1),size(kappa,1),size(contrast,1),GaborData.number_of_images,GaborData.stim_size,GaborData.stim_size);
% frame_categories = zeros(trials,size(correct_answer,1),size(ratio,1),size(kappa,1),size(contrast,1),GaborData.number_of_images);
% disp(size(image_array_pre));
k=1;
for i=1:trials
    for c_ans=1:size(correct_answer,1)
        for ra=1:size(ratio,1)
            for kp=1:size(kappa,1)
                for con=1:size(contrast,1)
                    im = image_array_pre{k};
                    correct_signals{k} = ...
                        bpg.getSignal(double(im) - 127, GaborData.left_category, max(kappa(kp), .04)) - ...
                        bpg.getSignal(double(im) - 127, GaborData.right_category, max(kappa(kp), .04));
                    faulty_signals{k} = ...
                        bpg.getSignal(uint8(im) - 127, GaborData.left_category, max(kappa(kp), .04)) - ...
                        bpg.getSignal(uint8(im) - 127, GaborData.right_category, max(kappa(kp), .04));
                    disp([i k])
                    k = k+1;
                end
            end
        end
        
    end
end
end