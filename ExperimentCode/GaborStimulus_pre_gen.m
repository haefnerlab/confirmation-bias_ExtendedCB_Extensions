function [image_array_pre, frame_categories] = GaborStimulus_pre_gen(GaborData,ratio,correct_answer,contrast,kappa,trials)
%GABORSTIMULUS(GaborData, trial) create (or recreate) stimulus frames based
%on parameters in GaborData and the seed, contrast, ratio, and noise on the
%given trial. If 'GaborData.iid(trial)' is true, each frame's category is
%drawn iid based on the 'ratio' parameter. Otherwise, exactly
%round(ratio*num_images) frames will match the 'true' category.
%
%This function makes no modifications to GaborData.

% Set RNG state to recreate stimulus for this trail.
% rng(GaborData.seed(trial), 'twister');
stim_fcn = @bpg.genImages;

% image_array_pre = zeros(trials,size(correct_answer,1),size(ratio,1),size(kappa,1),size(contrast,1),GaborData.number_of_images,GaborData.stim_size,GaborData.stim_size);
% frame_categories = zeros(trials,size(correct_answer,1),size(ratio,1),size(kappa,1),size(contrast,1),GaborData.number_of_images);
% disp(size(image_array_pre));
k=1;
for i=1:trials
    disp(i)
    for c_ans=1:size(correct_answer,1)
        for ra=1:size(ratio,1)
            for kp=1:size(kappa,1)
                for con=1:size(contrast,1)
                    rng(GaborData.seed(k), 'twister');
                    n_match = round(ratio(ra) * GaborData.number_of_images);
                    match_frames = [true(1, n_match) false(1, GaborData.number_of_images - n_match)];
                    match_frames = Shuffle(match_frames);
                    
                    % Choose frames based on whether correct answer this trial is Left or Right
                    if correct_answer(c_ans) == 1
                        frame_categories{k}(match_frames) = GaborData.left_category;
                        frame_categories{k}(~match_frames) = GaborData.right_category;
                    else
                        frame_categories{k}(~match_frames) = GaborData.left_category;
                        frame_categories{k}(match_frames) = GaborData.right_category;
                    end
                    
                    % Set random seed again to keep match_frames independent of pixel noise.
                    rng(GaborData.seed(k), 'twister');
                    [image_array,~,~,~] = stim_fcn(GaborData.number_of_images, GaborData.stim_size, ...
                        GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp, ...
                        squeeze(frame_categories{k}(:)), kappa(kp), GaborData.annulus);
                    
                    image_array = uint8(image_array * contrast + 127);
                    
                    image_array = min(image_array, 255);
                    image_array = max(image_array, 0);
                    image_array_pre{k} = image_array;
                    info_array(k,:) = [i correct_answer(c_ans) ratio(ra) kappa(kp) contrast(con)];
                    
                    k = k+1;
                end
            end
        end
        
    end
end
end