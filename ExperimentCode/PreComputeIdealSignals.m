
big_matched_all = [1 0];
for i=1:length(big_matched_all)
    newExperimentStimComputation;
    
    big_matched = big_matched_all(i); %1 if cortex area matched, 0 if cortex area not matched -1 if small stim
    blocks = 15;
    contrast = 25;
    trials_pre_gen = 10;
    
    if big_matched==1
        stim_size = big_pixels_out;
        annulus_size = big_pixels_in;
        stim_sf = sf_big_cycles_per_pixel;
    elseif big_matched==0
        stim_size = big_pixels_out;
        annulus_size = pixels_in;%pixels_in;
        stim_sf = sf_big_cycles_per_pixel1;
    else
        stim_size = pixels_out;
        annulus_size = pixels_in;
    end
    proto = load('bpgFinaltest-subject04-Session9-GaborDataNoiseOld.mat');
    newparams = newGaborData(proto.GaborData,'stim_size',stim_size,'annulus',annulus_size,'stim_sp_freq_cpp',stim_sf,'blocks',blocks,'eyelink_use',0,'contrast',contrast);
    ratio = unique(proto.GaborData.ratio)';
    correct_answer = [0 1]';
    contrast = [newparams.contrast(1)]';
    kappa = newparams.kappa_set';
    tr_val = trials_pre_gen;
    if big_matched==1
        
        frame_array_pre = load('PreGenFrameCat_big_matched.mat');
        frame_array_pre = frame_array_pre.frame_categories;
        ratios = load('PreGenRatios_big_matched.mat');
        ratios = ratios.ratio;
        contrasts = load('PreGenContrasts_big_matched.mat');
        contrasts = contrasts.contrast;
        image_array_pre = load('PreGenStim_big_matched.mat');
        image_array_pre = image_array_pre.image_array_pre;
        
        [correct_signals,faulty_signals] = CreateSignalsGaborStimulus_pre_gen(newparams,image_array_pre,ratio,correct_answer,contrast,kappa,tr_val);
        
        save('PreGenStim_big_matchedTrueSignals','correct_signals');
        save('PreGenStim_big_matchedFaultySignals','faulty_signals');
        
    else
        
        frame_array_pre = load('PreGenFrameCat_big.mat');
        frame_array_pre = frame_array_pre.frame_categories;
        ratios = load('PreGenRatios_big.mat');
        ratios = ratios.ratio;
        contrasts = load('PreGenContrasts_big.mat');
        contrasts = contrasts.contrast;
        image_array_pre = load('PreGenStim_big.mat');
        image_array_pre = image_array_pre.image_array_pre;
        
        [correct_signals,faulty_signals] = CreateSignalsGaborStimulus_pre_gen(newparams,image_array_pre,ratio,correct_answer,contrast,kappa,tr_val);
        
        save('PreGenStim_bigTrueSignals','correct_signals');
        save('PreGenStim_bigFaultySignals','faulty_signals');
        
        
    end
    
end
