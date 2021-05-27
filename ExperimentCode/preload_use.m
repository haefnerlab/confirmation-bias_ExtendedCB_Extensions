if big_matched==1
    
    frame_array_pre = load('PreGenFrameCat_big_matched.mat');
    frame_array_pre = frame_array_pre.frame_categories;
    ratios = load('PreGenRatios_big_matched.mat');
    ratios = ratios.ratio;
    contrasts = load('PreGenContrasts_big_matched.mat');
    contrasts = contrasts.contrast;
    image_array_pre = load('PreGenStim_big_matched.mat');
    image_array_pre = image_array_pre.image_array_pre;
    correct_signals_pre = load('PreGenStim_big_matchedTrueSignals.mat');
    correct_signals_pre = correct_signals_pre.correct_signals;
    faulty_signals_pre = load('PreGenStim_big_matchedFaultySignals.mat');
    faulty_signals_pre = faulty_signals_pre.faulty_signals;
    
else
    
    frame_array_pre = load('PreGenFrameCat_big.mat');
    frame_array_pre = frame_array_pre.frame_categories;
    ratios = load('PreGenRatios_big.mat');
    ratios = ratios.ratio;
    contrasts = load('PreGenContrasts_big.mat');
    contrasts = contrasts.contrast;
    image_array_pre = load('PreGenStim_big.mat');
    image_array_pre = image_array_pre.image_array_pre;
    correct_signals_pre = load('PreGenStim_big_TrueSignals.mat');
    correct_signals_pre = correct_signals_pre.correct_signals;
    faulty_signals_pre = load('PreGenStim_big_FaultySignals.mat');
    faulty_signals_pre = faulty_signals_pre.faulty_signals;
    
end
