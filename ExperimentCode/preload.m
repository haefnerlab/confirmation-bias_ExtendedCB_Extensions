if big_matched==1
    str_end = load('FilenameEnding.mat');
    str_ending = str_end.str_ending;
    frame_array_pre = load(['PreGenFrameCat_big_matched' str_ending '.mat']);
    frame_array_pre = frame_array_pre.frame_categories;
    ratios = load(['PreGenRatios_big_matched.mat' str_ending '.mat']);
    ratios = ratios.ratio;
    contrasts = load(['PreGenContrasts_big_matched' str_ending '.mat']);
    contrasts = contrasts.contrast;
    image_array_pre = load(['PreGenStim_big_matched' str_ending '.mat']);
    image_array_pre = image_array_pre.image_array_pre;
    correct_signals_pre = load(['PreGenStim_big_matchedTrueSignals' str_ending '.mat']);
    correct_signals_pre = correct_signals_pre.correct_signals_pre;
    faulty_signals_pre = load(['PreGenStim_big_matchedFaultySignals' str_ending '.mat']);
    faulty_signals_pre = faulty_signals_pre.faulty_signals_pre;
else
    str_end = load('FilenameEnding.mat');
    str_ending = str_end.str_ending;
    frame_array_pre = load(['PreGenFrameCat_big' str_ending '.mat']);
    frame_array_pre = frame_array_pre.frame_categories;
    ratios = load(['PreGenRatios_big' str_ending '.mat']);
    ratios = ratios.ratio;
    contrasts = load(['PreGenContrasts_big' str_ending '.mat']);
    contrasts = contrasts.contrast;
    image_array_pre = load(['PreGenStim_big' str_ending '.mat']);
    image_array_pre = image_array_pre.image_array_pre;
    correct_signals_pre = load(['PreGenStim_bigTrueSignals' str_ending '.mat']);
    correct_signals_pre = correct_signals_pre.correct_signals_pre;
    faulty_signals_pre = load(['PreGenStim_bigFaultySignals' str_ending '.mat']);
    faulty_signals_pre = faulty_signals_pre.faulty_signals_pre; 
end