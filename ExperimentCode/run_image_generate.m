function [] = run_image_generate(stim_size,annulus,stim_sp,blocks,contrast,big_matched,trials_pre_gen)
proto = load('bpgFinaltest-subject04-Session9-GaborDataNoiseOld.mat');
newparams = newGaborData(proto.GaborData,'stim_size',stim_size,'annulus',annulus,'stim_sp_freq_cpp',stim_sp,'blocks',blocks,'contrast',contrast);
ratio = unique(proto.GaborData.ratio)';
correct_answer = [0 1]';
contrast = [newparams.contrast(1)]';
kappa = newparams.kappa_set';
tr_val = trials_pre_gen;
[image_array_pre, frame_categories] = GaborStimulus_pre_gen(newparams,ratio,correct_answer,contrast,kappa,tr_val);
[faulty_signals, correct_signals] = CreateSignalsGaborStimulus_pre_gen(newparams,image_array_pre,ratio,correct_answer,contrast,kappa,tr_val);
%%
whos image_array_pre;
whos faulty_signals;
whos correct_signals;
%%
% if big_matched==1
%     save('PreGenStim_big_matched','image_array_pre','-v7.3');
%     save('PreGenFrameCat_big_matched','frame_categories');
%     save('PreGenRatios_big_matched','ratio');
%     save('PreGenContrasts_big_matched','contrast');
% else
%     save('PreGenStim_big','image_array_pre','-v7.3');
%     save('PreGenFrameCat_big','frame_categories');
%     save('PreGenRatios_big','ratio');
%     save('PreGenContrasts_big','contrast');
% end
str_ending = datestr(now);
if big_matched==1
    save(['PreGenStim_big_matched' str_ending],'image_array_pre','-v7.3');
    save(['PreGenStim_big_matchedTrueSignals' str_ending],'correct_signals','-v7.3');
    save(['PreGenStim_big_matchedFaultySignals' str_ending],'faulty_signals','-v7.3');
    save(['PreGenFrameCat_big_matched' str_ending],'frame_categories');
    save(['PreGenRatios_big_matched' str_ending],'ratio');
    save(['PreGenContrasts_big_matched' str_ending],'contrast');
    save('FilenameEnding','str_ending');
else
    save(['PreGenStim_big' str_ending],'image_array_pre','-v7.3');
    save(['PreGenStim_bigTrueSignals' str_ending],'correct_signals','-v7.3');
    save(['PreGenStim_bigFaultySignals' str_ending],'faulty_signals','-v7.3');
    save(['PreGenFrameCat_big' str_ending],'frame_categories');
    save(['PreGenRatios_big' str_ending],'ratio');
    save(['PreGenContrasts_big' str_ending],'contrast');
    save('FilenameEnding','str_ending');
end
end