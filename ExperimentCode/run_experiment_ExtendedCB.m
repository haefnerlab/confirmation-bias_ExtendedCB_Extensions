% clear all;close all;clc;
newExperimentStimComputation;

big_matched = 0; %1 if cortex area matched, 0 if cortex area not matched -1 if small stim
image_generated = 1;
blocks = 15;
track_eye = 1;
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
if big_matched==-1
    proto = load('bpgFinaltest-subject04-Session9-GaborDataNoiseOld.mat');
    newparams = newGaborData(proto.GaborData,'stim_size',stim_size,'annulus',annulus_size,'blocks',blocks,'eyelink_use',track_eye,'contrast',contrast);
    ExperimentGabor(newparams);
else
    if image_generated==0
        run_image_generate(stim_size,annulus_size,stim_sf,blocks,contrast,big_matched,trials_pre_gen);
    end
    preload_old
    proto = load('bpgFinaltest-subject04-Session9-GaborDataNoiseOld.mat');
    newparams = newGaborData(proto.GaborData,'stim_size',stim_size,'annulus',annulus_size,'stim_sp_freq_cpp',stim_sf,'blocks',blocks,'eyelink_use',track_eye,'contrast',contrast);
    ExperimentGabor_pre_gen(newparams,image_array_pre,frame_array_pre, correct_signals_pre, trials_pre_gen, ratios, contrasts);
end