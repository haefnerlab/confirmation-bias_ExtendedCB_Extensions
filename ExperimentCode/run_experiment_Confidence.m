clear all;close all;clc;

phase = 2;
image_generated = 0;
blocks = 15;
track_eye = 0;
contrast = 25;
confidence_report = 1;

if phase==2
    proto = load('bpgFinaltest-subject04-Session9-GaborDataNoiseOld.mat');
    newparams = newGaborData(proto.GaborData,'blocks',blocks,'eyelink_use',track_eye,'contrast',contrast,'confidence_report',confidence_report,'stair_fn', @Staircase.noise);
    ExperimentGabor(newparams);
elseif phase==1
    proto = load('bpgFinaltest-subject04-Session9-GaborDataNoiseOld.mat');
    newparams = newGaborData(proto.GaborData,'blocks',blocks,'eyelink_use',track_eye,'contrast',contrast,'confidence_report',confidence_report,'stair_fn', @Staircase.ratio);
    ExperimentGabor(newparams);
end