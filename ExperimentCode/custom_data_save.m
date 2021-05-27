datadir = 'C:\Users\HaefnerLab\Documents\MATLAB\Summer\RawData';
fileName = [datadir '\bpgFinaltest-subject17-Session1-GaborDataNoiseQuit.mat'];
contents = load(fileName);
point = 120;
GaborData = TruncateQuitDataGabor_data_save(contents.GaborData,point);
save(fileName, 'GaborData');