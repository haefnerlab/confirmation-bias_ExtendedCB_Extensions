clear all;close all;clc;
dir = 'RawData';

    subjects = {...
        'bpgshorterframes-subject02';...
%         'bpgshorterframes-subject04';
%         'bpgshorterframes-subject05';
%         'bpgshorterframes-subject03';
%         'bpgshorterframes-subject11';
%         'bpgshorterframes-subject08';%colleen
%         'bpgshorterframes-subject10';
%         'bpgshorterframes-subject07';
%         'bpgshorterframes-subject14';
%         'bpgshorterframes-subject13'%eimi
        };
[num_sub,~] = size(subjects);
expt_type = 2;
signal_stored = [];
signal_corrected = [];
signal_regenerated = [];
checksum_stored = [];
checksum_corrected = [];
checksum_regenerated = [];
for sub=1:num_sub
    subjectID = subjects{sub};
    datadir = ['/Users/achattoraj/Desktop/Projects/CB/confirmation-bias_ExtendedCB_Extensions/' dir];%fullfile(pwd, dir);
    data = LoadAllSubjectData(subjectID,expt_type,datadir);
    disp('Data loaded!');
    signal_stored = [signal_stored; data.ideal_frame_signals(:)];
    [frame_signals, checksum, faulty_signals] = ComputeFrameSignals(data, 0);
    signal_corrected = [signal_corrected; frame_signals(:)];
    signal_regenerated = [signal_regenerated; faulty_signals(:)];
    checksum_regenerated = [checksum_regenerated; checksum(:)];
    checksum_stored = [checksum_stored; data.checksum(:)];
end

%%
figure();
vals1 = linspace(min(min(checksum_regenerated),min(checksum_stored)),max(max(checksum_regenerated),max(checksum_stored)),100);
subplot(1,3,1)
scatter(checksum_stored,checksum_regenerated,'o');
hold on;
plot(vals1,vals1,'k');
xlabel('Saved mean of image array per trial')
ylabel('Recovered mean of image array per trial')
title('Matching mean of image array')

vals2 = linspace(min(min(signal_regenerated),min(signal_stored)),max(max(signal_regenerated),max(signal_stored)),100);
subplot(1,3,2)
scatter(signal_stored, signal_regenerated,'o');%signal_regenerated,'o');
hold on;
plot(vals2,vals2,'k');
xlabel('Saved faulty ideal signals')
ylabel('Recovered faulty ideal signals')
title('Matching faulty signals')

vals3 = linspace(min(min(signal_corrected),min(signal_stored)),max(max(signal_corrected),max(signal_stored)),100);
subplot(1,3,3)
scatter(signal_stored,signal_corrected,'o');
hold on;
plot(vals3,vals3,'k');
xlabel('Saved faulty signals')
ylabel('Recovered correct signals')
title('Comparing faulty stored and correct recovered signals')




