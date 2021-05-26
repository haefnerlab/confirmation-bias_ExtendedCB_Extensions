function [ AIC_separate, AIC_comb,LL_separate,LL_comb, params_boot1,best_hprs1,params_boot2,best_hprs2,params_boot3,best_hprs3 ] = return_AIC(subjectID,expt_type,boot_n,hpr1,hpr2)
% subjectID
case1 = subjectID(1);
case2 = subjectID(2);
%%
data1 = LoadAllSubjectData(case1,expt_type(1));
signal_raw1 = [];
choice_raw1 = [];

for k = 1:size(data1.ideal_frame_signals, 1)
    signal_raw1 = [signal_raw1; data1.ideal_frame_signals(k, :)];
    choice_raw1 = [choice_raw1 data1.choice(k)];
end

%%
data2 = LoadAllSubjectData(case2,expt_type(2));
signal_raw2 = [];
choice_raw2 = [];

for k = 1:size(data2.ideal_frame_signals, 1)
    signal_raw2 = [signal_raw2; data2.ideal_frame_signals(k, :)];
    choice_raw2 = [choice_raw2 data2.choice(k)];
end


%%
% data3 = ConcatGaborData(data1, data2);
% signal_raw3 = [];
% choice_raw3 = [];
% 
% for k = 1:size(data3.ideal_frame_signals, 1)
%     signal_raw3 = [signal_raw3; data3.ideal_frame_signals(k, :)];
%     choice_raw3 = [choice_raw3 data3.choice(k)];
% end

%%
trials1 = size(choice_raw1, 2);
trials2 = size(choice_raw2, 2);
for bt_num = 1:boot_n
    [signal1, choice1] = bootstrap(signal_raw1, choice_raw1,trials1);
    [signal2, choice2] = bootstrap(signal_raw2, choice_raw2,trials2);
    signal3 = [signal1' signal2']';
    choice3 = [choice1 choice2];
    [params_boot1(bt_num,:),best_hprs1(bt_num,:),LL1(bt_num)] = neglog_likelihood_for_PK(signal1,choice1,hpr1,hpr2);
    disp('done')
    [params_boot2(bt_num,:),best_hprs2(bt_num,:),LL2(bt_num)] = neglog_likelihood_for_PK(signal2,choice2,hpr1,hpr2);
    disp('done')
    [params_boot3(bt_num,:),best_hprs3(bt_num,:),LL3(bt_num)] = neglog_likelihood_for_PK(signal3,choice3,hpr1,hpr2);   
    disp('done')
end
    function [signal, choice] = bootstrap(signal_raw, choice_raw,trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signal = [];
        choice = [];
        for tr = 1:trials
            trial_num = sample_nums(tr);
            signal = [signal; signal_raw(trial_num, :)];
            choice = [choice choice_raw(trial_num)];
        end
    end
LL_separate = LL1 + LL2;
LL_comb = LL3;
AIC_separate = 2 * 24 + 2 * LL_separate;
AIC_comb = 2 * 12 + 2 * LL_comb;
end
