function [frame_signals, faulty_signals] = ComputeFrameSignals_bigstim(GaborData,correct_signals_pre, faulty_signals_pre, ratios, contrasts, tr_val)
% tr_val SHOULD be 10
% tr_val = 10 for Final Case and 5 for Case 2;
for trial=1:GaborData.current_trial
%     if (trial==1 || mod(trial,10)==0)
%         disp(['Trials:' num2str(trial) '/' num2str(trials)]);
%     end
    indx_ratio = find(ratios==GaborData.ratio(trial));
    indx_noise = find(GaborData.kappa_set==GaborData.noise(trial));
    indx_contrast = find(contrasts==GaborData.contrast(trial));
    for i=1:tr_val
        loc = compute_indx(i-1, GaborData.correct_answer(trial), indx_ratio-1, indx_noise-1, indx_contrast-1, 2, size(ratios,1), size(GaborData.kappa_set,2), size(contrasts,1));
        if all(abs(faulty_signals_pre{loc}' - GaborData.ideal_frame_signals(trial,:)) < 1e-10)
            frame_signals(trial,:) = correct_signals_pre{loc};
            faulty_signals(trial,:) = faulty_signals_pre{loc};
            break;
        end
    end
end
    function num  = compute_indx(i1, i2, i3, i4, i5, d2, d3, d4, d5)
        num = i1*d2*d3*d4*d5 + i2*d3*d4*d5 + i3*d4*d5 + i4*d5 + i5 + 1;
    end
end
