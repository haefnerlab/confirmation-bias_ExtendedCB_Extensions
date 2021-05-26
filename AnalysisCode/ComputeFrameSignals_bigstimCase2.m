function [frame_signals] = ComputeFrameSignals_bigstimCase2(GaborData, image_array_pre, ratios, contrasts, tr_val)
% tr_val SHOULD be 5
% tr_val =  5;
for trial=1:GaborData.current_trial
    temp_sig = [];
    if (trial==1 || mod(trial,100)==0)
        disp(['Trials:' num2str(trial) '/' num2str(GaborData.current_trial)]);
    end
    indx_ratio = find(ratios==GaborData.ratio(trial));
    indx_noise = find(GaborData.kappa_set==GaborData.noise(trial));
    indx_contrast = find(contrasts==GaborData.contrast(trial));
    kernel_kappa = max(0.04, GaborData.noise(trial));
    for i=1:tr_val
        loc = compute_indx(i-1, GaborData.correct_answer(trial), indx_ratio-1, indx_noise-1, indx_contrast-1, 2, size(ratios,1), size(GaborData.kappa_set,2), size(contrasts,1));
        image_array = image_array_pre{loc};
        image_array = uint8(image_array * GaborData.contrast(trial) + 127);
        image_array = min(image_array, 255);
        image_array = max(image_array, 0);
        faulty_signals_pre = bpg.getSignal(image_array-127, GaborData.left_category, kernel_kappa) - ...
            bpg.getSignal(image_array-127, GaborData.right_category, kernel_kappa);
        if all(abs(faulty_signals_pre - GaborData.ideal_frame_signals(trial,:)) < 1e-10)
             temp_sig = bpg.getSignal(double(image_array-127), GaborData.left_category, kernel_kappa) - ...
                bpg.getSignal(double(image_array-127), GaborData.right_category, kernel_kappa);
            break;
        end
    end
    if ~isempty(temp_sig)
        frame_signals(trial,:) = temp_sig;
    else
        frame_signals(trial,:) = GaborData.ideal_frame_signals(trial,:);
    end
end
    function num  = compute_indx(i1, i2, i3, i4, i5, d2, d3, d4, d5)
        num = i1*d2*d3*d4*d5 + i2*d3*d4*d5 + i3*d4*d5 + i4*d5 + i5 + 1;
    end
end