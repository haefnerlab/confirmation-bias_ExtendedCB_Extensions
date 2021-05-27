function [image_array,frame_categories, ideal_frame_signals, loc] = ChooseStimulus(GaborData,image_array_pre, frame_array_pre, correct_signals_pre, ratios, contrasts, trial, tr_val)

rndm = randi(tr_val);
% 
% for rr=1:size(ratios,1)
%     if GaborData.ratio(trial)==ratios(rr)
%         indx_ratio = rr;
%     end
% end
% for con=1:size(contrasts,1)
%     if GaborData.contrast(trial)==contrasts(con)
%         indx_contrast = con;
%     end
% end
indx_ratio = find(ratios==GaborData.ratio(trial));
indx_noise = find(GaborData.kappa_set==GaborData.noise(trial));
indx_contrast = find(contrasts==GaborData.contrast(trial));
tic
loc = compute_indx(rndm-1, GaborData.correct_answer(trial), indx_ratio-1, indx_noise-1, indx_contrast-1, 2, size(ratios,1), size(GaborData.kappa_set,2), size(contrasts,1));
image_array = squeeze(image_array_pre{loc}(:,:,:));
frame_categories = squeeze(frame_array_pre{loc}(:));
ideal_frame_signals = correct_signals_pre{loc};
toc

    function num  = compute_indx(i1, i2, i3, i4, i5, d2, d3, d4, d5)
        num = i1*d2*d3*d4*d5 + i2*d3*d4*d5 + i3*d4*d5 + i4*d5 + i5 + 1;
    end
end