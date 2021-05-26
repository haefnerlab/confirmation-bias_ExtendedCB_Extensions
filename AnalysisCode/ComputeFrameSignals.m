function [frame_signals, checksum, faulty_signals] = ComputeFrameSignals(GaborData, kernelKappa, spFreqCPP, spFreqStdCPP)
%COMPUTEFRAMESIGNALS wrapper around bpg.getSignal that computes the full [trials x frames] matrix of
%signal levels for the given data struct.
%
% If kernelKappa is 0, the per-trial 'true' value of kappa is used.
%
% spFreqCPP and spFreqStdCPP arguments are optional.

trials = length(GaborData.choice);

frame_signals = zeros(trials, GaborData.number_of_images);
faulty_signals = zeros(trials, GaborData.number_of_images);
checksum = zeros(1,trials);

if nargin == 2
    args = {kernelKappa};
elseif nargin == 4
    args = {kernelKappa, spFreqCPP, spFreqStdCPP};
else
    error('Expected 2 or 4 args');
end
if isfield(GaborData, 'wedge_chosen')
    GaborData.use_wedge_code = 1;
end
% NOTE: to use parfor and rng requires rng to set the generator type
% explicitly!
for t=1:trials
    if (t==1 || mod(t,10)==0)
        disp(['Trials:' num2str(t) '/' num2str(trials)]);
    end
    % Regenerate stimulus from seed (can be slow if error-recovery is needed)
    image_array = GaborStimulus_regenerate(GaborData, t, 2);
    checksum(t) = mean(image_array(:));
    
    args_copy = args;
    if kernelKappa == 0
        args_copy{1} = max(0.04, GaborData.noise(t));
    end
    faulty_signals(t, :) = ...
        bpg.getSignal(image_array - 127, GaborData.left_category, args_copy{:}) - ...
        bpg.getSignal(image_array - 127, GaborData.right_category, args_copy{:});
    
    % Center the image after converting out of uint8 type
    image_array = double(image_array)-127;
    frame_signals(t,:) = ...
        bpg.getSignal(image_array, GaborData.left_category, args_copy{:}) - ...
        bpg.getSignal(image_array, GaborData.right_category, args_copy{:});
end

end