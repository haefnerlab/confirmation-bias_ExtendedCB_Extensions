function [image_properties, tracker_info, broke_condition, quit] = trialStimuliGabor(Data, image_array, wPtr, tracker_info, settings)
% trialStimuli displays the animation of several gabor patches in quick
% succession to the subject, or runs through a single trial of the
% experiment.
if tracker_info.eyelink_use
    tracker_info = Eyelink.startTrialPrepare(tracker_info);
end
image_properties = struct();
quit = false;
broke_condition = false;
frame_duration = 1 / Data.stimulus_fps;


%% Make sure to have left/right patch to match the orientations used

% Create images to be displayed as left or right options
left_patch = uint8(squeeze(bpg.genImages(1, Data.stim_size, Data.stim_sp_freq_cpp, Data.stim_std_sp_freq_cpp, Data.left_category, 2)) * 64.0 + 127);
right_patch = uint8(squeeze(bpg.genImages(1, Data.stim_size, Data.stim_sp_freq_cpp, Data.stim_std_sp_freq_cpp, Data.right_category, 2)) * 64.0 + 127);

left_patch = min(left_patch, 255);
left_patch = max(left_patch, 0);
right_patch = min(right_patch, 255);
right_patch = max(right_patch, 0);


xc = settings.screenSize(3)/2; % Get the middle of the horizontal axis
yc = settings.screenSize(4)/2; % Get the middle of the vertical axis

black = [0 0 0];
gray = [127 127 127];

% Set up variables for keyboard functions
KbName('UnifyKeyNames');
exitKey = KbName(settings.keyExit);
leftKey = KbName(settings.keyLeft);
rightKey = KbName(settings.keyRight);
temp_length=0;

total_frames = Data.number_of_images + 1;

image_texture = zeros(1, total_frames);
for i = 1:Data.number_of_images
    image_texture(i) = Screen('MakeTexture', wPtr, squeeze(image_array(i, :, :)));
end
[~, h, w] = size(image_array);
noise_mask = uint8(squeeze(bpg.genImages(1, Data.stim_size, Data.stim_sp_freq_cpp, Data.stim_std_sp_freq_cpp, 0, 0, Data.annulus)) * Data.contrast(Data.current_trial)  + 127);
noise_mask = min(noise_mask, 255);
noise_mask = max(noise_mask, 0);
image_texture(end) = Screen('MakeTexture', wPtr, noise_mask);

show_left_patch = Screen('MakeTexture', wPtr, left_patch);
show_right_patch = Screen('MakeTexture', wPtr, right_patch);

stimulus_bbox = ptbCenteredRect([xc yc], [w h]);

Screen('FillRect', wPtr, gray);  % Make the background gray

Eyelink.drawFixationSymbol(tracker_info,xc,yc, wPtr);
[~, stimOnsetTime] = Screen('Flip', wPtr);
if tracker_info.eyelink_use
    tracker_info = Eyelink.startTrial(tracker_info);
    Eyelink.clearBuffer(tracker_info.drained);
end
% Get fixation (takes a variable amount of time).
% [tracker_info, broke_fixation, ~, temp_length] = Eyelink.track_eye(tracker_info,temp_length,[xc yc], tracker_info.fixationMinimumHold, 1,'is_fixating');
% % If the fixation broke, the trial aborts. Otherwise, continues.
% if broke_fixation==1
%     broke_condition = 1;
%     image_properties.choice=nan;
%     return;
% end

% Draw frame around where stimulus will appear as a timing cue (note:
% leaving fixation cue on the screen).
Screen('FillRect', wPtr, gray);
% drawTrialNo();
% Preparing for eye-tracking protocol when eyelink_use is true.

Eyelink.drawFixationSymbol(tracker_info, xc, yc, wPtr);
Screen('FrameOval', wPtr, black, stimulus_bbox);
[~, cueOnsetTime] = Screen('Flip', wPtr);
[tracker_info, broke_fixation, ~, temp_length ] = Eyelink.track_eye(tracker_info, temp_length, [xc yc], Data.cue_duration, 1, 'is_fixating');
if broke_fixation==1
    broke_condition = 1;
    image_properties.choice=nan;
    return;
end
% Prep for first stimulus frame by clearing the drawStimulusFrame.
Screen('FillRect', wPtr, gray);
nextStimTime = cueOnsetTime + Data.cue_duration;

% Present each image for 'frame_duration' seconds.
% TODO - track eyes at full temporal resolution rather than once per frame.
for i = 1:total_frames
    % Show stimulus.
    Screen('DrawTexture', wPtr, image_texture(i), [], stimulus_bbox); %Fill the buffer with the first texture
    Eyelink.drawFixationSymbol(tracker_info, xc,yc,wPtr);
    [~, stimOnsetTime] = Screen('Flip', wPtr, nextStimTime);
    [tracker_info, broke_fixation, ~, temp_length ] = Eyelink.track_eye(tracker_info, temp_length, [xc yc], frame_duration - Data.blank_duration, 1, 'is_fixating');
    if broke_fixation==1
        broke_condition = 1;
        image_properties.choice = nan;
        return;
    end
    nextStimTime = stimOnsetTime + frame_duration;
    
    % (Maybe) end stimulus frame with some blank frames.
    if Data.blank_duration > 0
        Screen('FillRect', wPtr, gray, stimulus_bbox);
        Eyelink.drawFixationSymbol(tracker_info,xc,yc, wPtr);
        Screen('Flip', wPtr, stimOnsetTime + frame_duration - Data.blank_duration);
        [tracker_info, broke_fixation, ~, temp_length ] = Eyelink.track_eye(tracker_info, temp_length, [xc yc], Data.blank_duration, 1, 'is_fixating');
        if broke_fixation==1
            broke_condition = 1;
            image_properties.choice=nan;
            return;
        end
    end
end

Screen('FillRect', wPtr, gray, stimulus_bbox);
Eyelink.drawFixationSymbol(tracker_info,xc,yc ,wPtr);
[~, endTime] = Screen('Flip', wPtr);
[tracker_info, broke_fixation, ~, ~ ] = Eyelink.track_eye(tracker_info, temp_length, [xc yc], Data.go_cue_time, 1, 'is_fixating');
if broke_fixation==1
    broke_condition = 1;
    image_properties.choice=nan;
    return;
end
Screen('Flip', wPtr, endTime + Data.go_cue_time);

Screen('DrawTexture', wPtr, show_left_patch, [], ptbCenteredRect([xc-w yc], [w h]));   % xc, yc indicates the coordinates of the middle of the screen
Screen('DrawTexture', wPtr, show_right_patch, [], ptbCenteredRect([xc+w yc], [w h]));
%Screen('DrawText', wPtr, sprintf('Current Trial - #%d', Data.current_trial), xc-600, yc+250, 0);   % Unobtrusive output to screen of the current trial number
Screen('Flip', wPtr);

[key, rt, timeout] = ptbWaitKey([leftKey, rightKey, exitKey], 1);

% Close textures to avoid memory problems.
for i = 1:total_frames
    Screen('Close', image_texture(i));
end
Screen('Close', show_left_patch);
Screen('Close', show_right_patch);

if key == exitKey
    quit = true;
    image_properties.choice = nan;
end

if timeout
    image_properties.choice = nan;
else
    image_properties.reaction = rt * 1000;
    if key == leftKey
        image_properties.choice = 1;
    elseif key == rightKey
        image_properties.choice = 0;
    end
end

    function drawTrialNo()
        Screen('DrawText', wPtr, sprintf('Current Trial - #%d', Data.current_trial), xc-900, yc+550, 0);
    end
end