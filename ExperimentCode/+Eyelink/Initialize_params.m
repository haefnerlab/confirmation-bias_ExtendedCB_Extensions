function tracker_info = Initialize_params(whichscreen,wPtr,varargin)

[winwidth, winheight] = WindowSize(wPtr);
resolution = Screen('Resolution', whichscreen);

tracker_info = struct('whichscreen', whichscreen, ...
    'pixelsPerGazeCoordinate', [resolution.width, resolution.height], ... % X, Y screen pixels per 'gaze unit'
    ... % parameters for getFixation()
    'fixationSymbol', '+', ... % 'r' for rect, 'c' for circle, 'b' for bullseye, or '+' for plus
    'fixationSymbolSize', [10, 10], ... % pixel size of fixation symbol, independent of the 'Rect' below
    'fixationSymbolColors', [255 255 255; 0 0 0], ... % primary/secondary color of fixation symbol
    'fixationTime', 1000, ... % ms. Max time allowed in getFixation()
    'fixationMinimumHold', 0.2, ... % Time required within fixation area to consider it held.
    ... % parameters for isFixation()
    'fixationCorrection', [0 0], ... % Add this to [gx, gy] to get corrected position (this is set automatically during getFixation)
    'fixationCenter', [resolution.width/2, resolution.height/2], ...
    'fixationRadius', 50, ... % true size for fixation requirement (separate from the symbol size above)
    'pre_fixationRadius', 120, ... % true size for fixation requirement for first few seconds of stimulus. 
    ... % parameters for calibration
    'calibration_matrix', [], ...
    'collectQueue', true, ...
    'custom_calibration', false, ...
    'custom_calibrationScale', 0.2500, ...
    'calibration_color', [127 127 127],...
    'calibrationtargetcolor' , [255 0 0],...
    'calibrationtargetsize', 30, ...
    ... % parameters for general purpose
    'saveEDF',false,...
    'eyelink_use', true,...
    'mouse_use', false,...
    'wPtr', wPtr, ...
    'display_ppd', 1, ...
    'sound_use', 1, ...
    'winRect', [0 0 winwidth winheight], ...
    'viewdist', 57, ...  % in cm
    'widthcm', 61, ... % in cm
    'heightcm', 35);


for val_idx=2:2:length(varargin)
    key = varargin{val_idx-1};
    if ~ischar(key)
        warning('invalid input to initEyeTracker. After whichscreen,wPtr all arguments should be (..., ''key'', value, ...)');
    elseif ~isfield(tracker_info, key)
        warning('unrecognized tracker_info field: ''%s''', key);
    else
        tracker_info.(key) = varargin{val_idx};
    end
end
if ~tracker_info.mouse_use
    tracker_info.fixationRadius = 1000;
    tracker_info.pre_fixationRadius = 1000;
end
end