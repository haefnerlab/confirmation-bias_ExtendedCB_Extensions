% function [im, imF, filterF, aperture, w_chosen] = genImages(frames, width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, annulusPix)%,w_chosen)
function [im, imF, filterF, aperture] = genImages_wedge(frames, width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa, annulusPix,w_chosen)

%BPG.GENIMAGES Create a sequence band-pass grating (bpg) stimuli.
%
%[im, imF] = BPG.GENIMAGES(frames, width, spFreqCPP, spFreqStdCPP, oriDEG, oriKappa)
% creates [frames x width x width] array of images. spFreqCPP sets the mean
% spatial frequency in cycles per pixel. spFreqStdCPP sets the range of
% spatial frequencies present. oriDEG sets the mean rotation, oriKappa
% sets the range of orientation energy present.
%
% oriDEG may be a vector of orientations of length 'frames'.

noise = randn(frames, width, width);
noiseF = framefun(@(f) fftshift(fft2(f)), noise);

[~, ~, rho, theta] = freq_coords(width);

% Create separate [-1, 1] range meshgrid for pixel-space filters.
[px, py] = meshgrid(linspace(-1, 1, width));
pr = sqrt(px.^2 + py.^2);

if length(oriDEG) == 1, oriDEG = oriDEG * ones(1, frames); end

im = zeros(frames, width, width);
imF = zeros(frames, width, width);

for i=1:2
    wedge(i,:,:) = wedge_block(width,i);
end
%% Create spatial frequency filter
spFreqFilter = pdf('rician', rho / width, spFreqCPP, spFreqStdCPP);

%% Create gaussian aperture
aperture = exp(-4 * pr.^2);
if nargin >= 7 && annulusPix > 0
    % Cut out annulus hole.
    aperture = aperture .* (1 + erf(10 * (pr - annulusPix / width)));
end
% w_chosen = randi([1,2],1);
%% Generate each frame.
for f=1:frames
    % Create orientation filters for each frame. Note that 'theta' is
    % doubled to create two symmetric filters in the Fourier domain
    % (bow-tie rather than cone shape). 'oriDEG' must also be doubled to
    % compensate.
    oriFilter = bpg.vmpdf(2 * theta, 2 * deg2rad(oriDEG(f)), oriKappa);
    oriFilter(isnan(oriFilter)) = 0;
    
    % Get full, normalized foureir-domain filter.
    filterF = spFreqFilter .* oriFilter;
    filterF = filterF / sum(filterF(:));
    
    % Apply fourier-domain filters on each frame.
    imF(f, :, :) = squeeze(noiseF(f, :, :)) .* filterF;
    if frames>1
        im(f, :, :) = aperture .* (real(ifft2(ifftshift(squeeze(imF(f, :, :)))))) .* squeeze(wedge(w_chosen,:,:));
    else
        im(f, :, :) = aperture .* (real(ifft2(ifftshift(squeeze(imF(f, :, :))))));
        
    end
end

%% Normalize range in pixel space to +/- 1
im = im / max(abs(im(:)));
end

function frames = framefun(fn ,frames)
%Helper to apply fn to each frame in frames.
for f=1:size(frames, 1)
    frames(f, :, :) = fn(squeeze(frames(f, :, :)));
end
end