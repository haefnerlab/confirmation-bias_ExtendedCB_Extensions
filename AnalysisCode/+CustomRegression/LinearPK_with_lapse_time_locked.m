function [sobl, postVal] = LinearPK_with_lapse_time_locked(data, responses, time, standardize)
%LinearPK Regress PK as a linear function.
%
% sob = LinearPK(data, responses) returns sob = [slope offset bias].

if nargin < 3, standardize = 0; end

% Standardize each regressor.
switch standardize
    case 0
        % do nothing
    case 1
        % assume 0 mean (nothing to subtact) and iid (std taken over all data)
        data = data / std(data(:));
    case 2
        data = zscore(data);
    otherwise
        error('Expected argument ''standardize'' to be one of [0, 1, 2]');
end

% convert boolean to float type
% assert(islogical(responses));
responses = 1.0 * responses(:);

[~, frames] = size(data);

    function NLL = neg_bernoulli_log_likelihood(sobl)
        weights = sobl(2) + ((0:frames-1) * time) * sobl(1);
        logits = data * weights(:) + sobl(3);
        lapse = 1e-4+(1-1e-4) * sigmoid(sobl(4));%sobl(4)^2;
        neg_log_bernoulli = -1*(log(0.5*lapse+(1-lapse)*sigmoid(logits(:))).*responses)-1*(log(1-0.5*lapse-(1-lapse)*sigmoid(logits(:))).*(1-responses));
        NLL = sum(neg_log_bernoulli);
    end

options = optimset('MaxIter', 1000000, 'Display', 'off', 'MaxFunEvals', 1000000);
[sobl, negPostVal] = fminunc(@neg_bernoulli_log_likelihood, zeros(1,4), options);
postVal = -negPostVal;
end
