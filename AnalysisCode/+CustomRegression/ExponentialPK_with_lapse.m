function [abbl, postVal, errors] = ExponentialPK_with_lapse(data, responses, standardize)
%EXPONENTIALPK Regress PK as an exponential function of time.
%
% abbl = LinearPK(data, responses) returns abbl = [alpha beta bias lapse] such that w=alpha*exp(beta*f)
% and choices~sigmoid(signal'*w+bias)

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

    function NLL = neg_bernoulli_log_likelihood(abbl)
        weights = abbl(1) * exp(abbl(2) * (0:frames-1));
        logits = data * weights(:) + abbl(3);
        lapse = 1e-4+(1-1e-4) * sigmoid(abbl(4));%abbl(4)^2;
        neg_log_bernoulli = -1*(log(0.5*lapse+(1-lapse)*sigmoid(logits(:))).*responses)-1*(log(1-0.5*lapse-(1-lapse)*sigmoid(logits(:))).*(1-responses));
        NLL = sum(neg_log_bernoulli);
        NLL = NLL + 1/2*abbl(2)^2/100;
    end

compute_error = nargout > 2;

glm_weights = glmfit(data, responses, 'binomial');
ab = expFit(glm_weights(2:end));
init_guess = [ab glm_weights(1) 0.0];
if sum(~isfinite(init_guess))>0 || sum(~isreal(init_guess))>0 || sum(~isnan(init_guess))>0 
    for t=1:4
       if ~isfinite(init_guess(t)) || ~isreal(init_guess(t)) || ~isnan(init_guess(t))
           init_guess(t) = 0.0;
       end
    end
end
options = optimset('MaxIter', 1000000, 'Display', 'off', 'MaxFunEvals', 1000000);

% Fit weights using 'fminunc', only computing the hessian (which is slow) if errors are requested.
if compute_error
    [abbl, negPostVal, ~, ~, ~, hessian] = fminunc(@neg_bernoulli_log_likelihood, init_guess, options);
    % attempt to invert the hessian for standard error estimate - this sometimes fails silently,
    % returning NaN.
    errors = sqrt(diag(abs(inv(-hessian))));
else
    [abbl, negPostVal] = fminunc(@neg_bernoulli_log_likelihood, init_guess, options);    
end
postVal = -negPostVal;
end
