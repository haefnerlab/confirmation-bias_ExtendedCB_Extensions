function [lpo, x_samples1, weights1,x_samples2, weights2] = isLogOddsUpdate_parallel(params, e, lpo)
%MODEL.ISLOGODDSUPDATE compute update to log odds of C using importance sampling model.

trials = size(e, 1);
oz = ones(trials, 1);

updates = params.updates;
noise = params.noise;
gamma = params.gamma;
sig_s = sqrt(params.var_s);
sig_x = sqrt(params.var_x);
p_match = params.p_match;
samples = floor(params.samples/4);
cognitive_load = params.cognitive_load;

% Create two distributions representing p(x|C=+1) and p(x|C=-1) in the generative model
p_x_Cp1 = mog.create([+1 -1], [sig_x sig_x], [p_match 1-p_match]);
p_x_Cm1 = mog.create([-1 +1], [sig_x sig_x], [p_match 1-p_match]);

p_x_Cp2 = mog.create([+1 -1], [sig_x sig_x], [p_match 1-p_match]);
p_x_Cm2 = mog.create([-1 +1], [sig_x sig_x], [p_match 1-p_match]);

x_samples1 = zeros(trials, samples, updates);
weights1 = zeros(trials, samples, updates);

x_samples2 = zeros(trials, samples, updates);
weights2 = zeros(trials, samples, updates);

for n=1:updates
    % Convert from lpo (log odds) to the probability that C=+1
    pC =    (1 ./ (1 + exp(-lpo)));
    
    % Create likelihoods. Format is a matrix where each row specifies triples of [mu, sigma, pi] of
    % a mixture of Gaussians. Only one mode in the likelihood, but it's useful to use the MoG
    % format. See @mog.create
    
    likelihoods1 = [e(:) sig_s*oz oz];
    likelihoods2 = [e(:) sig_s*oz oz];
    
    
    % Create the prior on x by marginalizing over the current posterior of C. The prior is also a
    % mixture of gaussians, but with 2 modes corresponding to C = +/-1
    p =  ((p_match * pC + (1 - p_match) * (1 - pC)));
    sig_x = sig_x + cognitive_load;
    priors = [+oz, sig_x*oz, p, -oz, sig_x*oz, (1-p)];
    
    % Q is the distribution from which samples of x are drawn; it is the current estimate of the
    % posterior over x using lpo as the prior over C
    Q1 = mog.prod(likelihoods1, priors);
    Q2 = mog.prod(likelihoods2, priors);
    
    % Draw samples from Q
    samp1 = mog.sample(Q1, samples);
    x_samples1(:, :, n) = samp1;
    flat_samples1 = samp1(:);
    
    samp2 = mog.sample(Q2, samples);
    x_samples2(:, :, n) = samp2;
    flat_samples2 = samp2(:);
    
    % Get unnormalized importance weights for each sample (note that this is vectorized over trials,
    % but we have to loop over samples. Typically number of samples << number of trials)
    for s=1:samples
        weights1(:, s, n) = 1 ./ mog.pdf(x_samples1(:, s, n), priors);
        weights2(:, s, n) = 1 ./ mog.pdf(x_samples2(:, s, n), priors);
    end
    
    % Normalize importance-sampling weights
    if params.importance_norm
        weights1(:, :, n) = weights1(:, :, n) ./ sum(weights1(:, :, n), 2);
        weights2(:, :, n) = weights2(:, :, n) ./ sum(weights2(:, :, n), 2);
    end
    
    % Compute p(x|C=+1) and p(x|C=-1), then take weighted sum for each trial.
    pCp1 = sum(reshape(mog.pdf(flat_samples1, p_x_Cp1), [trials samples]) .* weights1(:, :, n), 2);
    pCm1 = sum(reshape(mog.pdf(flat_samples1, p_x_Cm1), [trials samples]) .* weights1(:, :, n), 2);
    
    pCp2 = sum(reshape(mog.pdf(flat_samples2, p_x_Cp2), [trials samples]) .* weights2(:, :, n), 2);
    pCm2 = sum(reshape(mog.pdf(flat_samples2, p_x_Cm2), [trials samples]) .* weights2(:, :, n), 2);
    
    % Log likelihood odds is log(pCp/pCm)
    llo1 = (log(pCp1) - log(pCm1));
    llo2 = (log(pCp2) - log(pCm2));
    
    lpo = lpo  * (1 - gamma / updates) + llo1 / updates + llo2 / updates ;
    
    % Add zero-mean additive noise.
    lpo = lpo + randn(trials, 1) * noise;
end

end