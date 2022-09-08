function [sample] = sample_model(priors,N,Ncoeffs)

% Function to draw model from prior
sample = nan(N,Ncoeffs);
for ic = 1:Ncoeffs
    sample(:,ic) = priors(N,ic);
end

end

