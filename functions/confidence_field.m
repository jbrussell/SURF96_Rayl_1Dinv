function [conf_mat] = confidence_field(probability_mat)
%CONFIDENCE_FIELD
F = probability_mat;

% Sort the weighting function and take the cumsum
[FF,I] = sort(F(:),1,'descend');
FFsum = cumsum(FF);
% Load the values of FFsum into YY according to the original order of F
% given by index I
YY = ones(length(FF),1);
YY(I) = FFsum;
% Reshape YY to have the orignal form of F
YYY = reshape(YY,size(F));

conf_mat = YYY;

end

