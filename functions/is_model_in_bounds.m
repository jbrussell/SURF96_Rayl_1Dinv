function [is_acceptable] = is_model_in_bounds(mod,bounds)
% Check that model is acceptable within bounds of M

is_acceptable = 1; % initiate flag

nlays = size(mod,1);
for ilay = 1:nlays
    if mod(ilay,3)<bounds(ilay,1) || mod(ilay,3)>bounds(ilay,2)
        is_acceptable = 0;
        return
    end
end


end

