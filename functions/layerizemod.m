function [MODEL] = layerizemod(mod)
% Layerize surf96 mod file
%

THICK = mod(:,1);
[M,N] = size(mod);
L = length(THICK);

for i=1:L			%for all layers
    if i > 1
        dep(2*i-1,1) = dep(2*i-2,1);
    else
        dep(2*i-1,1) = 0;
    end
    dep(2*i,1) = dep(2*i-1,1)+THICK(i);
    for j = 1:N
        mod_lay(2*i-1,j) = mod(i,j);
        mod_lay(2*i,j) = mod(i,j);
    end
end

MODEL.z = dep;
MODEL.vp = mod_lay(:,2);
MODEL.vs = mod_lay(:,3);
MODEL.rho = mod_lay(:,4);

end
