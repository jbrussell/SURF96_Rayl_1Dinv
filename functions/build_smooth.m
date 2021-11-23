function [ F00 ] = build_smooth( nlayer )
% Second derivative smoothing matrix
%
% jbrussell 6/5/2020

F00 = 2*eye(1*nlayer);
Fup = -1*[zeros(1*nlayer-1,1) eye(1*nlayer-1) ;zeros(1,1*nlayer) ];
Fdown = -1*[zeros(1,1*nlayer); eye(1*nlayer-1) zeros(1*nlayer-1,1) ];
F00 = F00+Fup+Fdown;
F00(1,:) = 0; F00(end,:) = 0;

end

