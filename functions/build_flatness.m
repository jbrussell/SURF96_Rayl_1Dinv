function [ J00 ] = build_flatness( nlayer )
% Build first derivative "flatness" matrix
%
% jbrussell 6/5/2020

J00 = 1*eye(nlayer);
Jdown = -1*[zeros(1,nlayer); eye(nlayer-1) zeros(nlayer-1,1) ];
J00 = J00+Jdown;
J00(1,:) = 0; J00(end,:) = 0;

end

