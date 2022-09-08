function [ spbasis,spwts,spzz ] = make_splines( zknots,dz,zi,vi )
%  [ spbasis,spwts,spzz ] = make_splines( zknots,dz,zi,vi )
%   
% Function to make splines within a current model with knots at absolute
% depths zknots. The splines will be defined on an interpolated grid (spzz)
% between the top and bottom splines. We will find weights for these
% splines by interpolating the basis onto vectors zi and vi.
% 
% INPUTS:
%   zknots   - vector of depths of spline knots (km)
%   dz.mod.zi  - increment of depth grid (km)
%   zi   - vector of depths for input velocity profile (km)
%   vi   - vector of velocities for input velocity profile [optional] (km/s)
% OUTPUTS:
%   spbasis  - matrix with spline basis (Nz x Nsp)
%   spwts    - (Nsp x 1) vector of spline weights to fit [zi,vi] profile 
%   spzz     - (Nz x 1) vector of depths for splines (km)
% 
% N.B.  to move spline knots but not change velocities, just ignore the new
% weightings, use new spbasis and knots positions. 

if nargin<4 % no velocities input for interpolation, so just use dummies.
    vi = zi;
end

minz = zknots(1);
maxz = zknots(end);

allzknots = [repmat(minz,3,1);zknots;repmat(maxz,3,1)];

sp = fastBSpline.lsqspline(allzknots,2,zi,vi); % interpolate onto current model

% spzz = unique([minz:dz:maxz,maxz])';
spzz = zi;
spbasis = sp.getBasis(spzz); 
spbasis = spbasis(:,2:end-1);                

if nargin<4
    spwts = [];
else
    spwts = sp.weights(2:end-1); % pull back out spline coeff's from the interpolation        
end


end

