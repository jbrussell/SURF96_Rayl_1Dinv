function [phV,grV,phVq] = readMINEOS_qfile2(qfile,periods,mode)
% [phV,grV] = readMINEOS_qfile(qfile,swperiods)
%  
%  Function to read MINEOS qfile file (qfile) with the fundamental mode
%  phase (and group) velocities listed by period, and then interpolate to
%  get velocities at the desired periods (swperiods).

fid = fopen(qfile,'r');
line=fgetl(fid); % read line with # of lines in q model
dum=str2num(line); %#ok<ST2NM>
nQline=dum(1); % will skip this many lines
C = textscan(fid,'%d %d %f %f %f %f %f %f %f %f','Headerlines',nQline); % n,l,omega,Q,?,c,U
fclose(fid);

% find mode of interest
I_mode = find(C{1}==mode);

l = C{2}(I_mode); % mode degree
freq_all = C{3}(I_mode)/2/pi; % all the frequencies
Q = C{4}(I_mode); % Q at each frequency
phV_all = C{6}(I_mode); % phase velocity at each frequency
grV_all = C{7}(I_mode); % group velocity at each frequency
phV_all_q = C{8}(I_mode);
freq_all_q = 1./C{9}(I_mode);
T = C{10}(I_mode);

% desired frequencies
freq_want = 1./periods;

% interpolate to get velocities
phV = linterp(freq_all,phV_all,freq_want);
grV = linterp(freq_all,grV_all,freq_want);

phVq = linterp(freq_all_q,phV_all_q,freq_want);


end


function [ YI ] = linterp(X,Y,XI)
% YI = LINTERP(X,Y,XI) interpolates to find YI, the values of the
%     underlying function Y at the points in the array XI. X must be a
%     vector of length N.
% this function differs from the simple interp1 matlab function in that it
% can accept a vector x with multiple values of y (e.g. at the top of one
% layer and the bottom of another)
%
% Z. Eilon   May 2015
X = X(:); Y = Y(:); XI = XI(:);

YI = zeros(size(XI));

% use find_ilay to find where things fit - will not work properly for when
% a member of X matches a member of XI
[ilay] = find_ilay(XI,X);

% linearly interpolate
YI = (XI - X(ilay)).*(Y(ilay+1)-Y(ilay))./(X(ilay+1)-X(ilay)) + Y(ilay);
% now sort out coincident elements
olap = intersect(X,XI);
for i = 1:length(olap)
YI(XI==olap(i)) = mean(Y(X==olap(i)));
end

end


function [ilay] = find_ilay(r,Rb)
% [ilay] = find_ilay(r,Rb)
%
% Function to find the indices that each element of vector r would slot
% into vector Rb - originally conceived as a solution to the problem of
% having a series points at different radii and wanting to know which
% layers each of them were in, where the boundaries of the layers
% (including the top and bottom) are given by Rb. 
% 
% For example, if 3 layers were given by boundaries: [0;100;200;300] the
% point 50 would be in layer 1, and 217 would be in layer 3.
% thus, find_ilay([50;217],[0;100;200;300]) = [1;3]
%
% If a point in r is on a boundary, it is put into the upper layer, unless
% it is right at the max, then it is included in outermost layer

if any(r < min(Rb)) || any(r > max(Rb))
    error('r must be within extremes of Rb')
end

r = r(:);
Rb = Rb(:);

N = length(r);
Nlay = length(Rb);

[~,ilay] = max((ones(N,1)*Rb' - r*ones(1,Nlay))>0,[],2);
ilay = ilay-1;
ilay(ilay==0) = Nlay-1; % if any are on outer edge, say in outermost layer

end

