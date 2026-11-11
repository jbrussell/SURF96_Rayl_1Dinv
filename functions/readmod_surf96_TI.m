% Read a transversely isotropic (TI) velocity model from a CPS MODEL96 file
% (TI analog of readmod_surf96.m)
%
% model = readmod_surf96_TI(filename)
%
% model is nlayer x 9 matrix [H Vpv Vsv rho Vph Vsh eta Qp Qs]
%       model(:,1): H(KM)
%       model(:,2): VPV(KM/S)
%       model(:,3): VSV(KM/S)
%       model(:,4): DENSITY(g/cm^3)
%       model(:,5): VPH(KM/S)
%       model(:,6): VSH(KM/S)
%       model(:,7): ETA = F/(A-2L)  (computed from VPF in the file)
%       model(:,8): QP
%       model(:,9): QS
% Use model(:,1:7) for the elastic part only.
%
% If the file is an ISOTROPIC MODEL96 file it is read and converted to TI
% (Vph=Vpv=Vp, Vsh=Vsv=Vs, eta=1).
%
% jbrussell 11/2026
%
function model = readmod_surf96_TI(filename)

fid = fopen(filename,'r');
if(fid<0)
	error([filename,' does not exist']);
end

lines = {};
stemp = fgetl(fid);
while ischar(stemp)
	lines{end+1} = stemp; %#ok<AGROW>
	stemp = fgetl(fid);
end
fclose(fid);

istiso = any(~cellfun(@isempty,strfind(upper(lines(1:min(end,4))),'TRANSVERSE')));

% data begin after the 'H(KM)' header line (and for TI the VPH header line)
ihdr = find(~cellfun(@isempty,strfind(upper(lines),'H(KM)')),1);
if isempty(ihdr)
	ihdr = 12;
end
istart = ihdr + 1;
if istiso && istart <= numel(lines) && ~isempty(strfind(upper(lines{istart}),'VPH'))
	istart = istart + 1;
end

vals = [];
for ii = istart:numel(lines)
	t = sscanf(lines{ii},'%f');
	vals = [vals; t(:)]; %#ok<AGROW>
end

if istiso
	nper = 13; % H VPV VSV RHO QP QS ETAP ETAS FREFP FREFS VPH VSH VPF
else
	nper = 10; % H VP VS RHO QP QS ETAP ETAS FREFP FREFS
end
nlay = floor(numel(vals)/nper);
vals = reshape(vals(1:nlay*nper),nper,nlay)';

model = zeros(nlay,9);
model(:,1) = vals(:,1);
model(:,2) = vals(:,2);
model(:,3) = vals(:,3);
model(:,4) = vals(:,4);
model(:,8) = vals(:,5);
model(:,9) = vals(:,6);
if istiso
	model(:,5) = vals(:,11);
	model(:,6) = vals(:,12);
	vpf = vals(:,13);
	denom = model(:,5).^2 - 2*model(:,3).^2;
	model(:,7) = vpf.^2 ./ denom;
	model(denom==0,7) = 1;
else
	model(:,5) = vals(:,2);
	model(:,6) = vals(:,3);
	model(:,7) = 1;
end

end
