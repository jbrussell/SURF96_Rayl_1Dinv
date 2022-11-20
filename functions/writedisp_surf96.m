%function to write dispersion data to surf96 format file
% INPUT:
% data:nfreq*3 matrix: [vec_T, vec_v,vec_err]
% wavetype: 'R' for Rayleigh, L for Love; default is R
% veltype: 'C' for phase v and 'U' for group V; default is C
%
% Example:
% overwrite = 1;
% nmode = 0;
% writedisp_surf96(fakedata,'disp_obs.dsp','R','C',overwrite,nmode);
%
% Created by Yang Zha 09/28/2012
%
function fid=writedisp_surf96(data,filename,varargin)

if(nargin<2)
	error('need filename and data input!')
	exit -1
end
% determine wvetype and phase or group velocity
%
if(nargin==2)
	wavetype='R';
	veltype='C';
	overwrite=1;
    nmode=0;
elseif(nargin==3)
	wavetype=varargin{1};
	veltype='C';
	overwrite=1;
    nmode=0;
elseif(nargin==4)
	wavetype=varargin{1};
	veltype=varargin{2};
	overwrite=1;
    nmode=0;
elseif(nargin==5)
	wavetype=varargin{1};
	veltype=varargin{2};	
	overwrite=varargin{3};
    nmode=0;
elseif(nargin==6)
	wavetype=varargin{1};
	veltype=varargin{2};	
	overwrite=varargin{3};
    nmode=varargin{4};
end

% write data
%
if(overwrite)
	fid=fopen(filename,'w'); % open a file and overwrite (DEFAULT)
else
	fid=fopen(filename,'a'); % open a file and attach to end of exsiting files 
end	

[nt,nd]=size(data);
for it=1:nt
	fprintf(fid, 'SURF96 %s %s X   %g   %7.4f   %7.4f   %7.4f\n',...
    wavetype,veltype,nmode,data(it,1),data(it,2),data(it,3)); % assume same format as SURF96 output

end
fclose(fid);

return
end
