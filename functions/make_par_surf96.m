% make parameter file sobs.d for surf96 inversion/ forward calculations 
% function fid=make_par_surf96(datatype),
% datatype:
%    R: rayleigh
%    L: Love
%    J: Joint
%
% jbrussell 11/20/2022: Include ability to calculate overtones.
% nmode=0 for fund. mode; nmode=1 for 1st overtone, etc...
%
function fid=make_par_surf96(datatype,varargin)

if nargin == 1
    nmode = 0; % default fundamental mode
elseif nargin == 2
    nmode = varargin{1};
end
num_modes_to_calc = nmode + 1;

	if (datatype=='R') 
		% Rayleigh wave
		fid=fopen('sobs.d','w');
		fprintf(fid,'0.005 '); %df
		fprintf(fid,'0.005 '); % dcr
		fprintf(fid,'0. ');
		fprintf(fid,'0.005 '); % dcl
		fprintf(fid,'0. \n'); % dcr
		fprintf(fid,'0 '); % error based on redisual(1) or std of data(0)  
		fprintf(fid,'0  0  0  0  %d  %d',num_modes_to_calc,num_modes_to_calc); %  if use gam_love/c_love/u_love/gama_ray/c_ray/u_ray
		fprintf(fid,' 0  1  0\n'); % don't know what they mean
		fprintf(fid,'start.mod'); % name of starting model file
		fprintf(fid,'\n'); 
		fprintf(fid,'disp_obs.dsp');% name of dispersion data file
		fprintf(fid,'\n');
		fclose(fid);
	elseif (datatype=='L')
		% love wave
		fid=fopen('sobs.d','w');
		fprintf(fid,'0.005 '); %df
		fprintf(fid,'0.005 '); % dcr
		fprintf(fid,'0. ');
		fprintf(fid,'0.005 '); % dcl
		fprintf(fid,'0. \n'); % dcr
		fprintf(fid,'0 '); % error based on redisual(1) or std of data(0)  
		fprintf(fid,'0  %d  %d  0  0  0',num_modes_to_calc,num_modes_to_calc); %  if use gam_love/c_love/u_love/gama_ray/c_ray/u_ray
		fprintf(fid,' 0  1  0\n'); % don't know what they mean
		fprintf(fid,'start.mod'); % name of starting model file
		fprintf(fid,'\n'); 
		fprintf(fid,'disp_obs.dsp');% name of dispersion data file
		fprintf(fid,'\n');
		fclose(fid);
	elseif  (datatype=='J')
		% Joint inversion
		fid=fopen('sobs.d','w');
		fprintf(fid,'0.005 '); %df
		fprintf(fid,'0.005 '); % dcr
		fprintf(fid,'0. ');
		fprintf(fid,'0.005 '); % dcl
		fprintf(fid,'0. \n'); % dcr
		fprintf(fid,'0 '); % error based on redisual(1) or std of data(0)  
		fprintf(fid,'0  %d  %d  0  %d  %d',num_modes_to_calc,num_modes_to_calc,num_modes_to_calc,num_modes_to_calc); %  if use gam_love/c_love/u_love/gama_ray/c_ray/u_ray
		fprintf(fid,' 0  1  0\n'); % don't know what they mean
		fprintf(fid,'start.mod'); % name of starting model file
		fprintf(fid,'\n'); 
		fprintf(fid,'disp_obs.dsp');% name of dispersion data file
		fprintf(fid,'\n');
		fclose(fid);
	end
return
end



	
		


