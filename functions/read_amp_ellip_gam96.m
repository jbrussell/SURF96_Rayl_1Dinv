% read amplitude, ellipticity, and gamma
%
% jbrussell 8/2023
% 
function disper = read_amp_ellip_gam96(filename,wavetype,varargin)

if nargin==2
    nmode = 0; % default, fundamental model
elseif nargin==3
	nmode = varargin{1};
end

fid = fopen(filename,'r');
if(fid>=0)	
    % skip header line
    %  RMODE NFREQ    PERIOD(S) FREQUENCY(Hz)  C(KM/S)      U(KM/S)       AR          GAMMA(1/KM)  ELLIPTICITY
    stemp=fgetl(fid);
	ip = 0;
	while ~feof(fid)
	    stemp=fgetl(fid);
        C = textscan(stemp,'%f');
        mode_k = C{1}(1);
        if mode_k ~= nmode
            continue
        end
        
        ip = ip + 1;
        
        nfreq(ip) = C{1}(2);
        disper.period(ip) = C{1}(3);
        freq(ip) = C{1}(4);
        disper.phv(ip) = C{1}(5); % phase velocity km/s
        disper.grv(ip) = C{1}(6); % group velocity km/s
        disper.gamma(ip) = C{1}(8); % attenuation coefficient (1/km)        
        if strcmp(wavetype,'R')
            disper.A_R(ip) = C{1}(7); % Rayleigh amplification
            disper.RZ(ip) = C{1}(9); % ellipticity            
        elseif strcmp(wavetype,'L')
            disper.A_L(ip) = C{1}(7); % Love amplification
        end
	end
	fclose(fid);
else
	error([filename,' does not exist']);% return empty data if file is empty 
	eig=[];
return
end

