% read eigenfunction file written by sdpder96
%
% jbrussell 8/2023
% 
function eig = read_eig96(filename,wavetype,varargin)

if nargin==2
    nmode = 0; % default, fundamental model
elseif nargin==3
	nmode = varargin{1};
end

fid = fopen(filename,'r');
if(fid>=0)	
	it = 0;
    ip = 0;
	while ~feof(fid)
    	stemp=fgetl(fid);
	    if (contains(stemp,'MODE'))    % find where eigenfunctions begin
        	C = textscan(stemp,'%s WAVE      MODE #  %d' );
            rayl_love = C{1};
            mode_k = C{2};
            
            stemp=fgetl(fid);
            C = textscan(stemp,'T = %f C =    %f U   = %f');
            period_k = C{1};
            phv_k = C{2};
            grv_k = C{3};
            
            stemp=fgetl(fid);
            C = textscan(stemp,'AR= %f GAMMA= %f ZREF= %f');
            AR_k = C{1};
            gamma_k = C{2};
            zref_k = C{3};
            
            if mode_k ~= nmode
            % wrong mode, keep going...
                continue
            end
            
            ip = ip + 1; % increment period
            
             % start reading
             
            stemp = fgetl(fid);	% skip a line	
	
            il=1;
            while ~feof(fid) % when line hasn't hit the dash line
                stemp = fgetl(fid);		
                temp = sscanf(stemp,'%f');

                if ~isempty(temp) 
                    if(wavetype=='L')
%                         M       UT         TT       DC/DH      DC/DB      DC/DR
                        eig.ut(il,ip)  = temp(2);% Horizontal displacement
                        
                        eig.tt(il,ip)  = temp(3);% Horizontal stress

                    elseif(wavetype=='R')
        %                 M       UR         TR        UZ         TZ        DC/DH      DC/DA      DC/DB      DC/DR
                        eig.ur(il,ip)  = temp(2);% Horizontal displacement
                        eig.uz(il,ip) = temp(4);% Vertical displacement
                        
                        eig.tr(il,ip)  = temp(3);% Horizontal stress				 
                        eig.tz(il,ip) = temp(5);% Vertical stress
               
                    end
                    
                    eig.mode = mode_k;
                    eig.periods(ip) = period_k;

                    il=il+1;
                else
                    break
                end
            end
        end
        it = it + 1;
        if it > 1e6
            kernel = [];
            return
        end
	end
	fclose(fid);
else
	error([filename,' does not exist']);% return empty data if file is empty 
	eig=[];
return
end

