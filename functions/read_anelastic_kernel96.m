% read surface wave sensitivity kernels wrote by srfker96 command
%	kernel=[dcdvs(:) dcdvp(:) dudvs(:) dudvp(:) dcdrho(:) dudrho(:)];
%
% jbrussell 11/20/2022: Include density kernels in addition to Vp and Vs.
% Also, load kernels for all available mode branches.
% 
% jbrussell 6/12/2023: Read anelastic kernels
function [kern, dispersion] = read_anelastic_kernel96(filename,wavetype,varargin)

if nargin==2
    nmode = 0; % default, fundamental model
elseif nargin==3
	nmode = varargin{1};
end

fid = fopen(filename,'r');
if(fid>=0)	
	it = 0;
	while 1
    	stemp=fgetl(fid);
	    if (strmatch('Anelastic',stemp))    % find where Elastic Love wave kernel begins
%             Anelastic Rayleigh wave: Period=    47.000 Mode =    0 C=     3.931 U=     3.960 GAMMA=  1.578E-04
        	C = textscan(stemp,'%s%s wave: Period=%f Mode =%d C=%f U=%f GAMMA=%f' );
            elastic_anelastic = C{1};
            rayl_love = C{2};
            period_k = C{3};
            mode_k = C{4};
            phv_k = C{5};
            grv_k = C{6};
            gamma_k = C{7};
            if mode_k == nmode
                dispersion.period = period_k;
                dispersion.mode = mode_k;
                dispersion.phv = phv_k;
                dispersion.grv = grv_k;
                dispersion.gamma = gamma_k;
                break
            else % wrong mode, keep going...
                continue
            end
        end
        it = it + 1;
        if it > 1e6
            kern = [];
            data = [];
            return
        end
	end
	stemp = fgetl(fid);	% skip a line	

 % start reading	
	il=1;
	while 1 % when line hasn't hit the dash line
		stemp = fgetl(fid);		
	    temp = sscanf(stemp,'%f');

		if ~isempty(temp) 
			thick(il) = temp(2);% thickness H(i)
			if(wavetype=='L')
				kern.dcdvs(il) = temp(3);% 
				kern.dudvs(il)  = temp(4);% 
	
				kern.dcdvp(il)  = 0;% 
				kern.dudvp(il)  = 0;% 
                
                kern.dcdrho(il)  = temp(5);% 
				kern.dudrho(il)  = temp(6);% 
				
			elseif(wavetype=='R')
                % dc/da      dc/db      dU/da      dU/db      dc/dh      dU/dh     dc/dQai    dc/dQbi    dU/dQai    dU/dQbi    dg/dQai    dg/dQbi
				kern.dcdvp(il)  = temp(3);% 
				kern.dcdvs(il)  = temp(4);% 			
	
				kern.dudvp(il) = temp(5);% 
				kern.dudvs(il) = temp(6);% 
                
                kern.dcdrho(il)  = temp(7);% 
				kern.dudrho(il)  = temp(8);% 
                
                kern.dcdQp(il) = temp(9);%
                kern.dcdQs(il) = temp(10);%
                
                kern.dudQp(il) = temp(11);
                kern.dudQs(il) = temp(12);
                
                kern.dgdQp(il) = temp(13);
                kern.dgdQs(il) = temp(14);
			end
			
			il=il+1;
		else
			break
		end
	end
% 	kern=[dcdvs(:) dcdvp(:) dudvs(:) dudvp(:) dcdrho(:) dudrho(:) dcdQs(:) dcdQp(:) dudQs(:) dudQp(:) dgdQs(:) dgdQp(:)];
	fclose(fid);
else
	error([filename,' does not exist']);% return empty data if file is empty 
	kern=[];
return
end

