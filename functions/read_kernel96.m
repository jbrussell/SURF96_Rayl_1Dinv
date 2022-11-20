% read surface wave sensitivity kernels wrote by srfker96 command
%	kernel=[dcdvs(:) dcdvp(:) dudvs(:) dudvp(:) dcdrho(:) dudrho(:)];
%
% jbrussell 11/20/2022: Include density kernels in addition to Vp and Vs.
% Also, load kernels for all available mode branches.
% 
function kernel= read_kernel96(filename,wavetype,varargin)

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
	    if (strmatch('Elastic',stemp))    % find where Elastic Love wave kernel begins
        	C = textscan(stemp,'%s%s wave: Period=%f Mode =%d C=%f U=%f' );
            elastic_anelastic = C{1};
            rayl_love = C{2};
            period_k = C{3};
            mode_k = C{4};
            phv_k = C{5};
            grv_k = C{6};
            if mode_k == nmode
                break
            else % wrong mode, keep going...
                continue
            end
        end
        it = it + 1;
        if it > 1e6
            kernel = [];
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
				dcdvs(il) = temp(3);% 
				dudvs(il)  = temp(4);% 
	
				dcdvp(il)  = 0;% 
				dudvp(il)  = 0;% 
                
                dcdrho(il)  = temp(5);% 
				dudrho(il)  = temp(6);% 
				
			elseif(wavetype=='R')
				dcdvp(il)  = temp(3);% 
				dcdvs(il)  = temp(4);% 			
	
				dudvp(il) = temp(5);% 
				dudvs(il) = temp(6);% 
                
                dcdrho(il)  = temp(7);% 
				dudrho(il)  = temp(8);% 
			end
			
			il=il+1;
		else
			break
		end
	end
	kernel=[dcdvs(:) dcdvp(:) dudvs(:) dudvp(:) dcdrho(:) dudrho(:)];
	fclose(fid);
else
	error([filename,' does not exist']);% return empty data if file is empty 
	kernel=[];
return
end

