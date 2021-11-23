% read surface wave sensitivity kernels wrote by srfker96 command
%	kernel=[dcdvs(:) dcdvp(:) dudvs(:) dudvp(:)];
% 
function kernel= read_kernel96(filename,wavetype)


fid = fopen(filename,'r');
if(fid>=0)	
	it = 0;
	while 1
    	stemp=fgetl(fid);
	    if (strmatch('Elastic',stemp));    % find where Elastic Love wave kernel begins
        	break
        end
        it = it + 1;
        if it > 100
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
				
			elseif(wavetype=='R')
				dcdvp(il)  = temp(3);% 
				dcdvs(il)  = temp(4);% 			
	
				dudvp(il) = temp(5);% 
				dudvs(il) = temp(6);% 
			end
			
			il=il+1;
		else
			break
		end
	end
	kernel=[dcdvs(:) dcdvp(:) dudvs(:) dudvp(:)];
	fclose(fid);
else
	error([filename,' does not exist']);% return empty data if file is empty 
	kernel=[];
return
end

