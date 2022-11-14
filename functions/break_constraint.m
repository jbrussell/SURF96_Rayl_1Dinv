function [X00] = break_constraint(X00, z, z_brks)
% Break constraint equation at specified depths
%
% jbrussell 6/5/2020
    for ibrk = 1:length(z_brks)
        if z_brks(ibrk) > max(z)
            continue
        end
        [~,I_brk] = min(abs(z-z_brks(ibrk)));
        I_brk = I_brk + 1;
        if length(find(X00(I_brk,:)~=0))==3 % second derivative
            if I_brk-1<1 || I_brk+1>length(z)
                continue
            end
            X00(I_brk-1:I_brk+1,:) = 0;
        end
        if length(find(X00(I_brk,:)~=0))==2 % first derivative
            X00(I_brk,:) = 0;
        end
    end
    
end

