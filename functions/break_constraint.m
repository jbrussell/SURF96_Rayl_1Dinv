function [X00] = break_constraint(X00, z, z_brks)
% Break constraint equation at specified depths
%
% z: depth to the bottom of each layer. A break at z_brks lies between
% layer I_brk-1 (whose bottom is closest to z_brks) and layer I_brk.
%
% First derivative (row i: m(i) - m(i-1)): only row I_brk crosses the
% break.
% Second derivative (row i: -m(i-1) + 2m(i) - m(i+1)): rows I_brk-1 and
% I_brk cross the break. Row I_brk+1 (layers I_brk..I_brk+2) lies entirely
% below it and must be kept. (Zeroing it, as was done before 11/2026, left
% layer I_brk in no second-derivative equation, so the layer just below a
% break -- e.g. the uppermost mantle -- was nearly unconstrained and could
% take spurious, large velocity changes.)
%
% jbrussell 6/5/2020
    for ibrk = 1:length(z_brks)
        if z_brks(ibrk) > max(z)
            continue
        end
        [~,I_brk] = min(abs(z-z_brks(ibrk)));
        I_brk = I_brk + 1;
        if I_brk > size(X00,1)
            continue
        end
        if length(find(X00(I_brk,:)~=0))==3 || length(find(X00(I_brk-1,:)~=0))==3 % second derivative
            X00(I_brk-1:I_brk,:) = 0;
            continue
        end
        if length(find(X00(I_brk,:)~=0))==2 % first derivative
            X00(I_brk,:) = 0;
        end
    end
    
end

