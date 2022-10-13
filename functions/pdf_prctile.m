function [vals] = pdf_prctile(pdf_mat,xvec,prctile)

Nlays = size(pdf_mat,1);
vals = nan(Nlays,1);
for ii = 1:Nlays
    cumrelfreq = cumsum(pdf_mat(ii,:)) ./ sum(pdf_mat(ii,:));
    if isempty(find(cumrelfreq~=1))
        continue
    end
    
    if prctile/100 < min(cumrelfreq)
        vals(ii,1) = min(xvec);
    else    
        vals(ii,1) = interp1(cumrelfreq+(0:length(cumrelfreq)-1)*1e-10,xvec,prctile/100);
    end

end

end
