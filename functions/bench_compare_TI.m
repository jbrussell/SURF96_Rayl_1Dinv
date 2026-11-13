% Compare a TI result to the corresponding isotropic result and print a
% one-line PASS/FAIL summary (helper for the benchmark_TI_*.m scripts)
%
% [pass, err] = bench_compare_TI(label, a_TI, b_iso, tol, type, logfile)
%
% type: 'rel'  - max over elements of |a-b|/|b|          (dispersion etc.)
%       'peak' - max over columns of max|a-b| / max|b|   (kernels, eigenfunctions)
% Elements that are NaN in either input are ignored (the number of NaN
% mismatches is reported).
% logfile (OPTIONAL): if given, the summary line is appended to this file
%
% jbrussell 11/2026
%
function [pass, err] = bench_compare_TI(label, a, b, tol, type, varargin)

if nargin < 5 || isempty(type)
    type = 'rel';
end
a = double(a); b = double(b);
if isvector(a); a = a(:); end
if isvector(b); b = b(:); end
if ~isequal(size(a),size(b))
    error('bench_compare_TI: %s: size mismatch',label);
end
nanmis = sum(isnan(a(:)) ~= isnan(b(:)));
mask = ~isnan(a) & ~isnan(b);

switch type
    case 'rel'
        d = abs(a(mask)-b(mask)) ./ max(abs(b(mask)), eps);
        err = max([d; 0]);
    case 'peak'
        err = 0;
        for ic = 1:size(a,2)
            m = mask(:,ic);
            if ~any(m); continue; end
            nrm = max(abs(b(m,ic)));
            if nrm == 0
                nrm = max(max(abs(a(m,ic))), eps);
            end
            err = max(err, max(abs(a(m,ic)-b(m,ic)))/nrm);
        end
end
if ~any(mask(:))
    err = NaN;
end
pass = ~isnan(err) && err <= tol;
if pass; str = 'PASS'; else; str = 'FAIL'; end
line = sprintf('[%s] %-55s err = %9.3e  (tol %8.1e)', str, label, err, tol);
if nanmis > 0
    line = [line, sprintf('  [%d NaN mismatches ignored]', nanmis)];
end
fprintf('%s\n', line);
if nargin >= 6 && ~isempty(varargin{1})
    fid = fopen(varargin{1},'a');
    fprintf(fid,'%s\n',line);
    fclose(fid);
end

end
