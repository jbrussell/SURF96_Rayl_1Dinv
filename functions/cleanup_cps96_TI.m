% Remove temporary files written by the CPS transversely isotropic programs
% (tprep96, tdisp96, tregn96, tlegn96, tdpegn96, tdpder96) and by the *_TI
% wrapper functions. Missing files are silently ignored.
%
% jbrussell 11/2026
%
function cleanup_cps96_TI()

files = {'periods.txt','start.mod','tdisp96.dat','tdisp96.ray','tdisp96.lov', ...
    'tregn96.der','tlegn96.der','tregn96.egn','tlegn96.egn', ...
    'TREGN.ASC','TLEGN.ASC','TRDER.TXT','TLDER.TXT','TRDER.PLT','TLDER.PLT', ...
    'TREGNC.PLT','TREGNU.PLT','TREGNG.PLT','TREGNA0.PLT','TREGNE.PLT', ...
    'TLEGNC.PLT','TLEGNU.PLT','TLEGNG.PLT','TLEGNA0.PLT','TLEGNE.PLT'};
for ii = 1:length(files)
    if exist(files{ii},'file')
        delete(files{ii});
    end
end

end
