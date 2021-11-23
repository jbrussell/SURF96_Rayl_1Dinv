% write_MINE_mod
% Takes in an input file and writes out finely sampled model card for
% MIENOS
% We only care about the top 250 km because the rest of the model will be
% identical to the anisotropic PREM model we already have
%
% NJA, 2014]
%
% BEWARE -- modifications were made to simplify things for a linear test
%
% JOSH 4/2/17 - Simplified to just write out a card file
%
% COLUMNS OF CARD FILE: 1   2  3   4     5     6     7   8   9
%                       R,RHO,VPV,VSV,QKAPPA,QSHEAR,VPH,VSH,ETA
%
function [] = write_MINEOS_mod(ncard,CARD)

% CARD : Name of output card file
% ncard : input card structure to write out
%        ncard.rad 
%        ncard.rho 
%        ncard.vpv 
%        ncard.vsv 
%        ncard.qkap 
%        ncard.qmu 
%        ncard.vph 
%        ncard.vsh 
%        ncard.eta

warning('off','all');

card = ncard;

%% Now begin to write things out to the new model card

fid=fopen(CARD,'w');

% First the header information (3 lines)

% Line 1 - Model Card Name
fprintf(fid,'%s\n',CARD);

% Line 2 - ifanis, tref, ifdeck
ifanis=1;
trec = -1;
ifdeck = 1;

fprintf(fid, '%i\t%f\t%i\n',[ifanis trec ifdeck]);

% Line 3 - N, nic, noc
N = length(card.rho);
ind_oc = find(card.vsv==0 & card.vpv>7000); % index liquid outer core
nic = ind_oc(1)-1; % top of inner core
noc = ind_oc(end); % top of outer core

fprintf(fid,'%3i\t%2i\t%2i\n',[N nic noc]);

%% Loop through layers
% Loop until end - r, rho, vpv, vsv, qkapp, qshear, vph, vsh, eta
% Remember that r means radius (not depth!)
% Everything else is in m/s (not km/s!)

count = 0;

% Loop through each model layer and write out
for id = 1:length(card.rho)
    del = ' ';
    fprintf(fid,'%8.0f%9.2f%9.2f%9.2f%9.1f%9.1f%9.2f%9.2f%9.5f\n'...
        ,[card.rad(id) card.rho(id) card.vpv(id) card.vsv(id) card.qkap(id) card.qmu(id) card.vph(id) card.vsh(id) card.eta(id)]);
    count = count+1;
end


fclose(fid);

%% Check that number of layers match the count
if count ~= N
    disp(['N : ',num2str(N),' COUNT : ',num2str(count)]);
    error('Mismatch in layer count!');
end

%turn the warnings back on b/c they can be useful in some cases
warning('on','all')
