% Function to read in the initial model card and plot it to double check
% that everything looks ok
%
% varargin used to take in name of model card for iteration step.
% if no name is given then the model is named init_model
% NJA, 2014

function card=read_model_card(CARD)

% datapath = 'data/';


isfigure = 0;
% Parameters spcefic to the format of these files
hlines1 = 3;


fid = fopen(CARD,'r');
A=textscan...
    (fid,'%f %f %f %f %f %f %f %f %f','headerlines',hlines1);

card.fname=CARD;
card.z = 6371-A{1}/1000;
card.rad = A{1};
card.rho = A{2};
card.vpv = A{3};
card.vph = A{7};
card.vsv = A{4};
card.vsh = A{8};
card.eta = A{9};
card.qmu = A{6};
card.qkap = A{5};

fclose(fid);

end