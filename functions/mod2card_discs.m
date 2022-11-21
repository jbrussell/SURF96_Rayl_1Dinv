function [card_out] = mod2card_discs(model,discs,cardn)
% Convert surf96 model file to a MINEOS card file with discontinuities
% properly inserted

% Add discontinuities to model
[vs, z] = mod2knot_discs(model(:,3),model(:,1),discs);
[vp, z] = mod2knot_discs(model(:,2),model(:,1),discs);
[rho, z] = mod2knot_discs(model(:,4),model(:,1),discs);

card_out = cardn;
rad = (6371-z)*1000;
Ireplace = find(card_out.z < max(z));
card_out.rho(Ireplace)=[]; card_out.rho =[card_out.rho; flip(rho)*1000];
card_out.vpv(Ireplace)=[]; card_out.vpv =[card_out.vpv; flip(vp)*1000];
card_out.vph(Ireplace)=[]; card_out.vph =[card_out.vph; flip(vp)*1000];
card_out.vsv(Ireplace)=[]; card_out.vsv =[card_out.vsv; flip(vs)*1000];
card_out.vsh(Ireplace)=[]; card_out.vsh =[card_out.vsh; flip(vs)*1000];
%     z_card = card_out.z(Ireplace) + flip([0:length(Ireplace)-1]*1e-10);
z_card = card_out.z + flip([0:length(card_out.z)-1]'*1e-10);
eta_card = interp1(z_card,card_out.eta,flip(z));
qmu_card = interp1(z_card,card_out.qmu,flip(z));
qkap_card = interp1(z_card,card_out.qkap,flip(z));
card_out.eta(Ireplace)=[]; card_out.eta =[card_out.eta; eta_card];
card_out.qmu(Ireplace)=[]; card_out.qmu =[card_out.qmu; qmu_card];
card_out.qkap(Ireplace)=[]; card_out.qkap =[card_out.qkap; qkap_card];

card_out.z(Ireplace)=[]; card_out.z =[card_out.z; flip(z)];
card_out.rad(Ireplace)=[]; card_out.rad =[card_out.rad; flip(rad)];


end

