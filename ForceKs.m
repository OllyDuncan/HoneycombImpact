function [ksn,kh,kf,khf,khb,kfb,khfb] = ForceKs(Es,w,t,l,q,typ)

ks=Es.*w.*t./l;
kf=Es.*w.*(t.^3)./(((l.^3)).*typ); %after buckling, there are effectively three buckling oblique ribs
kh=Es.*w.*(t.^3)./(((l.^2).*q).*typ); %

ksn=ks./ks;
kf=kf./(ks.*typ);
kh=kh./(ks.*typ);
khf=(t.^2)./(l.^2 + 6.*l.*q.*typ);

khb=Es.*w.*(t.^3)./(6.*(l.^2).*q.*typ); %after hinging is set along whole rib
khb=khb./(ks.*typ);
kfb=khb;
khfb=(t.^2)./(7.*(q.^2).*typ);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by O. Duncan,  Department of Engineering,   %
% Manchester Metropolitan University                                       %
%  M15 6BG, UK.                                                            %
% Please sent your comments to: o.duncan@mmu.ac.uk                         %
%                                                                          %
% The code is described in the paper:                                               %
% Duncan O. Fast optimisation of honeycombs for impact protection                                %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, I shall not be liable in any event        %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
