function [G_xy1_s, G_xy1_fh] = Shearmodxy_components(h,l,ks,t,theta, khf)

% stretch

noma=ks.*cos(theta);
nomb=h./l+sin(theta);
nom=noma.*nomb;

a= (cos(theta)).^2;
b=(h./l+sin(theta)).*sin(theta);

denom=t.*((a+b).^2);

G_xy1_s=nom./denom;

% flex hinge

nom2=khf.*((h./l)+sin(theta));

c= t.*((h./l).^2);
d=(1+(2.*h./l)).*cos(theta);

denom2=c.*d;

G_xy1_fh=nom2./denom2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by O. Duncan,  Department of Engineering,   %
% Manchester Metropolitan University                                       %
%  M15 6BG, UK.                                                            %
% Please sent your comments to: o.duncan@mmu.ac.uk                         %
%                                                                          %
% The code is described in the paper:                                                     %
% Duncan O. Fast optimisation of honeycombs for impact protection                                %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, I shall not be liable in any event        %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%