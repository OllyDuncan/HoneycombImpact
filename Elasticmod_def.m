function [Ex, Ey] = Elasticmod_def(theta,t,h,l,khf,ks)

%Axial mod

nom(:,1)=cos(theta(:,1));

a= t;
b=h./l(:,1)+sin(theta(:,1));
c =(1/ks.*cos(theta(:,1)).^2)+(1/khf.*sin(theta(:,1)).^2);

Ex(:,1)=nom(:,1)./(a.*b.*c);

% Lateral mod

nom2(:,1) =h./l(:,1)+sin(theta(:,1));

d(:,1)= (cos(theta(:,:)).^2)/(khf);
e(:,1)=((sin(theta(:,1)).^2)+(2*h./l))./(ks);
% f(:,1)=2*h./l./khfb;

denom2(:,1)=t.*cos(theta(:,1)).*(e(:,1)+d(:,1));

Ey(:,1)=nom2(:,1)./denom2(:,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by O. Duncan,  Department of Engineering,   %
% Manchester Metropolitan University                                       %
%  M15 6BG, UK.                                                            %
% Please sent your comments to: o.duncan@mmu.ac.uk                         %
%                                                                          %
% The code is described in the paper:                                                   %
% Duncan O. Fast optimisation of honeycombs for impact protection                                %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, I shall not be liable in any event        %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%