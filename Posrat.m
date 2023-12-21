function [NU_axial, NU_trans] = Posrat(khf,ks,h,l,theta)

% Axial - load parallel to point of hexagon

noma(:,1)=sin(theta(:,1));
nomb(:,1)=(cos(theta(:,1)).^2);
nomc=1/khf-1/ks;
% 
nom(:,1)=noma(:,1).*nomb(:,1).*nomc;

% a= w*(cos(theta1));
a(:,1)=(sin(theta(:,1))+h./l(:,1));
b(:,1)=(sin(theta(:,1)).^2).*(1/khf);
c(:,1)=(cos(theta).^2)/ks;

denom(:,1)=a(:,1).*(b(:,1)+c(:,1));

NU_axial(:,1)=nom(:,1)./denom(:,1);

% Transvrse

nomd(:,1)=sin(theta(:,1));
nome(:,1)=h./l(:,1)+sin(theta(:,1));
nomf=1/khf-1/ks;
% 
nom2(:,1)=nomd(:,1).*nome(:,1).*nomf;

% a= w*(cos(theta1));
c(:,1)=(cos(theta(:,1)).^2).*(1/khf);
d(:,1)=((sin(theta(:,1)).^2)+2*h./l(:,1))/ks;

denom2(:,1)=c(:,1)+d(:,1);


NU_trans(:,1)=nom2./denom2;

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