function [Ey_offaxis] = Offset_Axialmod_def(axis,Ey,Gxy,NU_yx,Ex)

%axis=offset angle

a= ((cos(axis)).^4)./Ey;
b= (cos(axis).^2).*(sin(axis).^2);
c =(1./Gxy)-(2.*NU_yx./Ey);
d= ((sin(axis)).^4)./Ex;

Ey_offaxis(:,1)=1./(a+b.*c+d);
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