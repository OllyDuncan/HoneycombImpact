function [E_offaxis] = Offset_Linear_Mod(axis,E_parallel,G,NU_parallel,E_perpendicular)

a= ((cos(axis)).^4)./E_parallel;
b= (cos(axis).^2).*(sin(axis).^2);
c =(1./G)-(2.*NU_parallel./E_parallel);
d= ((sin(axis)).^4)./E_perpendicular;

E_offaxis(:,1)=1./(a+b.*c+d);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by O. Duncan,  Department of Engineering,   %
% Manchester Metropolitan University                                       %
%  M15 6BG, UK.                                                            %
% Please sent your comments to: o.duncan@mmu.ac.uk                         %
%                                                                          %
% The code is described in the paper:                                                    %
% Duncan O. Fast optimisation of honeycombs for impact protection                                %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, I shall not be liable in any event        %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%