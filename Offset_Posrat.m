function [Nu_offaxis] = Offset_Posrat(axis,E_parallel,G,NU_parallel,E_perpendicular,Eperp_offaxis)

a= (((cos(axis)).^4)+((sin(axis)).^4)).*(NU_parallel./E_parallel);
b= (cos(axis).^2).*(sin(axis).^2).*((1./E_perpendicular)+(1./E_parallel)-(1./G));


Nu_offaxis(:,1)=Eperp_offaxis.*(a-b);


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