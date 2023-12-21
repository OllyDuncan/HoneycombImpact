function [l] = Length(kh,ks,theta,theta_1,l_1)

a=-2*kh/ks;
b=log(sin(theta(:,1))/sin(theta_1));

c=l_1^2;

l(:,1)=((a*b)+c).^0.5;
l(:,1)=real(l(:,1));
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

