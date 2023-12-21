function [ strainx, strainy,r] = Strain( theta, delta,l,h,i,deltah)

% lat
ll=length(theta(:,1));


a(:,1)=l(:,1).*cos(theta(:,1))-(delta(:,1).*sin(theta(:,1)));
b=l(i,1).*cos(theta(i,1))-(delta(i,1).*sin(theta(i,1)));

stretchx=a(:,1)./b;

if stretchx(1,1)<=0
x(1)=find(stretchx(1:i,1)<=0,1,'last');
else x(1)=1;
end
if stretchx(ll,1)<=0
x(2)=find(stretchx(i:end,1)<=0,1,'first');
x(2)=x(2)+i;
else x(2)=ll;
end

% ax
c(:,1)=l(:,1).*sin(theta(:,1))+(delta(:,1).*cos(theta(:,1))+h)+deltah(:,1);
d=l(i,1).*sin(theta(i,1))+(delta(i,1).*cos(theta(i,1))+h(i,1));

stretchy=c(:,1)./d;
if stretchy(1,1)<=0
y(1)=find(stretchy(1:i,1)<=0,1,'last');
else y(1)=1;
end
if stretchy(ll,1)<=0
y(2)=find(stretchy(i:end,1)<=0,1,'first');
y(2)=y(2)+i;
else y(2)=ll;
end


r(1)=max([x(1),y(1)]);
r(2)=min([x(2),y(2)]);

strainx(r(1):r(2),1)=real(log(stretchx(r(1):r(2),1)));
strainx(1:r(1),1)=NaN;
strainx(r(2):ll,1)=NaN;
strainy(r(1):r(2),1)=real(log(stretchy(r(1):r(2),1)));
strainx(1:r(1),1)=NaN;
strainy(r(2):ll,1)=NaN;
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
