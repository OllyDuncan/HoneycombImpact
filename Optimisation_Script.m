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
% clear all; close all; clc
%% Problem formulation
EIn=12.5; %joules
SHeight=48;
SWidth=40;
SDepth=40;
Ess = 518*10^6; % Pa, of polyester
%% Initial Geometries
h0=1; %relative oblique rib length
t=0.1; %rib width
theta0=-20; %initial rib angle
Offset=0; %offset angle (degrees) - set between 0 & 90
%% change cell wall angle - find best angle to nearest 5 degrees
i=1;
c=0.9^i;
[out{1,1},Buckle(1),Densif(1),Completed(i),Absorbed(i)] = Impact_function(EIn,SHeight,SWidth,SDepth,Ess,[c,c,c],h0,t,theta0(1),Offset,1,3);
Plat(1)=(SHeight-(out{1,1}(Densif(1),1) - out{1,1}(Buckle(1),1)))./SHeight;
df=5;
int=5;
 for i = 2:14 
 c=0.9^i;
 theta0(i)=theta0(i-1)+int;
 [out{i,1},Buckle(i,1),Densif(i,1),Completed(i,1),Absorbed(i,1)] = Impact_function(EIn,SHeight,SWidth,SDepth,Ess,[c,c,c],h0,t,theta0(i),Offset,1,3);
 Plat(i,1)=out{i,1}(Densif(i,1),1) - out{i,1}(Buckle(i,1),1);
 end
[Plat1,ii]=max(Plat(:,1));
%% change h length - find best value to nearest 0.05
inth=0.05;
h02(1)=0.75;
 for i = 2:16
 c=0.9^i;
 h02(i)=h02(i-1)+inth;
 [out{i,2},Buckle(i,2),Densif(i,2),Completed(i,2),Absorbed(i,2)] = Impact_function(EIn,SHeight,SWidth,SDepth,Ess,[c,c,c],h02(i),t,theta0(ii),Offset,2,3);
 Plat(i,2)=out{i,2}(Densif(i,2),1) - out{i,2}(Buckle(i,2),1);
 end
[Plat2,iii]=max(Plat(:,2));
%% tune wall thickness, so that the energy is absorbed at the point of densification
Completed(1,3)=Completed(iii,2);
Densif(1,3)=Densif(iii,2);
i=1;
if Completed(1,3) > Densif(i,3)
    while Completed(i,3) > Densif(i,3)
    c=0.9^i;
    [out{i,3},Buckle(i,3),Densif(i+1,3),Completed(i+1,3),Absorbed(i,3)] = Impact_function(EIn,SHeight,SWidth,SDepth,Ess,[c,c,c],h02(iii),t(i),theta0(ii),Offset,3,3);
    i=i+1;
     t(i) = t(i-1) + 0.01;
    end
elseif Completed(1,3) < Densif(i,3)    
    while Completed(i,3) > Densif(i,3)
    c=0.9^i;
    [out{i,3},Buckle(i,3),Densif(i+1,3),Completed(i+1,3),Absorbed(i,3)] = Impact_function(EIn,SHeight,SWidth,SDepth,Ess,[c,c,c],h02(iii),t(i),theta0(ii),Offset,3,3);
    i=i+1;
     t(i) = t(i-1) - 0.01;
    end
    if i>1
    i=i-1;
    else
    end
else
end
%%

t(i)
h02(iii)
theta0(ii)