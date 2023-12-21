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
%% clear all; close all; clc
EIn=10; %joules
%Sample dimensions
SHeight=48;
SWidth=40;
SDepth=40;
%Material properties
Ess = 518*10^6; % Pa
Ess=((Ess-114).*48./SHeight)+114;
col=[0,0,0];
maxrun=1;
run=1;
%% Geometries
h0=1; %relative oblique rib length
t=0.1; %rib width
theta0=30; %initial rib angle (must be between +/- 70 (not inclusive))
Offset=0; %offset angle (degrees) - set between 0 & 90
%% Geometries
l0=1; %set
w=1; %relative rib width
theta0=(theta0+0.1)*pi/180; %marginal offset to avoid indeterminate responses
inc=0.01;
thetadeg(:,1)= -70:inc:70; %fills a vector with potential rib angles
theta(:,1)=(thetadeg(:,1)+0.1)*pi/180;
delta0=0; %length change
Offset=Offset*pi/180;
Chir=sqrt((6.*Offset./pi-round(6.*Offset./pi)).^2);
q=(0.045/(((h0.^2))*(t)))*(cos((pi*Chir)^3)^2);  %relative length over which bending happend
if q > 1
    q=1;
end
zer=find(theta(:,1)==theta0);
%% Properties
Es=1; %base material stiffness (normalised)
nus=0.45; %base material Poisson's ratio
Gs=Es./(2.*(1+nus)); %rib shear modulus
itt=10; % Newton-rapson method iterations
it=10; %makes line black
%% Force constants
[ks,kh,kf,khf,~,~,~] = ForceKs(Es,w,t,l0,q,1); %select values related to local bending at hinges
%% Models - linear elastic (theta dependent)
% Length 
[l(:,1)] = Length(kh,ks,theta,theta0,l0);
[h(:,1)] = LengthH(kh,ks,theta,theta0,l,h0); 
%delta
[ delta(:,1)] = Delta(kh,kf,theta,theta0,l,delta0);
% Strain
[ strainx, strainy,~] = Strain( theta,delta,l,h,zer,0);
% Poissons ratio
[NU_xy, NU_yx] = Posrat(khf,ks,h,l,theta);
% Elastic Modulus
[Ex, Ey] = Elasticmod_def(theta,t,h,l,khf,ks);
% Shear Modulus
[G_xy1_s, G_xy1_fh] = Shearmodxy_components(h,l,ks,t,theta,khf);
Gxy(:,1) = 1./((1./G_xy1_s)+(1./G_xy1_fh)); 
% Relative Density
As= (h*t + 2*l*t)-2*t^2/sqrt(3); %solid area
Hz=sqrt((l.*cos(theta)).^2); 
Ve=l.*sin(theta);
Ai= 2.*Hz.*h + 2.*Hz.*Ve; %cell area (varies)
Rho=As./Ai;
dens=1-(Rho(zer));
endd=log((1-dens)/10);
%% Offset
[Ex_offaxis] = Offset_Linear_Mod(Offset,Ex,Gxy,NU_xy,Ey);
[Ey_offaxis] = Offset_Linear_Mod(Offset,Ey,Gxy,NU_yx,Ex);
[Nuxy_offaxis] = Offset_Posrat(Offset,Ex,Gxy,NU_xy,Ey,Ex_offaxis);
[Nuyx_offaxis] = Offset_Posrat(Offset,Ey,Gxy,NU_yx,Ex,Ey_offaxis);
[Strainx_Off] = Strain_Offset (strainx, strainy, Offset);
[Strainy_Off] = Strain_Offset (strainy, strainx, Offset);
%% Set up virtual test (no buckling or densification)
[~,e]=min((Strainy_Off(:,1)-endd).^2); 
    dir = (e-zer)/(sqrt((e-zer)^2));    
%% Stress & critical stress
VirStraina(:,1)=Strainy_Off(zer:dir:e,1);
VirStraint(:,1)=Strainx_Off(zer:dir:e,1);
Dd = length(VirStraina(:,1));   
VirSt(Dd,1)=zeros;
VirSt(2:Dd,1)=cumsum(Ey_offaxis(zer+(1*dir):dir:e,1).*(VirStraina(2:Dd,1)-VirStraina(1:Dd-1,1)));     
Str(:,1)=sqrt((VirSt(:,1).*cos(2.*Offset)./(t.^2)).^2);
Str(:,2)=sqrt((VirSt(:,1).*cos(2.*(Offset+theta(zer:dir:e)))./(t.^2)).^2);
Str(:,3)=sqrt((VirSt(:,1).*cos(2.*(Offset-theta(zer:dir:e)))./(t.^2)).^2);
ql(:,1:3)=sqrt(3.*sqrt(3.*Str(:,1:3)));
St(:,1)=ql(:,1).*tan(ql(:,1)./2)-ql(:,2).*cot(ql(:,2)./2)-ql(:,3).*cot(ql(:,3)./2);
St(:,2)=ql(:,2).*tan(ql(:,2)./2)-ql(:,1).*cot(ql(:,1)./2)-ql(:,3).*cot(ql(:,3)./2); %3rd variant not required - rotations only defined in one axis (symmetery used otherwise)
%% Critical Uniaxial Buckling 
Bcc(1:2,1:2)=length(VirSt);
for i=1:2
[Bcc(:,i)]=find(islocalmin(St(:,i).^2)==1,2);
if  St(Bcc(1,i),i)^2 > 1
    Bcc(1,i) = Bcc(2,i) ;   
else
end
end
[Bc,typ]=min(Bcc(1,:));
BcN=zer+(Bc.*dir);
if typ == 1 % set rib for h buckling
    dum=h;
else % set rib for l buckling
    dum=l;
end
%% find angle of buckled rib imediatly after buckling
thi_n(1,1)=20.*pi()./180;
 %Newton-Rapson method
    for i=1:itt
    f(i)=sin(thi_n(i)) -  dum(BcN)./dum(zer).*thi_n(i);
    df(i)=cos(thi_n(i)) - dum(BcN)./dum(zer);
    thi_n(i+1,1)=thi_n(i) - f(i)./df(i);
    con(i)=sqrt((thi_n(i+1)-thi_n(i)).^2);
    end
Bc=Bc+1;
BcN=zer+((1+Bc).*dir);
thi(Bc-1,1)=sqrt((thi_n(itt)).^2);
Off2(1:Bc-1,1)=Offset;
%% Adjust strains
e=e+2*dir;
delta_b(length(strainy),1)=zeros;
    [~,kh2,~,khf2,khb,kfb,khfb] = ForceKs(Es,w,t,dum(zer),dum(zer),3);
    ksb=1./((1./ks)+(1./khfb));

thi(Bc:Dd,1)=zeros;
    thi(Bc:Dd,1)=cumsum((theta(Bc:Dd)-theta(Bc-1:Dd-1)).*khf2./khfb);
    lum(Bc:Dd,1)=(dum(zer)./(2.*thi(Bc:Dd))).*(1-cos(thi(Bc:Dd)));
    lim(Bc:Dd,1)=(dum(zer)./(2.*thi(Bc:Dd))).*sin(thi(Bc:Dd));
    Off2(Bc:Dd,1)= Offset + atan(2*lum(Bc:Dd,1)./lim(Bc:Dd,1));

    [delta_b(BcN:dir:e,1)] = -1.*(Delta(khb,kfb,thi(Bc:Dd),thi(Bc),lum(Bc:Dd,1),delta(zer,1)));
 
     if typ == 1 % set rib for h buckling
     [l(BcN:dir:e,1)] = Length(kh2,ks,theta(BcN:dir:e),theta0,l0);
     [h(BcN:dir:e,1)] = LengthH(kh2,ksb,theta(BcN:dir:e),theta0,l(BcN:dir:e,1),h0);
     else% set rib for l buckling
     [l(BcN:dir:e,1)] = Length(kh2,ksb,theta(BcN:dir:e),theta0,l0);
     [h(BcN:dir:e,1)] = LengthH(kh2,ks,theta(BcN:dir:e),theta0,l(BcN:dir:e,1),h0);
     end

Ju=100; %fill jump created by delta_b with Ju data points, then iterate cell paramemters (& delta)
Lev=Bc+Ju;
LevN=zer+((1+Lev).*dir);
e2=e+Ju*dir;
Dd2=Dd+Ju;

thi(Lev:Dd2)=thi(Bc:Dd); %inserts data points to interpolate accross buckling region
lim(Lev:Dd2)=lim(Bc:Dd);
lum(Lev:Dd2)=lum(Bc:Dd);
Off2(Lev:Dd2)=Off2(Bc:Dd);
delta(LevN:dir:e2,1)=delta(BcN:dir:e);
delta_b(LevN:dir:e2,1)=delta_b(BcN:dir:e);
theta(LevN:dir:e2,1)=theta(BcN:dir:e);
l(LevN:dir:e2,1)=l(BcN:dir:e);
h(LevN:dir:e2,1)=h(BcN:dir:e);

thi(Bc-2:Lev)=thi(Bc-2):((thi(Bc)-thi(Bc-2))/(Ju+2)):thi(Bc); %iterpolate linearly
lim(Bc-2:Lev)=lim(Bc-2):((lim(Bc)-lim(Bc-2))/(Ju+2)):lim(Bc);
lum(Bc-2:Lev)=lum(Bc-2):((lum(Bc)-lum(Bc-2))/(Ju+2)):lum(Bc);
Off2(Bc-2:Lev)=Off2(Bc-2):((Off2(Bc)-Off2(Bc-2))/(Ju+2)):Off2(Bc);
delta(BcN-(2*dir):dir:LevN,1)=delta(BcN-2*dir):((delta(BcN)-delta(BcN-2*dir))/(Ju+2)):delta(BcN);
delta_b(BcN-(2*dir):dir:LevN)=delta_b(BcN):(delta_b(BcN+1*dir)-delta_b(BcN))/(Ju+2):delta_b(BcN+1*dir);
theta(BcN-(2*dir):dir:LevN)=theta(BcN-2*dir):((theta(BcN)-theta(BcN-2*dir))/(Ju+2)):theta(BcN);
l(BcN-(2*dir):dir:LevN)=l(BcN-2*dir):((l(BcN)-l(BcN-2*dir))/(Ju+2)):l(BcN);
h(BcN-(2*dir):dir:LevN)=h(BcN-2*dir):((h(BcN)-h(BcN-2*dir))/(Ju+2)):h(BcN);
if typ==1
    [strainx, strainy,~] = Strain(theta,delta,l,h,zer,delta_b);
else
    [strainx, strainy,~] = Strain(theta,delta,l,h,zer,delta_b);
end  
    [VirStraina(Bc:Dd2,1)] = Strain_Offset (strainy(BcN:dir:e2,1), strainx(BcN:dir:e2,1), Off2(Bc:Dd2,1));
    [VirStraint(Bc:Dd2,1)] = Strain_Offset (strainx(BcN:dir:e2,1), strainy(BcN:dir:e2,1), Off2(Bc:Dd2,1));
%% Adjust tensors
% Poissons ratio
NU_yx(BcN:dir:e2)= -strainx(BcN:dir:e2)./strainy(BcN:dir:e2);
NU_yx(NU_yx > 10) = 10;
NU_yx(NU_yx < -10) = -10;
% Elastic Moduli
[Ex(BcN:dir:e2), Ey(BcN:dir:e2)] = Elasticmod_def(theta(BcN:dir:e2),t,h(BcN:dir:e2),l(BcN:dir:e2),khf2,ksb);
NU_xy(BcN:dir:e2) = NU_yx(BcN:dir:e2).*Ex(BcN:dir:e2)./Ey(BcN:dir:e2);
% Shear Modulus
[G_xy1_s(BcN:dir:e2,1), G_xy1_fh(BcN:dir:e2,1)] = Shearmodxy_components(h(BcN:dir:e2),l(BcN:dir:e2),ksb,t,theta(BcN:dir:e2),khf2);
Gxy(BcN:dir:e2,1) = 1./((1./G_xy1_s(BcN:dir:e2,1))+(1./G_xy1_fh(BcN:dir:e2,1))); 
% Offset
[Ex_offaxis(BcN:dir:e2)] = Offset_Linear_Mod(Off2(Bc:Dd2,1),Ex(BcN:dir:e2),Gxy(BcN:dir:e2),NU_xy(BcN:dir:e2),Ey(BcN:dir:e2));
[Ey_offaxis(BcN:dir:e2)] = Offset_Linear_Mod(Off2(Bc:Dd2,1),Ey(BcN:dir:e2),Gxy(BcN:dir:e2),NU_yx(BcN:dir:e2),Ex(BcN:dir:e2));
[Nuxy_offaxis(BcN:dir:e2)] = Offset_Posrat(Off2(Bc:Dd2,1),Ex(BcN:dir:e2),Gxy(BcN:dir:e2),NU_xy(BcN:dir:e2),Ey(BcN:dir:e2),Ex_offaxis(BcN:dir:e2));
[Nuyx_offaxis(BcN:dir:e2)] = Offset_Posrat(Off2(Bc:Dd2,1),Ey(BcN:dir:e2),Gxy(BcN:dir:e2),NU_yx(BcN:dir:e2),Ex(BcN:dir:e2),Ey_offaxis(BcN:dir:e2));
%% Post Densification 
yden=(h./2)+l.*sin(theta)+(delta+delta_b).*(cos(Offset)) - 2.*t.*sqrt((sin(theta).^2)); 
xden=(l.*cos(theta)-delta.*sin(Offset)-2.*t.*sqrt((cos(theta).^2))); 
Densify = find(yden(zer:dir:e) <= 0,1);
Densifx = find(xden(zer:dir:e) <= 0,1);
 if isempty(Densify)
    Densify = inf;
 end
 if isempty(Densifx)
    Densifx = inf;
 end
Densif=(min(Densify,Densifx));
DensifN=zer+Densif*dir;
    Nuxy_offaxis(DensifN:dir:e2)= nus;
    Nuyx_offaxis(DensifN:dir:e2)= nus;
    % Elastic Modulus
    Ex_offaxis(DensifN:dir:e2)= Es.*Rho(DensifN:dir:e2,1);
    Ey_offaxis(DensifN:dir:e2)= Es.*Rho(DensifN:dir:e2,1);
    % Shear Modulus
    Gxy(DensifN:dir:e2)= Gs.*Rho(DensifN:dir:e2,1);
%% Recalculate Stress
if Bc > Densif
    Bc1=Densif;
elseif Bc < Densif
    Bc1=Bc+1;
end
%failure here suggests t is too large (i.e. densification happens at very low strain)
VirSt(Bc1:Dd2,1)=cumsum(Ey_offaxis(zer+(Bc1*dir):dir:zer+(Dd2*dir)).*(VirStraina(Bc1:Dd2,1)-VirStraina(Bc1-1:Dd2-1,1)));
VirSt(Bc1:Dd2,1) = VirSt(Bc1:Dd2,1)+VirSt(Bc1-1,1);     
VirSt(:,1)=VirSt(:,1).*(SDepth/1000); %in N/m^2, after accounting for w in ks calcaulation
%%
figure(1)
subplot(1,maxrun,run)
[~,Ddi]=min((VirStraina(1:Dd2,1)-endd).^2);
hh(1)=plot(-VirStraina(1:Ddi,1),-VirSt(1:Ddi),'Color',col,'LineWidth',2); hold on;   
xlabel('True Strain','FontSize',24);
ylabel('Stress/E','FontSize',24);
axis('square')
set(gca,'FontSize',24);
set(gca,'FontName','Palatino Linotype');
%% energy - in J/(Es*m^3)
VirStrainEng=exp(VirStraina(1:Densif,1))-1;
VirStrainEngt=exp(VirStraint(1:Densif,1))-1;
VirStT(1:Densif,1) = -VirSt(1:Densif,1).*0;%assumed frictionless contact
VirStraing=(VirStraina(1:Densif,1) - VirStraint(1:Densif,1))./2;
VirStG=(VirSt(1:Densif,1) - VirStT(1:Densif,1))./2;
U1=trapz(-VirStraina(1:Densif,1),-VirSt(1:Densif,1));
U2=trapz(-VirStraint(1:Densif,1),-VirStT(1:Densif,1));
U3=trapz(-VirStraing(1:Densif,1),-VirStG(1:Densif,1));
U=(U1+U2+U3);
MaxStrainT=exp(VirStrainEngt(Densif,1))-1;
MaxStress=(VirSt(Densif,1))./(1-(MaxStrainT)); %engineering stress to obtain applied force
MaxStrain=exp(VirStrainEng(Densif))-1;
%% Force Disp
EngStrainA=exp(VirStraina(1:end,1))-1;
EngStrainT=exp(VirStraint(1:end,1))-1;
Stress=Ess.*VirSt(1:end,1);
Force=-Stress(1:end,1).*SWidth.*SDepth./(10^6);
Disp=-EngStrainA(1:end,1).*SHeight;
NRG(1:Dd2,1)=cumtrapz(Disp(1:Dd2,1),Force(1:Dd2,1));

if max(NRG(:,1))>(EIn*1000)
num=find(NRG(:,1)>(EIn*1000));
disp('energy was absorbed');
else
disp('energy not absorbed - try editing cell parameters, dens (and so endd), or extend theta range');
num=Dd;
end
nn=num(1); 
%%
Disp(nn:end,1)=NaN;
Force(nn:end,1)=NaN;
figure(2)
subplot(1,maxrun,run)
plot(Disp(1:end,1),Force(1:end),'Color',col,'LineWidth',2); hold on;
xlabel('Displacement (mm)','FontSize',24);
ylabel('Force (N)','FontSize',24);
axis('square')
axis([0,48,0,2250]);
set(gca,'FontSize',24);
set(gca,'FontName','Palatino Linotype');
%%
out(:,1)=Disp(:,1);
out(:,2)=Force(:,1);
Buckle=Bc;
Completed = nn;
Absorbed=Dd2-nn;


