% reproduce 1D model in "2017-Coherent terahertz spin-wave emission assoc
%iated with ferrimagnetic domain wall dynamics-PRBr-Se-Hyeok Oh"
clear all;clc;close all;tic
%control parameter
conf_file();
DMIenable=0;
pcs=2;%1,optilex 7040 2,acer laptop 3,Landauer 4,gold-c01
rk4=1;%1:rk4,0:heun Method,2:4th predictor-corrector
gpusave=2e-12;%how often saving gpu data
debugg=0;
loadstartm=0;%1:load mat file; 0:direct calculate
systemselec=2;%1:50%-50% TM-RE,2:FM,3:random FiM
%gpuDevice(1)
%constant
if pcs==1
    addpath('E:\dropbox\Dropbox\phd\code\general\gitcontrol\constant');
elseif pcs==2
    addpath('E:\Dropbox\Dropbox\phd\code\general\gitcontrol\constant');
    addpath(genpath('E:\Dropbox\Dropbox\phd\code\general\nogitcontrol'));
elseif pcs==3
    addpath('/home/a0132576/code/general/general/constant');
else
    addpath('/home/svu/a0132576/code/general/constant');
end
constantfile;
clear gam
%params from paper
Asim=7.5e-3*ele;%[J], exchange
Ksim=0.4e-3*ele;%[J], easy-axis anisotropy
kksim=0.2e-6*ele;%[J], DW hard-axis anisotropy
if DMIenable
    Dsim=250e-6*ele;%[J], DMI
else
    Dsim=0;
end
d=0.4e-9;%[m],distance between two spin
alp=0.002;%Gilbert damping
gTM=2.2;gRE=2;%g-factor
T_=linspace(0,400,1);
mTMy_=(-0.6855*T_+1242)*1e3;%[A/m]check "Dropbox\phd\code\papers_books\2017-
%Coherent terahertz spin-wave emission associated with ferrimagnetic
%domain wall dynamics-PRBr-Se-Hyeok Oh\fig1b\fig1b.m"
mREy_=(-0.462*T_+1105)*1e3;%[A/m]
delta_sy_=0.009204*T_-1.8;
%% electrical parameters
jc=0;%[A/m2]
Hext=[0,0,0e-3];%corresponding to runtime
%% FiM parameters
tz=d;%[m],thickness of FiM
TA_paper=201.66;%from Fig. 1(b) in paper
TA=149.4642;%calculated from muTM,muRE,gamTM,gamRE in code
T_=[108.6484,130.3781,152.1078,173.8375,217.2968,239.0265,260.7562,282.4859];
%corresponds to delta_s=[-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8]e-7 [J.s/m^3]
T=108;%[K]
msTM=(-0.6855*T+1242)*(1e3);%[A/m]
msRE=(-0.462*T+1105)*(1e3);
muTM=msTM*d^3;%[A.m^2=J/T]
muRE=msRE*d^3;%[A.m^2=J/T]
gamTM=gTM*mub/(hbar*ele);%1/(s.T)refer to "PRL 97, 217202 (2006), Jiang Xin"
gamRE=gRE*mub/(hbar*ele);
%% SOT parameters
SOT_DLT=1;%1(0),enable(disable) SOT damping torque
SOT_FLT=0;%1(0),enable(disable) SOT field-like torque
psjSHE=[0,1,0];%spin flux polarization
psjSHEx=psjSHE(1);
psjSHEy=psjSHE(2);
psjSHEz=psjSHE(3);
thetaSH=0.2;%spin hall angle
if SOT_FLT
    chi=0;%ratio of FLT/DLT
else
    chi=0;
end
%% time control
tstep=2e-15;
runtime=1*gpusave;%second run for dw motion
savetstep=1;%to reduce data size
gpusteps=round(gpusave/tstep);
totstep=round(runtime/tstep);
%% other parameters
%natom=1000;
natom=10;
delta_sy=0.009204*T-1.8;%[e-7J.s/m^3]
bc=1;%0.periodic condition;1,not periodic
loc_=linspace(0,(natom-1)*d,natom);%atom location
% if (SOT_DLT || SOT_FLT || STT_DLT || STT_FLT) && ~(rk4==1)
%     error('only rk4 is implemented for spin torque driven')
% end
systemgeneration();
if loadstartm
    clear m_
    load('startm_natom1000.mat')
    m_=[mmxstart;mmystart;mmzstart];
end
%dynamic calc
rk4_4llg(); toc
if pcs==3 || pcs==4
    save('finalgpu_tt8.mat')
else
    if (1)%compare with previous result
        figure
        plot(tt,mmx(:,1)','-b','LineWidth',2);
        xlabel('time(ns)','fontsize',15);ylabel('mx','fontsize',15)
        set(gca,'fontsize',20)
        %xlim([0,15]);ylim([-1.05,1.05]);
%         hold on
%         clear all
%         load('test.mat');
%         plot(tt,mmx(:,1)','-r','LineWidth',1);
%         legend('new','old')
    end
end





