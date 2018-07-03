% reproduce 1D model in "2017-Coherent terahertz spin-wave emission assoc
%iated with ferrimagnetic domain wall dynamics-PRBr-Se-Hyeok Oh"
clear all;clc;close all;tic
%control parameter
conf_file();
DMIenable=0;
rk4=1;%1:rk4,0:heun Method,2:4th predictor-corrector
gpusave=5e-12;%how often saving gpu data
debugg=1;
loadstartm=0;%1:load mat file; 0:direct calculate
systemselec=1;%1:50%-50% FiM,2:FM,3:random FiM 4:AFM
%gpuDevice(1)
constantfile;
clear gam
switch systemselec
    case 1
        A_TMRE=7.5e-3*ele;%[J], exchange, prefer AFM
        A_TMTM=A_TMRE;
        A_RERE=A_TMRE;
    case 2 
        A_TMTM=-1.5e-21;%[J], exchange,prefer FM
        A_RERE=-0.98e-21;
    case 3
        compositionn=0.3;%[J]
        A_TMTM=-1.5e-21;%[J], same as PRB 97, 184410 (2018)
        A_RERE=-0.98e-21;%[J]
        A_TMRE=7.63e-21;%[J]
    case 4
        A_TMTM=1.5e-21;%[J], exchange,prefer AFM
        A_RERE=0.98e-21;
end
%params from paper
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
runtime=3*gpusave;%second run for dw motion
savetstep=10;%to reduce data size
gpusteps=round(gpusave/tstep);
totstep=round(runtime/tstep);
%% other parameters
%natom=1000;
natom=20;
delta_sy=0.009204*T-1.8;%[e-7J.s/m^3]
bc=1;%0.periodic condition;1,not periodic
loc_=linspace(0,(natom-1)*d,natom);%atom location
systemgeneration();
if loadstartm
    clear m_
    if systemselec==3
        if compositionn==0.3
        load('startm_TT108_natom200_pc4_0.3.mat')
        elseif compositionn==0.4
            load('0.5.mat')
        elseif compositionn==0.5
            load('startm_TT108_natom200_pc4_0.5.mat')
        elseif compositionn==0.6
            load('startm_TT108_natom200_pc4_0.6.mat')
        else
           error('this composition is not inilitized') 
        end
    end
    m_=[mmxstart;mmystart;mmzstart];
end
systemgenerationB();
if (0)%view initial state
   dwplotstep=1;
    figure%initial magnetization
    hold on
    for ct=1:dwplotstep:natom
        if mark_(ct)==1
            quiver3(0,loc_(ct)*1e9,0,m_(1,ct),m_(2,ct),m_(3,ct),'r');
        else
            quiver3(0,loc_(ct)*1e9,0,m_(1,ct),m_(2,ct),m_(3,ct),'b');
        end
    end
    xlim([-1 1]);ylim([-5 d*1e9*natom]);zlim([-2 2]);
    xlabel('x axis');ylabel('y axis');zlabel('z axis');
    view(-27,20) 
end
%dynamic calc
rk4_4llg(); toc
if ~debugg
    save('final.mat')
else
    if (1)%compare with previous result
        %save('test.mat');
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





