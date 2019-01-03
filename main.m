% reproduce 1D model in "2017-Coherent terahertz spin-wave emission assoc
%iated with ferrimagnetic domain wall dynamics-PRBr-Se-Hyeok Oh"
clear all;clc;close all;tic
%control parameter
conf_file();
DMIenable=1;
rk4=1;%1:rk4,0:heun Method,2:4th predictor-corrector
gpusave=20e-12;%how often saving gpu data
debugg=0;
loadstartm=1;%1:load mat file; 0:direct calculate
systemselec=1;%1:50%-50% FiM,2:FM,3:random FiM 4:AFM
gpuDevice(1)
constantfile;
clear gam
switch systemselec
    case 1
        A_TMRE=15e-3*ele;%[J], exchange, prefer AFM
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
Ksim=0.08e-3*ele;%[J], easy-axis anisotropy
kksim=0.08e-6*ele;%[J], DW hard-axis anisotropy
if DMIenable
    Dsim=0.05e-3*ele;%[J], DMI
else
    Dsim=0;
end
d=0.4e-9;%[m],distance between two spin
alp=0.02;%Gilbert damping
gTM=2.2;gRE=2;%g-factor
%% electrical parameters
jc=256e9;%[A/m2]
Hext=[0,0,0e-3];%corresponding to runtime
%% FiM parameters
tz=d;%[m],thickness of FiM
T0=300;%[K] initial temperature
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
tstep=2000e-18;
runtime=200*gpusave;%second run for dw motion
gpusteps=round(gpusave/tstep);
ct3max=round(runtime/gpusave);
numpersave=5;
savetstep=gpusteps/numpersave;%to reduce data size
totstep=numpersave*ct3max;
%% other parameters
natom=12500;
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
    load('startm_natom12500_x234.mat')
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





