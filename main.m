% reproduce 1D model in "2017-Coherent terahertz spin-wave emission assoc
%iated with ferrimagnetic domain wall dynamics-PRBr-Se-Hyeok Oh"
clear all;clc;close all;tic
%control parameter
conf_file();
DMIenable=1;
gpuc=1;%using GPU
pcs=2;%1,optilex 7040 2,acer laptop 3,Landauer 4,gold-c01
rk4=1;%1:rk4,0:heun Method,2:4th predictor-corrector
paramselec=1;%1:2017PRBr-SeHyeok 0:2011nature-I.Radu
gpusave=5e-12;%how often saving gpu data
debugg=0;
loadstartm=0;%1:load mat file; 0:direct calculate
%gpuDevice(1)
%constant
if pcs==1
    addpath('E:\dropbox\Dropbox\phd\code\general\gitcontrol\constant');
elseif pcs==2
    addpath('D:\Dropbox\phd\code\general\gitcontrol\constant');
    addpath(genpath('D:\Dropbox\phd\code\general\nogitcontrol'));
elseif pcs==3
    addpath('/home/a0132576/code/general/general/constant');
else
    addpath('/home/svu/a0132576/code/general/constant');
end
constantfile;
clear gam
%unknown params
if paramselec
    %     muRE=7.63*mub;%[J/T], magnetic moment, refer to my paper "Theoretical
    %     %Proposal for Determining Angular Momentum Compensation in Ferrimagnets"
    %     muTM=2.217*mub;%[J/T], magnetic moment
else
    d=0.4e-9;%[m],distance between two spin, no use when no dipolar field
    if DMIenable
        Dsim=25e-6*ele;%[J], DMI
    else
        Dsim=0;
    end
end
%params from paper
if paramselec
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
else
    muRE=7.63*mub;%[J/T]
    muTM=1.92*mub;
    Asim=1.09e-21;%[J]
    Ksim=0.807246e-23;%[J], easy-axis anisotropy
    kksim=0;%[J], DW hard-axis anisotropy
    alp=0.05;%Gilbert damping
    gTM=2;gRE=2;%g-factor
    tsteppaper=0.1e-15;%[s]
end
%% electrical parameters
jc=0;%[A/m2]
jc2=1e5;%[A/m2]
Hext=[0,0,0e-3];%corresponding to runtime
Hext2=[0,0,400e-3];% corresponding to runtime2
%% FiM parameters
tz=d;%[m],thickness of FiM
mmu=3.45*mub;%[J/T],refer to 2016-Antiferromagnetic Domain Wall Motion Driven by Spin-Orbit Torques-PRL-Takayuki Shiino
ms=mmu/d^3;%[A/m], saturation magnetization
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
BD=hbar/2*thetaSH*jc/(ms*tz);%[T]
BD2=hbar/2*thetaSH*jc2/(ms*tz);
BF=chi*BD;
BF2=chi*BD2;
%% other parameters
%natom=1000;
natom=10;
dmdt_stop_count=20;%continuously 20 times is recgnized as complete.
dmdt_stop_count_tmp=0;
TA=201.66;%from Fig. 1(b) in paper
TA=149.4642;%calculated from muTM,muRE,gamTM,gamRE in code
T_=[108.6484,130.3781,152.1078,173.8375,217.2968,239.0265,260.7562,282.4859];
%corresponds to delta_s=[-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8]e-7 [J.s/m^3]
T=108;%[K]
%T=260.7562;
muTM=(-0.6855*T+1242)*(1e3)*d^3;%[A.m^2=J/T]
muRE=(-0.462*T+1105)*(1e3)*d^3;%[A.m^2=J/T]
delta_sy=0.009204*T-1.8;%[e-7J.s/m^3]
if paramselec
    tstep=2e-15;
    if loadstartm
        runtime=0;
    else
        runtime=2*gpusave;%first run for relaxation
    end
    runtime2=1*gpusave;%second run for dw motion
    dmdt_stop=1e-6;%reference value when relaxation completes.
else
    runtime=80e-14;%first run for relaxation
    runtime2=0e-14;%second run for dw motion
    %tstep=tsteppaper;
    tstep=4e-18;
    dmdt_stop=1e-8;%reference value when relaxation completes.
end
savetstep=400;%to reduce data size
gpusteps=round(gpusave/tstep);
bc=1;%0.periodic condition;1,not periodic
%calc
if runtime2==0
    relaxflag=0;%0:no relax, 1: have relax
else
    relaxflag=1;
end
gamTM=gTM*mub/(hbar*ele);%1/(s.T)refer to "PRL 97, 217202 (2006), Jiang Xin"
gamRE=gRE*mub/(hbar*ele);
totstep1=round(runtime/tstep);
totstep2=round(runtime2/tstep);
m_=zeros(3,natom);%initial magnetization
mark_=0.5*ones(1,natom);%used as a mark
loc_=linspace(0,(natom-1)*d,natom);%atom location
if (SOT_DLT || SOT_FLT || STT_DLT || STT_FLT) && ~(rk4==1)
    error('only rk4 is implemented for spin torque driven')
end
if ~(gpuc==1)
    error('only gpu version is implemented')
end
if (1)%
    for ct=1:natom/2%initization for left half spin
        if mod(ct,2)==1%the atom is TM
            thet_=10/180*pi;
            phi_=0;
            m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
            mark_(ct)=1;
        else
            thet_=(10+180)/180*pi;
            phi_=0;
            m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
            mark_(ct)=0;
        end
        
    end
    for ct=natom/2+1:natom%initization for right half spin
        if mod(ct,2)==1%the atom is TM
            thet_=(10+180)/180*pi;
            phi_=0;
            m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
            mark_(ct)=1;
        else
            thet_=10/180*pi;
            phi_=0;
            m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
            mark_(ct)=0;
        end
        
    end
else
    for ct=1:natom
        if mod(ct,2)==1%the atom is TM
            thet_=10/180*pi;
            phi_=0;
            m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
            mark_(ct)=1;
        else
            thet_=(10+180)/180*pi;
            phi_=0;
            m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
            mark_(ct)=0;
        end
    end
end
if loadstartm
    clear m_
    load('startm_natom1000.mat')
    m_=[mmxstart;mmystart;mmzstart];
end
%clear ct
%dynamic calc
rk4_4llg(); toc
if pcs==3 || pcs==4
    save('finalgpu_tt8.mat')
else
    %save('finalgpu.mat')
end
%plot_proc();

%save('t1.mat')





