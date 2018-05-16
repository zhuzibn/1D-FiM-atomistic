% data process
clear all;clc;close all
if (1)%Eq.5&6 with DMI=0
    addpath('D:\Dropbox\phd\code\general\gitcontrol\constant');
    constantfile;
    clear gam
    %parameters from paper
    d=0.4e-9;%[m]
    kappaa=2*0.2e-6*ele/d^3;%[J/m3]
    A=(4*7.5e-3)/d;%[eV/m]
    K=(2*0.4e-3)/d^3;%[eV/m3]
    alph=0.002;
    tselc=1;%1:140K,2:149.4642K,3:160K 4:T sweep 5:test
    switch tselc
        case 1
            T_=130;%[K]
        case 2
            T_=149.4642;%[K]
        case 3
            T_=170;%[K]
        case 4
            T_=[149.1:0.1:149.9];
        case 5
            T_=140;
    end
    Hwb_=zeros(size(T_,2),1);
    for ctT=1:size(T_,2)
        T=T_(ctT);
        msTM=(-0.6855*T+1242)*(1e3);%[A/m]
        msRE=(-0.462*T+1105)*(1e3);%[A/m]
        gTM=2.2;gRE=2;
        %calc
        lambd=sqrt(A/K);%[m]
        gamTM=gTM*mub/(hbar*ele);%1/(s.T)refer to "PRL 97, 217202 (2006), Jiang Xin"
        gamRE=gRE*mub/(hbar*ele);
        sTM=msTM/gamTM;%[AsT/m]
        sRE=msRE/gamRE;%[AsT/m]
        Mnet=msTM-msRE;
        if (1)
            H_=[1:1:15]*1e-3;%[T]
            Vdw_=zeros(size(H_,2),1);
            for ct1=1:size(H_,2)
                H=H_(ct1);
                Vdw_(ct1)=Mnet*lambd*H/(alph*(sTM+sRE));%[m/s] eq.5
            end
        end
        Hwb_(ctT)=kappaa*alph*(sTM+sRE)/(2*(sTM-sRE)*Mnet)*1000%[mT] Eq.6
    end
end