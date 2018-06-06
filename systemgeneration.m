m_=zeros(3,natom);%initial magnetization
mark_=0.5*ones(1,natom);%used as a mark
muigpu=zeros(1,natom,'gpuArray');
scalgpu=zeros(1,natom,'gpuArray');
BD=zeros(1,natom,'gpuArray');
Asim_next=zeros(1,natom,'gpuArray');
Asim_previous=zeros(1,natom,'gpuArray');
thet_TM_left=5/180*pi;
thet_RE_left=(5+180)/180*pi;
thet_TM_right=(5+180)/180*pi;
thet_RE_right=5/180*pi;
phi_=0;
switch systemselec
    case 1
        for ct=1:natom/2%initization for left half spin
            if mod(ct,2)==1%the atom is TM
                m_(:,ct)=[sin(thet_TM_left)*cos(phi_),sin(thet_TM_left)*sin(phi_),cos(thet_TM_left)];
                mark_(ct)=1;
            else
                m_(:,ct)=[sin(thet_RE_left)*cos(phi_),sin(thet_RE_left)*sin(phi_),cos(thet_RE_left)];
                mark_(ct)=0;
            end
        end
        for ct=natom/2+1:natom%initization for right half spin
            if mod(ct,2)==1%the atom is TM
                m_(:,ct)=[sin(thet_TM_right)*cos(phi_),sin(thet_TM_right)*sin(phi_),cos(thet_TM_right)];
                mark_(ct)=1;
            else
                m_(:,ct)=[sin(thet_RE_right)*cos(phi_),sin(thet_RE_right)*sin(phi_),cos(thet_RE_right)];
                mark_(ct)=0;
            end
        end      
    case 2
        for ct=1:natom/2%initization for left half spin
            thet_=5/180*pi;
            m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
            mark_(ct)=1;
        end
        for ct=natom/2+1:natom%initization for right half spin
            thet_=(5+180)/180*pi;
            m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
            mark_(ct)=1;
        end
    case 3
        rng('shuffle');
        tmp=randperm(natom,round(natom*compositionn));
        mark_=10*ones(natom,1);
        for ct=1:natom        
                mark_(ct)=~(ismember(ct,tmp));
        end
        clear tmp
        for ct=1:natom/2%initization for left half spin
            if mark_(ct)==1%the atom is TM
                m_(:,ct)=[sin(thet_TM_left)*cos(phi_),sin(thet_TM_left)*sin(phi_),cos(thet_TM_left)];
            else
                m_(:,ct)=[sin(thet_RE_left)*cos(phi_),sin(thet_RE_left)*sin(phi_),cos(thet_RE_left)];
            end
        end
        for ct=natom/2+1:natom%initization for right half spin
            if mark_(ct)==1%the atom is TM
                m_(:,ct)=[sin(thet_TM_right)*cos(phi_),sin(thet_TM_right)*sin(phi_),cos(thet_TM_right)];
            else
                m_(:,ct)=[sin(thet_RE_right)*cos(phi_),sin(thet_RE_right)*sin(phi_),cos(thet_RE_right)];
            end
        end 
end
clear ct
for ct2=1:natom
    muigpu(ct2)=(mark_(ct2)==1)*muTM+(mark_(ct2)==0)*muRE;
    scalgpu(ct2)=((mark_(ct2)==1)*gamTM+(mark_(ct2)==0)*gamRE)/(1+alp^2);%scale parameter
    BD(ct2)=hbar/2*thetaSH*jc/(((mark_(ct2)==1)*msTM+(mark_(ct2)==0)*msRE)*tz);%[T]
    if mark_(ct2)==1%local atom is TM
        if ct2==natom
            Asim_next(ct2)=0;
        else
            if mark_(ct2+1)==1%next atom is TM
                Asim_next(ct2)=A_TMTM;
            else
                Asim_next(ct2)=A_TMRE;
            end
        end
        
        if ct2==1
            Asim_previous(ct2)=0;
        else
            if mark_(ct2-1)==1%previous atom is TM
                Asim_previous(ct2)=A_TMTM;
            else
                Asim_previous(ct2)=A_TMRE;
            end
        end
    else%local atom is RE
        if ct2==natom
            Asim_next(ct2)=0;
        else
            if mark_(ct2+1)==1%next atom is TM
                Asim_next(ct2)=A_TMRE;
            else
                Asim_next(ct2)=A_RERE;
            end
        end
        
        if ct2==1
            Asim_previous(ct2)=0;
        else
            if mark_(ct2-1)==1%previous atom is TM
                Asim_previous(ct2)=A_TMRE;
            else
                Asim_previous(ct2)=A_RERE;
            end
        end
    end
end
clear ct2
BF=chi*BD;