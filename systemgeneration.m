m_=zeros(3,natom);%initial magnetization
mark_=0.5*ones(1,natom);%used as a mark
muigpu=zeros(1,natom,'gpuArray');
scalgpu=zeros(1,natom,'gpuArray');
BD=zeros(1,natom,'gpuArray');
switch systemselec
    case 1
        for ct=1:natom/2%initization for left half spin
            if mod(ct,2)==1%the atom is TM
                thet_=5/180*pi;
                phi_=0;
                m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
                mark_(ct)=1;
            else
                thet_=(5+180)/180*pi;
                phi_=0;
                m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
                mark_(ct)=0;
            end
        end
        for ct=natom/2+1:natom%initization for right half spin
            if mod(ct,2)==1%the atom is TM
                thet_=(5+180)/180*pi;
                phi_=0;
                m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
                mark_(ct)=1;
            else
                thet_=5/180*pi;
                phi_=0;
                m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
                mark_(ct)=0;
            end
        end      
    case 2
        for ct=1:natom/2%initization for left half spin
            thet_=5/180*pi;
            phi_=0;
            m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
            mark_(ct)=1;
        end
        for ct=natom/2+1:natom%initization for right half spin
            thet_=(5+180)/180*pi;
            phi_=0;
            m_(:,ct)=[sin(thet_)*cos(phi_),sin(thet_)*sin(phi_),cos(thet_)];
            mark_(ct)=1;
        end
    case 3
        
end
for ct2=1:natom
    muigpu(ct2)=(mark_(ct2)==1)*muTM+(mark_(ct2)==0)*muRE;
    scalgpu(ct2)=((mark_(ct2)==1)*gamTM+(mark_(ct2)==0)*gamRE)/(1+alp^2);%scale parameter
    BD(ct2)=hbar/2*thetaSH*jc/(((mark_(ct2)==1)*msTM+(mark_(ct2)==0)*msRE)*tz);%[T]
end
clear ct2
BF=chi*BD;