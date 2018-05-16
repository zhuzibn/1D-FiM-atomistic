%% 4th order Runge Kutta method, for LLG calaulation in both IMA, PMA MTJ with FLT and DLT
% usage: add path which contain this file, call the function
% don't create the same function in new project

% call this function using: [,]=rk4_4llg(,)

%zzf,March.18,19.2016;
%1.changed based on PMA;2.add in FL torque
%% input
% Demag_, 3 by 3 matrix
% tstep is time step, unit [s]
% totstep is total number of steps
% m_init is initial magnetization, it is a 1-by-3 matrix, unit vector
% Ms: saturation magnetization, unit [emu/cm3]
% Hk: uniaxial anisotropy field, one value unit [tesla]
% Hext: applied field, 1-by-3 vector, unit [tesla]
% alp: damping constant
% P: polarization of FL and PL, currently only support same for both layer

% psj: unit 1-by-3 vector, spin flux polarization,
% note in STT the reflection type is opposite to m_pin_layer

% dimension FL_length,FL_width,FL_thickness, unit [nm]

%% output
%mmx,mmy,mmz: magnetization component, unit vector
%tt: simulation time list, unit [ns]
%Icri: critical current for switching unit:[Ampere]
if dimensionlessLLG
    Hk_=Hk;
    Hk=[1*(FL_width<FL_length)*Hk,1*(FL_width>FL_length)*Hk,0];
    tau_c=(g*Hk(2*(FL_width>FL_length)+1*(FL_width<FL_length)))/(1+alp^2); %natural time constant 1/s
    scal=1;
else
    %Hk=[1,1,1];%normalization purpose
    tau_c=1;
    %scal=gam/(1+alp^2);%scale parameter
end
ts1=tstep*tau_c; %time step
t=[linspace(tstep,runtime,totstep1),linspace(runtime+tstep,runtime+runtime2,totstep2)];

mmx_=zeros(round((runtime+runtime2)/tstep),natom);
mmy_=zeros(round((runtime+runtime2)/tstep),natom);
mmz_=zeros(round((runtime+runtime2)/tstep),natom);

muigpu=zeros(1,natom,'gpuArray');
scalgpu=zeros(1,natom,'gpuArray');
for ct2=1:natom
    muigpu(ct2)=(mod(ct2,2)==1)*muTM+(mod(ct2,2)==0)*muRE;
    scalgpu(ct2)=((mod(ct2,2)==1)*gamTM+(mod(ct2,2)==0)*gamRE)/(1+alp^2);%scale parameter
end
clear ct2

ct3=1;
ct3relax=round((runtime)/gpusave);
ct3run=round((runtime2)/gpusave);
ct3max=ct3relax+ct3run;
while ~(ct3>ct3max)
    if ct3>ct3relax
        ct0=2;
    else
        ct0=1;
    end
    %tmp=((ct0==1)*runtime+(ct0==2)*runtime2)-(ct3-1)*gpusave;
    mmx=zeros(gpusteps,natom,'gpuArray');
    mmy=zeros(gpusteps,natom,'gpuArray');
    mmz=zeros(gpusteps,natom,'gpuArray');
    
    if ~(ct3==1)
        mmx(1,:)=tmp2xn0;mmy(1,:)=tmp2yn0;mmz(1,:)=tmp2zn0;
    else
        mmx(1,:)=m_(1,:);mmy(1,:)=m_(2,:);mmz(1,:)=m_(3,:);
    end
    clear tmpx tmpy tmpz
    ct1=1; %count 1
    while ct1<gpusteps
        
        if gpuc%gpu calculation
            mmxtmp=mmx(ct1,:);
            mmytmp=mmy(ct1,:);
            mmztmp=mmz(ct1,:);
            if bc%not periodic condition
                mmxnext=circshift(mmxtmp,[0,-1]);
                mmynext=circshift(mmytmp,[0,-1]);
                mmznext=circshift(mmztmp,[0,-1]);
                mmxprevious=circshift(mmxtmp,[0,1]);
                mmyprevious=circshift(mmytmp,[0,1]);
                mmzprevious=circshift(mmztmp,[0,1]);
                mmxnext(end)=0;mmynext(end)=0;mmznext(end)=0;
                mmxprevious(1)=0;mmyprevious(1)=0;mmzprevious(1)=0;
            else%periodic condition
                mmxnext=circshift(mmxtmp,[0,-1]);
                mmynext=circshift(mmytmp,[0,-1]);
                mmznext=circshift(mmztmp,[0,-1]);
                mmxprevious=circshift(mmxtmp,[0,1]);
                mmyprevious=circshift(mmytmp,[0,1]);
                mmzprevious=circshift(mmztmp,[0,1]);
            end
            
            hex_x=-Asim./muigpu.*(mmxnext+mmxprevious);
            hex_y=-Asim./muigpu.*(mmynext+mmyprevious);
            hex_z=-Asim./muigpu.*(mmznext+mmzprevious);
            
            hani_x=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');%anisotropy
            hani_y=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');
            hani_z=2*Ksim./muigpu.*mmztmp;
            
            hDWani_x=-2*kksim./muigpu.*mmxtmp;%anisotropy
            hDWani_y=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');
            hDWani_z=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');
            
            hdmi_x=-Dsim./muigpu.*(-mmznext-mmzprevious);
            hdmi_y=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');
            hdmi_z=-Dsim./muigpu.*(mmxnext+mmxprevious);
            
            hext_x=((ct0==1)*Hext(1)+(ct0==2)*Hext2(1))*ones(size(hex_x,1),size(hex_x,2),'gpuArray');
            hext_y=((ct0==1)*Hext(2)+(ct0==2)*Hext2(2))*ones(size(hex_x,1),size(hex_x,2),'gpuArray');
            hext_z=((ct0==1)*Hext(3)+(ct0==2)*Hext2(3))*ones(size(hex_x,1),size(hex_x,2),'gpuArray');
            
            hhx=hex_x+hani_x+hdmi_x+hDWani_x+hext_x;
            hhy=hex_y+hani_y+hdmi_y+hDWani_y+hext_y;
            hhz=hex_z+hani_z+hdmi_z+hDWani_z+hext_z;
            if rk4==2%4th predictor-corrector
                if ct3==1 && ~(ct1>3)
                    [sxx,syy,szz]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,scalgpu,alp,...
                        tstep,hhx,hhy,hhz);
                else   
                [sxx,syy,szz]=arrayfun(@atomgpupc4,tmpxn0,tmpyn0,tmpzn0,...
                    tmpxn1,tmpyn1,tmpzn1,tmpxn2,tmpyn2,tmpzn2,tmpxn3,tmpyn3,tmpzn3,...
                    scalgpu,alp,tstep,hhx,hhy,hhz);  
                end
            elseif rk4==1 %rk4
                [sxx,syy,szz]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,scalgpu,alp,...
                    tstep,hhx,hhy,hhz);               
            else%heun
                [sxx,syy,szz]=arrayfun(@atomgpu,mmxtmp,mmytmp,mmztmp,scalgpu,alp,...
                    tstep,hhx,hhy,hhz);%
            end
            
            mmx(ct1+1,:)=sxx; mmy(ct1+1,:)=syy; mmz(ct1+1,:)=szz;
        else%cpu calculation, not maintained.
            for ct2=1:natom
                mm1=[mmx(ct1,ct2),mmy(ct1,ct2),mmz(ct1,ct2)];
                if mod(ct2,2)==1%the atom is TM
                    mui=muTM;
                    scal=gamTM/(1+alp^2);%scale parameter
                else
                    mui=muRE;
                    scal=gamRE/(1+alp^2);%scale parameter
                end
                hcalc();
                if rk4==1 %RK4
                    mmm=mm1;
                    dmdt=LLG_solver(alp,mmm,hh,psjSHE,BD,BF);
                    kk1=scal*dmdt;
                    
                    mmm=mm1+kk1*ts1/2;
                    dmdt=LLG_solver(alp,mmm,hh,psjSHE,BD,BF);
                    kk2=scal*dmdt;
                    
                    mmm=mm1+kk2*ts1/2;
                    dmdt=LLG_solver(alp,mmm,hh,psjSHE,BD,BF);
                    kk3=scal*dmdt;
                    
                    mmm=mm1+kk3*ts1;
                    dmdt=LLG_solver(alp,mmm,hh,psjSHE,BD,BF);
                    kk4=scal*dmdt;
                    
                    mn1=mm1+ts1/6*(kk1+2*kk2+2*kk3+kk4);
                    mn1=mn1/norm(mn1);
                    mmx(ct1+1,ct2)=mn1(1);mmy(ct1+1,ct2)=mn1(2);mmz(ct1+1,ct2)=mn1(3);
                else%Heun Method
                    mmm=mm1;
                    dmdt=LLG_solver(alp,mmm,hh);
                    kk1=scal*dmdt;
                    
                    mmm=mm1+kk1*ts1;
                    dmdt=LLG_solver(alp,mmm,hh);
                    kk2=scal*dmdt;
                    
                    mn1=mm1+ts1/2*(kk1+kk2);
                    mn1=mn1/norm(mn1);
                    mmx(ct1+1,ct2)=mn1(1);mmy(ct1+1,ct2)=mn1(2);mmz(ct1+1,ct2)=mn1(3);
                end
            end
        end
        if debugg
            if ct0==1 && relaxflag
                tmpx=max(abs(mmx(ct1+1,:)-mmx(ct1,:)));
                tmpy=max(abs(mmy(ct1+1,:)-mmy(ct1,:)));
                tmpz=max(abs(mmz(ct1+1,:)-mmz(ct1,:)));
                if max([tmpx,tmpy,tmpz])<dmdt_stop
                    dmdt_stop_count_tmp=dmdt_stop_count_tmp+1;
                    if dmdt_stop_count_tmp==dmdt_stop_count
                        fprintf('relaxation completes at %i ps', ct1*tstep*1e12)
                        relaxx=1;
                        break;
                        
                    end
                else
                    dmdt_stop_count_tmp=0;%clear the counting value
                end
                clear tmpx tmpy tmpz
                if ct1==totstep-1
                    fprintf('relaxation is not completed, stops due to time limit\n')
                    relaxx=0;%0:not completed 1:completed
                end
            end
        end
        ct1=ct1+1;
        if ~(ct3==1 && ~(ct1>3)) && ct1>3
            tmpxn0=mmx(ct1,:);tmpyn0=mmy(ct1,:);tmpzn0=mmz(ct1,:);
            tmpxn1=mmx(ct1-1,:);tmpyn1=mmy(ct1-1,:);tmpzn1=mmz(ct1-1,:);
            tmpxn2=mmx(ct1-2,:);tmpyn2=mmy(ct1-2,:);tmpzn2=mmz(ct1-2,:);
            tmpxn3=mmx(ct1-3,:);tmpyn3=mmy(ct1-3,:);tmpzn3=mmz(ct1-3,:);
        elseif ~(ct3==1 && ~(ct1>3)) && ct1==2
            tmpxn0=mmx(ct1,:);tmpyn0=mmy(ct1,:);tmpzn0=mmz(ct1,:);
            tmpxn1=tmp2xn0;tmpyn1=tmp2yn0;tmpzn1=tmp2zn0;
            tmpxn2=tmp2xn1;tmpyn2=tmp2yn1;tmpzn2=tmp2zn1;
            tmpxn3=tmp2xn2;tmpyn3=tmp2yn2;tmpzn3=tmp2zn2;
        elseif ~(ct3==1 && ~(ct1>3)) && ct1==3
            tmpxn0=mmx(ct1,:);tmpyn0=mmy(ct1,:);tmpzn0=mmz(ct1,:);
            tmpxn1=mmx(ct1-1,:);tmpyn1=mmy(ct1-1,:);tmpzn1=mmz(ct1-1,:);
            tmpxn2=tmp2xn0;tmpyn2=tmp2yn0;tmpzn2=tmp2zn0;
            tmpxn3=tmp2xn1;tmpyn3=tmp2yn1;tmpzn3=tmp2zn1;
        end
    end
    tmp2xn0=mmx(end,:);tmp2yn0=mmy(end,:);tmp2zn0=mmz(end,:);
    tmp2xn1=mmx(end-1,:);tmp2yn1=mmy(end-1,:);tmp2zn1=mmz(end-1,:);
    tmp2xn2=mmx(end-2,:);tmp2yn2=mmy(end-2,:);tmp2zn2=mmz(end-2,:);
    mmx_((ct3-1)*gpusteps+1:ct3*gpusteps,:)=gather(mmx);
    mmy_((ct3-1)*gpusteps+1:ct3*gpusteps,:)=gather(mmy);
    mmz_((ct3-1)*gpusteps+1:ct3*gpusteps,:)=gather(mmz);
    ct3=ct3+1;
end

clear mmx mmy mmz tmp2xn0 tmp2yn0 tmp2zn0 tmp2xn1 tmp2yn1 tmp2zn1
clear tmp2xn2 tmp2yn2 tmp2zn2 
mmx=mmx_(1:savetstep:end,:);
mmy=mmy_(1:savetstep:end,:);
mmz=mmz_(1:savetstep:end,:);
clear mmx_ mmy_ mmz_ 
t=t(1:savetstep:end);
if dimensionlessLLG
    tt=t/tau_c*1e9;%unit[ns]
else
    tt=t;
end
