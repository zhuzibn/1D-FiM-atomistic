t=linspace(tstep,runtime,totstep);
mmx_=zeros(totstep,natom);
mmy_=zeros(totstep,natom);
mmz_=zeros(totstep,natom);
muigpu=zeros(1,natom,'gpuArray');
scalgpu=zeros(1,natom,'gpuArray');
BD=zeros(1,natom,'gpuArray');
for ct2=1:natom
    muigpu(ct2)=(mod(ct2,2)==1)*muTM+(mod(ct2,2)==0)*muRE;
    scalgpu(ct2)=((mod(ct2,2)==1)*gamTM+(mod(ct2,2)==0)*gamRE)/(1+alp^2);%scale parameter
    BD(ct2)=hbar/2*thetaSH*jc/(((mod(ct2,2)==1)*msTM+(mod(ct2,2)==0)*msRE)*tz);%[T]
end
BF=chi*BD;
clear ct2
ct3=1;
ct3max=round((runtime)/gpusave);
while ~(ct3>ct3max)
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
        
        hext_x=Hext(1)*ones(size(hex_x,1),size(hex_x,2),'gpuArray');
        hext_y=Hext(2)*ones(size(hex_x,1),size(hex_x,2),'gpuArray');
        hext_z=Hext(3)*ones(size(hex_x,1),size(hex_x,2),'gpuArray');
        
        hhx=hex_x+hani_x+hdmi_x+hDWani_x+hext_x;
        hhy=hex_y+hani_y+hdmi_y+hDWani_y+hext_y;
        hhz=hex_z+hani_z+hdmi_z+hDWani_z+hext_z;
        if rk4==2%4th predictor-corrector
            if ct3==1 && ~(ct1>3)
                [sxx,syy,szz]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,0,0,0,scalgpu,alp,...
                    tstep,hhx,hhy,hhz,0,0);
            else
                [sxx,syy,szz]=arrayfun(@atomgpupc4,tmpxn0,tmpyn0,tmpzn0,...
                    tmpxn1,tmpyn1,tmpzn1,tmpxn2,tmpyn2,tmpzn2,tmpxn3,tmpyn3,tmpzn3,...
                    scalgpu,alp,tstep,hhx,hhy,hhz);
            end
        elseif rk4==1 %rk4
            [sxx,syy,szz]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,psjSHEy,psjSHEz,scalgpu,alp,...
                tstep,hhx,hhy,hhz,BD,BF);
        else%heun
            [sxx,syy,szz]=arrayfun(@atomgpu,mmxtmp,mmytmp,mmztmp,scalgpu,alp,...
                tstep,hhx,hhy,hhz);%
        end
        
        mmx(ct1+1,:)=sxx; mmy(ct1+1,:)=syy; mmz(ct1+1,:)=szz;
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
tt=t;
