t=linspace(tstep,runtime,totstep);
mmx_=zeros(totstep,natom);
mmy_=zeros(totstep,natom);
mmz_=zeros(totstep,natom);
Eex_=zeros(totstep,1);
Eani_=zeros(totstep,1);
EDWani_=zeros(totstep,1);
Edmi_=zeros(totstep,1);

ct3=1;
ct3max=round((runtime)/gpusave);
scalgpu=((mark_==1)*gamTM+(mark_==0)*gamRE)/(1+alp^2);%scale parameter
T_=tempera(10,t*1e9,0)+T0;
[MsT_TM,MsT_RE]=MsTemper(T_);
if (0)%view the interpolated "Ms vs T"
    figure;
    plot(T_+T0,MsT_TM,'-*');
    xlabel('T(K)','Ms_{TM}(A/m)')
    figure;
    plot(T_+T0,MsT_RE,'-*');
    xlabel('T(K)','Ms_{TM}(A/m)')
end
muTM=MsT_TM*d^3;%[A.m^2=J/T]
muRE=MsT_RE*d^3;%[A.m^2=J/T]
while ~(ct3>ct3max)
    mmx=zeros(gpusteps,natom,'gpuArray');
    mmy=zeros(gpusteps,natom,'gpuArray');
    mmz=zeros(gpusteps,natom,'gpuArray');
    Eex=zeros(gpusteps,1,'gpuArray');
    Eani=zeros(gpusteps,1,'gpuArray');
    EDWani=zeros(gpusteps,1,'gpuArray');
    Edmi=zeros(gpusteps,1,'gpuArray');
    if ~(ct3==1)
        mmx(1,:)=tmp2xn0;mmy(1,:)=tmp2yn0;mmz(1,:)=tmp2zn0;
    else
        mmx(1,:)=m_(1,:);mmy(1,:)=m_(2,:);mmz(1,:)=m_(3,:);
    end
    clear tmpx tmpy tmpz
    ct1=1; %count 1
    
    while ct1<gpusteps
        muigpu=(mark_==1)*muTM((ct3-1)*gpusteps+ct1)+(mark_==0)*muRE((ct3-1)*gpusteps+ct1);
        BD=hbar/2*thetaSH*jc./(((mark_==1)*MsT_TM((ct3-1)*gpusteps+ct1)+(mark_==0)*MsT_RE((ct3-1)*gpusteps+ct1))*tz);%[T]
        BF=chi*BD;
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
        
        hex_x=-(Asim_next.*mmxnext+Asim_previous.*mmxprevious)./muigpu;
        hex_y=-(Asim_next.*mmynext+Asim_previous.*mmyprevious)./muigpu;
        hex_z=-(Asim_next.*mmznext+Asim_previous.*mmzprevious)./muigpu;
        Eex(ct1)=sum(Asim_next.*(mmxtmp.*mmxnext+mmytmp.*mmynext+mmztmp.*mmznext)+...
            Asim_previous.*(mmxtmp.*mmxprevious+mmytmp.*mmyprevious+mmztmp.*mmzprevious));
        
        hani_x=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');%anisotropy
        hani_y=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');
        hani_z=2*Ksim./muigpu.*mmztmp;
        Eani(ct1)=sum(-Ksim.*(mmztmp.^2));
        
        hDWani_x=-2*kksim./muigpu.*mmxtmp;%anisotropy
        hDWani_y=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');
        hDWani_z=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');
        EDWani(ct1)=sum(kksim.*(mmxtmp.^2));
        
        hdmi_x=-Dsim./muigpu.*(-mmznext+mmzprevious);
        hdmi_y=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');
        hdmi_z=-Dsim./muigpu.*(mmxnext-mmxprevious);
        Edmi(ct1)=sum(Dsim.*(mmztmp.*mmxnext-mmxtmp.*mmznext+...
            mmzprevious.*mmxtmp-mmxprevious.*mmztmp));
        
        hext_x=Hext(1)*ones(size(hex_x,1),size(hex_x,2),'gpuArray');
        hext_y=Hext(2)*ones(size(hex_x,1),size(hex_x,2),'gpuArray');
        hext_z=Hext(3)*ones(size(hex_x,1),size(hex_x,2),'gpuArray');
        
        hhx=hex_x+hani_x+hdmi_x+hDWani_x+hext_x;
        hhy=hex_y+hani_y+hdmi_y+hDWani_y+hext_y;
        hhz=hex_z+hani_z+hdmi_z+hDWani_z+hext_z;
        if rk4==2%4th predictor-corrector
            if ct3==1 && ~(ct1>3)
                [sxx,syy,szz]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,psjSHEy,psjSHEz,scalgpu,alp,...
                    tstep,hhx,hhy,hhz,BD,BF);
            else
                [sxx,syy,szz]=arrayfun(@atomgpupc4,tmpxn0,tmpyn0,tmpzn0,...
                    tmpxn1,tmpyn1,tmpzn1,tmpxn2,tmpyn2,tmpzn2,tmpxn3,tmpyn3,tmpzn3,...
                    psjSHEx,psjSHEy,psjSHEz,scalgpu,alp,tstep,hhx,hhy,hhz,BD,BF);
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
    Eex_((ct3-1)*gpusteps+1:ct3*gpusteps)=gather(Eex);
    Eani_((ct3-1)*gpusteps+1:ct3*gpusteps)=gather(Eani);
    EDWani_((ct3-1)*gpusteps+1:ct3*gpusteps)=gather(EDWani);
    Edmi_((ct3-1)*gpusteps+1:ct3*gpusteps)=gather(Edmi);
    ct3=ct3+1;
end
clear mmx mmy mmz tmp2xn0 tmp2yn0 tmp2zn0 tmp2xn1 tmp2yn1 tmp2zn1 Eex Eani EDWani Edmi
clear tmp2xn2 tmp2yn2 tmp2zn2
mmx=mmx_(1:savetstep:end,:);
mmy=mmy_(1:savetstep:end,:);
mmz=mmz_(1:savetstep:end,:);
Eex=Eex_(1:savetstep:end);
Eani=Eani_(1:savetstep:end);
EDWani=EDWani_(1:savetstep:end);
Edmi=Edmi_(1:savetstep:end);
clear mmx_ mmy_ mmz_ Eex_ Eani_ EDWani_ Edmi_
t=t(1:savetstep:end);
tt=t;
