%Heff calc
if ct2==1
mmsecond=[mmx(ct1,2),mmy(ct1,2),mmz(ct1,2)];
mmend=[mmx(ct1,end),mmy(ct1,end),mmz(ct1,end)];
elseif ct2==natom
mmendn1=[mmx(ct1,end-1),mmy(ct1,end-1),mmz(ct1,end-1)];
mmfirst=[mmx(ct1,1),mmy(ct1,1),mmz(ct1,1)];
else
mmnext=[mmx(ct1,ct2+1),mmy(ct1,ct2+1),mmz(ct1,ct2+1)];
mmprevious=[mmx(ct1,ct2-1),mmy(ct1,ct2-1),mmz(ct1,ct2-1)];
end
    %exchange interaction
    if ct2==1 %for exchange interaction, assume perodic boundry condition
        hex_=-Asim/mui*(mmend+mmsecond);%[T]
    elseif ct2==natom
        hex_=-Asim/mui*(mmendn1+mmfirst);
    else
        hex_=-Asim/mui*(mmnext+mmprevious);
    end
    %easy-axis anisotropy
    hani_=2*Ksim/mui*[0,0,mmz(ct1,ct2)];
    %DW hard-axis anisotropy
    hdwani_=-2*kksim/mui*[mmx(ct1,ct2),0,0];
    %DMI, only consider nearest neighbor
    if ct2==1 %for DMI interaction, assume perodic boundry condition
        tmp1=[-mmsecond(3),0,mmsecond(1)];
        tmp2=[-mmend(3),0,mmend(1)];
        hdmi_=-Dsim/mui*(tmp1+tmp2);
    elseif ct2==natom
        tmp1=[-mmfirst(3),0,mmfirst(1)];
        tmp2=[-mmendn1(3),0,mmendn1(1)];
        hdmi_=-Dsim/mui*(tmp1+tmp2);
    else
        tmp1=[-mmnext(3),0,mmnext(1)];
        tmp2=[-mmprevious(3),0,mmprevious(1)];
        hdmi_=-Dsim/mui*(tmp1+tmp2);
    end
    %Hext
    hext_=Hext;
    
hh=hex_+hani_+hdwani_+hdmi_+hext_;

