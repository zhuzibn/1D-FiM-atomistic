%LLG solver for gpu calculation using Heun Method
% "clear xx" is not allowed in gpu version of arrayfun
function [sxx,syy,szz]=atomgpu(ssx,ssy,ssz,scal,alph,ts,hhx,hhy,hhz)
%cross(u,v)=(u2v3-u3v2)i+(u3v1-u1v3)j+(u1v2-u2v1)k
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=ssx;v2=ssy;v3=ssz;
dsdt1x=u2*v3-u3*v2;
dsdt1y=u3*v1-u1*v3;
dsdt1z=u1*v2-u2*v1;
%cross(cross(sss,hh),sss) Beff DLT
u1=-dsdt1x;u2=-dsdt1y;u3=-dsdt1z;
v1=ssx;v2=ssy;v3=ssz;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%
dsdtx=dsdt1x+alph*dsdt2x;
dsdty=dsdt1y+alph*dsdt2y;
dsdtz=dsdt1z+alph*dsdt2z;
%
kk1x=ts*scal*dsdtx;kk1y=ts*scal*dsdty;kk1z=ts*scal*dsdtz;
%sss=ss1+ts*kk1;%y[i+1]
sxtmp=ssx+kk1x;
sytmp=ssy+kk1y;
sztmp=ssz+kk1z;
%cross(u,v)=(u2v3-u3v2)i+(u3v1-u1v3)j+(u1v2-u2v1)k
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=sxtmp;v2=sytmp;v3=sztmp;
dsdt1x=u2*v3-u3*v2;
dsdt1y=u3*v1-u1*v3;
dsdt1z=u1*v2-u2*v1;
%cross(cross(sss,hh),sss) Beff DLT
u1=-dsdt1x;u2=-dsdt1y;u3=-dsdt1z;
v1=sxtmp;v2=sytmp;v3=sztmp;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%
dsdtx=dsdt1x+alph*dsdt2x;
dsdty=dsdt1y+alph*dsdt2y;
dsdtz=dsdt1z+alph*dsdt2z;
%
kk2x=ts*scal*dsdtx;kk2y=ts*scal*dsdty;kk2z=ts*scal*dsdtz;

snx=ssx+1/2*(kk1x+kk2x);
sny=ssy+1/2*(kk1y+kk2y);
snz=ssz+1/2*(kk1z+kk2z);
normsn=sqrt(snx^2+sny^2+snz^2);
snx=snx/normsn;
sny=sny/normsn;
snz=snz/normsn;

sxx=snx;syy=sny;szz=snz;
end