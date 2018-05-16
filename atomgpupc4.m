%LLG solver for gpu calculation using 4th Predictor-corrector method
% "clear xx" is not allowed in gpu version of arrayfun
function [sxx,syy,szz]=atomgpupc4(ssx,ssy,ssz,ssxn1,ssyn1,sszn1,ssxn2,ssyn2,sszn2,ssxn3,ssyn3,sszn3,scal,alph,ts,hhx,hhy,hhz)
%cross(u,v)=(u2v3-u3v2)i+(u3v1-u1v3)j+(u1v2-u2v1)k
%% predictor
%------------------f(ti)--------------------------
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
kk1x=scal*dsdtx;kk1y=scal*dsdty;kk1z=scal*dsdtz;
%------------------f(t_{i-1})--------------------------
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=ssxn1;v2=ssyn1;v3=sszn1;
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
kk1xn1=scal*dsdtx;kk1yn1=scal*dsdty;kk1zn1=scal*dsdtz;
%------------------f(t_{i-2})--------------------------
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=ssxn2;v2=ssyn2;v3=sszn2;
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
kk1xn2=scal*dsdtx;kk1yn2=scal*dsdty;kk1zn2=scal*dsdtz;
%------------------f(t_{i-3})--------------------------
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=ssxn3;v2=ssyn3;v3=sszn3;
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
kk1xn3=scal*dsdtx;kk1yn3=scal*dsdty;kk1zn3=scal*dsdtz;
%-----------------predictor final---------------------
ssxp1=ssx+ts/24*(55*kk1x-59*kk1xn1+37*kk1xn2-3*kk1xn3);
ssyp1=ssy+ts/24*(55*kk1y-59*kk1yn1+37*kk1yn2-3*kk1yn3);
sszp1=ssz+ts/24*(55*kk1z-59*kk1zn1+37*kk1zn2-3*kk1zn3);
%% corrector
%------------------f(t_{i+1})--------------------------
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=ssxp1;v2=ssyp1;v3=sszp1;
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
kk1xp1=scal*dsdtx;kk1yp1=scal*dsdty;kk1zp1=scal*dsdtz;
%-----------------corrector final---------------------
snx=ssx+ts/24*(9*kk1xp1+19*kk1x-5*kk1xn1+kk1xn2);
sny=ssy+ts/24*(9*kk1yp1+19*kk1y-5*kk1yn1+kk1yn2);
snz=ssz+ts/24*(9*kk1zp1+19*kk1z-5*kk1zn1+kk1zn2);
normsn=sqrt(snx^2+sny^2+snz^2);
snx=snx/normsn;
sny=sny/normsn;
snz=snz/normsn;

sxx=snx;syy=sny;szz=snz;
end