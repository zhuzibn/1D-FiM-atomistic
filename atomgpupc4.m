%LLG solver for gpu calculation using 4th Predictor-corrector method
% "clear xx" is not allowed in gpu version of arrayfun
function [sxx,syy,szz]=atomgpupc4(ssx,ssy,ssz,ssxn1,ssyn1,sszn1,ssxn2,ssyn2,sszn2,ssxn3,ssyn3,sszn3,...
    psjSHEx,psjSHEy,psjSHEz,scal,alph,ts,hhx,hhy,hhz,Bd,Bf)
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
%cross(-sss,cross(sss,ey)) Bd DLT
u1tmp=ssx;u2tmp=ssy;u3tmp=ssz;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-ssx;u2=-ssy;u3=-ssz;
dbdx=u2*v3-u3*v2;
dbdy=u3*v1-u1*v3;
dbdz=u1*v2-u2*v1;
%cross(sss,ey) Bd FLT
u1=ssx;u2=ssy;u3=ssz;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdx=u2*v3-u3*v2;
fbdy=u3*v1-u1*v3;
fbdz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) Bf DLT
u1tmp=ssx;u2tmp=ssy;u3tmp=ssz;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=ssx;u2=ssy;u3=ssz;
dbfx=u2*v3-u3*v2;
dbfy=u3*v1-u1*v3;
dbfz=u1*v2-u2*v1;
%cross(sss,ey) Bf FLT
u1=ssx;u2=ssy;u3=ssz;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfx=u2*v3-u3*v2;
fbfy=u3*v1-u1*v3;
fbfz=u1*v2-u2*v1;
%
dsdtx=dsdt1x+alph*dsdt2x+Bd*dbdx+alph*Bd*fbdx+alph*Bf*dbfx+Bf*fbfx;
dsdty=dsdt1y+alph*dsdt2y+Bd*dbdy+alph*Bd*fbdy+alph*Bf*dbfy+Bf*fbfy;
dsdtz=dsdt1z+alph*dsdt2z+Bd*dbdz+alph*Bd*fbdz+alph*Bf*dbfz+Bf*fbfz;
%
kk1x=ts*scal*dsdtx;kk1y=ts*scal*dsdty;kk1z=ts*scal*dsdtz;
%------------------f(t_{i-1})--------------------------
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=ssxn1;v2=ssyn1;v3=sszn1;
dsdt1x=u2*v3-u3*v2;
dsdt1y=u3*v1-u1*v3;
dsdt1z=u1*v2-u2*v1;
%cross(cross(sss,hh),sss) Beff DLT
u1=-dsdt1x;u2=-dsdt1y;u3=-dsdt1z;
v1=ssxn1;v2=ssyn1;v3=sszn1;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) Bd DLT
u1tmp=ssxn1;u2tmp=ssyn1;u3tmp=sszn1;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-ssxn1;u2=-ssyn1;u3=-sszn1;
dbdx=u2*v3-u3*v2;
dbdy=u3*v1-u1*v3;
dbdz=u1*v2-u2*v1;
%cross(sss,ey) Bd FLT
u1=ssxn1;u2=ssyn1;u3=sszn1;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdx=u2*v3-u3*v2;
fbdy=u3*v1-u1*v3;
fbdz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) Bf DLT
u1tmp=ssxn1;u2tmp=ssyn1;u3tmp=sszn1;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=ssxn1;u2=ssyn1;u3=sszn1;
dbfx=u2*v3-u3*v2;
dbfy=u3*v1-u1*v3;
dbfz=u1*v2-u2*v1;
%cross(sss,ey) Bf FLT
u1=ssxn1;u2=ssyn1;u3=sszn1;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfx=u2*v3-u3*v2;
fbfy=u3*v1-u1*v3;
fbfz=u1*v2-u2*v1;
%
dsdtx=dsdt1x+alph*dsdt2x+Bd*dbdx+alph*Bd*fbdx+alph*Bf*dbfx+Bf*fbfx;
dsdty=dsdt1y+alph*dsdt2y+Bd*dbdy+alph*Bd*fbdy+alph*Bf*dbfy+Bf*fbfy;
dsdtz=dsdt1z+alph*dsdt2z+Bd*dbdz+alph*Bd*fbdz+alph*Bf*dbfz+Bf*fbfz;
%
kk1xn1=ts*scal*dsdtx;kk1yn1=ts*scal*dsdty;kk1zn1=ts*scal*dsdtz;
%------------------f(t_{i-2})--------------------------
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=ssxn2;v2=ssyn2;v3=sszn2;
dsdt1x=u2*v3-u3*v2;
dsdt1y=u3*v1-u1*v3;
dsdt1z=u1*v2-u2*v1;
%cross(cross(sss,hh),sss) Beff DLT
u1=-dsdt1x;u2=-dsdt1y;u3=-dsdt1z;
v1=ssxn2;v2=ssyn2;v3=sszn2;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) Bd DLT
u1tmp=ssxn2;u2tmp=ssyn2;u3tmp=sszn2;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-ssxn2;u2=-ssyn2;u3=-sszn2;
dbdx=u2*v3-u3*v2;
dbdy=u3*v1-u1*v3;
dbdz=u1*v2-u2*v1;
%cross(sss,ey) Bd FLT
u1=ssxn2;u2=ssyn2;u3=sszn2;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdx=u2*v3-u3*v2;
fbdy=u3*v1-u1*v3;
fbdz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) Bf DLT
u1tmp=ssxn2;u2tmp=ssyn2;u3tmp=sszn2;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=ssxn2;u2=ssyn2;u3=sszn2;
dbfx=u2*v3-u3*v2;
dbfy=u3*v1-u1*v3;
dbfz=u1*v2-u2*v1;
%cross(sss,ey) Bf FLT
u1=ssxn2;u2=ssyn2;u3=sszn2;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfx=u2*v3-u3*v2;
fbfy=u3*v1-u1*v3;
fbfz=u1*v2-u2*v1;
%
dsdtx=dsdt1x+alph*dsdt2x+Bd*dbdx+alph*Bd*fbdx+alph*Bf*dbfx+Bf*fbfx;
dsdty=dsdt1y+alph*dsdt2y+Bd*dbdy+alph*Bd*fbdy+alph*Bf*dbfy+Bf*fbfy;
dsdtz=dsdt1z+alph*dsdt2z+Bd*dbdz+alph*Bd*fbdz+alph*Bf*dbfz+Bf*fbfz;
%
kk1xn2=ts*scal*dsdtx;kk1yn2=ts*scal*dsdty;kk1zn2=ts*scal*dsdtz;
%------------------f(t_{i-3})--------------------------
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=ssxn3;v2=ssyn3;v3=sszn3;
dsdt1x=u2*v3-u3*v2;
dsdt1y=u3*v1-u1*v3;
dsdt1z=u1*v2-u2*v1;
%cross(cross(sss,hh),sss) Beff DLT
u1=-dsdt1x;u2=-dsdt1y;u3=-dsdt1z;
v1=ssxn3;v2=ssyn3;v3=sszn3;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) Bd DLT
u1tmp=ssxn3;u2tmp=ssyn3;u3tmp=sszn3;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-ssxn3;u2=-ssyn3;u3=-sszn3;
dbdx=u2*v3-u3*v2;
dbdy=u3*v1-u1*v3;
dbdz=u1*v2-u2*v1;
%cross(sss,ey) Bd FLT
u1=ssxn3;u2=ssyn3;u3=sszn3;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdx=u2*v3-u3*v2;
fbdy=u3*v1-u1*v3;
fbdz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) Bf DLT
u1tmp=ssxn3;u2tmp=ssyn3;u3tmp=sszn3;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=ssxn3;u2=ssyn3;u3=sszn3;
dbfx=u2*v3-u3*v2;
dbfy=u3*v1-u1*v3;
dbfz=u1*v2-u2*v1;
%cross(sss,ey) Bf FLT
u1=ssxn3;u2=ssyn3;u3=sszn3;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfx=u2*v3-u3*v2;
fbfy=u3*v1-u1*v3;
fbfz=u1*v2-u2*v1;
%
dsdtx=dsdt1x+alph*dsdt2x+Bd*dbdx+alph*Bd*fbdx+alph*Bf*dbfx+Bf*fbfx;
dsdty=dsdt1y+alph*dsdt2y+Bd*dbdy+alph*Bd*fbdy+alph*Bf*dbfy+Bf*fbfy;
dsdtz=dsdt1z+alph*dsdt2z+Bd*dbdz+alph*Bd*fbdz+alph*Bf*dbfz+Bf*fbfz;
%
kk1xn3=ts*scal*dsdtx;kk1yn3=ts*scal*dsdty;kk1zn3=ts*scal*dsdtz;
%-----------------predictor final---------------------
ssxp1=ssx+1/24*(55*kk1x-59*kk1xn1+37*kk1xn2-3*kk1xn3);
ssyp1=ssy+1/24*(55*kk1y-59*kk1yn1+37*kk1yn2-3*kk1yn3);
sszp1=ssz+1/24*(55*kk1z-59*kk1zn1+37*kk1zn2-3*kk1zn3);
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
v1=ssxp1;v2=ssyp1;v3=sszp1;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) Bd DLT
u1tmp=ssxp1;u2tmp=ssyp1;u3tmp=sszp1;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-ssxp1;u2=-ssyp1;u3=-sszp1;
dbdx=u2*v3-u3*v2;
dbdy=u3*v1-u1*v3;
dbdz=u1*v2-u2*v1;
%cross(sss,ey) Bd FLT
u1=ssxp1;u2=ssyp1;u3=sszp1;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdx=u2*v3-u3*v2;
fbdy=u3*v1-u1*v3;
fbdz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) Bf DLT
u1tmp=ssxp1;u2tmp=ssyp1;u3tmp=sszp1;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=ssxp1;u2=ssyp1;u3=sszp1;
dbfx=u2*v3-u3*v2;
dbfy=u3*v1-u1*v3;
dbfz=u1*v2-u2*v1;
%cross(sss,ey) Bf FLT
u1=ssxp1;u2=ssyp1;u3=sszp1;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfx=u2*v3-u3*v2;
fbfy=u3*v1-u1*v3;
fbfz=u1*v2-u2*v1;
%
dsdtx=dsdt1x+alph*dsdt2x+Bd*dbdx+alph*Bd*fbdx+alph*Bf*dbfx+Bf*fbfx;
dsdty=dsdt1y+alph*dsdt2y+Bd*dbdy+alph*Bd*fbdy+alph*Bf*dbfy+Bf*fbfy;
dsdtz=dsdt1z+alph*dsdt2z+Bd*dbdz+alph*Bd*fbdz+alph*Bf*dbfz+Bf*fbfz;
%
kk1xp1=ts*scal*dsdtx;kk1yp1=ts*scal*dsdty;kk1zp1=ts*scal*dsdtz;
%-----------------corrector final---------------------
snx=ssx+1/24*(9*kk1xp1+19*kk1x-5*kk1xn1+kk1xn2);
sny=ssy+1/24*(9*kk1yp1+19*kk1y-5*kk1yn1+kk1yn2);
snz=ssz+1/24*(9*kk1zp1+19*kk1z-5*kk1zn1+kk1zn2);
normsn=sqrt(snx^2+sny^2+snz^2);
snx=snx/normsn;
sny=sny/normsn;
snz=snz/normsn;

sxx=snx;syy=sny;szz=snz;
end