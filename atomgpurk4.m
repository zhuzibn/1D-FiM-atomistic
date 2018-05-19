%LLG solver for gpu calculation using Heun Method
% "clear xx" is not allowed in gpu version of arrayfun
function [sxx,syy,szz]=atomgpurk4(ssx,ssy,ssz,psjSHEx,psjSHEy,psjSHEz,scal,alph,ts,hhx,hhy,hhz,Bd,Bf)
%cross(u,v)=(u2v3-u3v2)i+(u3v1-u1v3)j+(u1v2-u2v1)k
%------------------kk1--------------------------
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
kk1x=scal*dsdtx;kk1y=scal*dsdty;kk1z=scal*dsdtz;
%-------------------kk2------------------------
%sss=ss1+ts/2*kk1;%y[i+1]
sxtmp=ssx+ts/2*kk1x;
sytmp=ssy+ts/2*kk1y;
sztmp=ssz+ts/2*kk1z;
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
%cross(-sss,cross(sss,ey)) Bd DLT
u1tmp=sxtmp;u2tmp=sytmp;u3tmp=sztmp;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-sxtmp;u2=-sytmp;u3=-sztmp;
dbdx=u2*v3-u3*v2;
dbdy=u3*v1-u1*v3;
dbdz=u1*v2-u2*v1;
%cross(sss,ey) Bd FLT
u1=sxtmp;u2=sytmp;u3=sztmp;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdx=u2*v3-u3*v2;
fbdy=u3*v1-u1*v3;
fbdz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) Bf DLT
u1tmp=sxtmp;u2tmp=sytmp;u3tmp=sztmp;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=sxtmp;u2=sytmp;u3=sztmp;
dbfx=u2*v3-u3*v2;
dbfy=u3*v1-u1*v3;
dbfz=u1*v2-u2*v1;
%cross(sss,ey) Bf FLT
u1=sxtmp;u2=sytmp;u3=sztmp;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfx=u2*v3-u3*v2;
fbfy=u3*v1-u1*v3;
fbfz=u1*v2-u2*v1;
%
dsdtx=dsdt1x+alph*dsdt2x+Bd*dbdx+alph*Bd*fbdx+alph*Bf*dbfx+Bf*fbfx;
dsdty=dsdt1y+alph*dsdt2y+Bd*dbdy+alph*Bd*fbdy+alph*Bf*dbfy+Bf*fbfy;
dsdtz=dsdt1z+alph*dsdt2z+Bd*dbdz+alph*Bd*fbdz+alph*Bf*dbfz+Bf*fbfz;
%
kk2x=scal*dsdtx;kk2y=scal*dsdty;kk2z=scal*dsdtz;
%-------------------kk3------------------------
%sss=ss1+ts/2*kk1;%y[i+1]
sxtmp2=ssx+ts/2*kk2x;
sytmp2=ssy+ts/2*kk2y;
sztmp2=ssz+ts/2*kk2z;

%cross(u,v)=(u2v3-u3v2)i+(u3v1-u1v3)j+(u1v2-u2v1)k
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=sxtmp2;v2=sytmp2;v3=sztmp2;
dsdt1x=u2*v3-u3*v2;
dsdt1y=u3*v1-u1*v3;
dsdt1z=u1*v2-u2*v1;
%cross(cross(sss,hh),sss) Beff DLT
u1=-dsdt1x;u2=-dsdt1y;u3=-dsdt1z;
v1=sxtmp2;v2=sytmp2;v3=sztmp2;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) Bd DLT
u1tmp=sxtmp2;u2tmp=sytmp2;u3tmp=sztmp2;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-sxtmp2;u2=-sytmp2;u3=-sztmp2;
dbdx=u2*v3-u3*v2;
dbdy=u3*v1-u1*v3;
dbdz=u1*v2-u2*v1;
%cross(sss,ey) Bd FLT
u1=sxtmp2;u2=sytmp2;u3=sztmp2;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdx=u2*v3-u3*v2;
fbdy=u3*v1-u1*v3;
fbdz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) Bf DLT
u1tmp=sxtmp2;u2tmp=sytmp2;u3tmp=sztmp2;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=sxtmp2;u2=sytmp2;u3=sztmp2;
dbfx=u2*v3-u3*v2;
dbfy=u3*v1-u1*v3;
dbfz=u1*v2-u2*v1;
%cross(sss,ey) Bf FLT
u1=sxtmp2;u2=sytmp2;u3=sztmp2;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfx=u2*v3-u3*v2;
fbfy=u3*v1-u1*v3;
fbfz=u1*v2-u2*v1;
%
dsdtx=dsdt1x+alph*dsdt2x+Bd*dbdx+alph*Bd*fbdx+alph*Bf*dbfx+Bf*fbfx;
dsdty=dsdt1y+alph*dsdt2y+Bd*dbdy+alph*Bd*fbdy+alph*Bf*dbfy+Bf*fbfy;
dsdtz=dsdt1z+alph*dsdt2z+Bd*dbdz+alph*Bd*fbdz+alph*Bf*dbfz+Bf*fbfz;
%
kk3x=scal*dsdtx;kk3y=scal*dsdty;kk3z=scal*dsdtz;
%-------------------kk4------------------------
%sss=ss1+ts/2*kk1;%y[i+1]
sxtmp3=ssx+ts*kk3x;
sytmp3=ssy+ts*kk3y;
sztmp3=ssz+ts*kk3z;

%cross(u,v)=(u2v3-u3v2)i+(u3v1-u1v3)j+(u1v2-u2v1)k
%-cross(sss,hh)=cross(hh,sss) Beff FLT
u1=hhx;u2=hhy;u3=hhz;
v1=sxtmp3;v2=sytmp3;v3=sztmp3;
dsdt1x=u2*v3-u3*v2;
dsdt1y=u3*v1-u1*v3;
dsdt1z=u1*v2-u2*v1;
%cross(cross(sss,hh),sss) Beff DLT
u1=-dsdt1x;u2=-dsdt1y;u3=-dsdt1z;
v1=sxtmp3;v2=sytmp3;v3=sztmp3;
dsdt2x=u2*v3-u3*v2;
dsdt2y=u3*v1-u1*v3;
dsdt2z=u1*v2-u2*v1;
%cross(-sss,cross(sss,ey)) Bd DLT
u1tmp=sxtmp3;u2tmp=sytmp3;u3tmp=sztmp3;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=-sxtmp3;u2=-sytmp3;u3=-sztmp3;
dbdx=u2*v3-u3*v2;
dbdy=u3*v1-u1*v3;
dbdz=u1*v2-u2*v1;
%cross(sss,ey) Bd FLT
u1=sxtmp3;u2=sytmp3;u3=sztmp3;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbdx=u2*v3-u3*v2;
fbdy=u3*v1-u1*v3;
fbdz=u1*v2-u2*v1;
%cross(sss,cross(sss,ey)) Bf DLT
u1tmp=sxtmp3;u2tmp=sytmp3;u3tmp=sztmp3;
v1tmp=psjSHEx;v2tmp=psjSHEy;v3tmp=psjSHEz;
v1=u2tmp*v3tmp-u3tmp*v2tmp;
v2=u3tmp*v1tmp-u1tmp*v3tmp;
v3=u1tmp*v2tmp-u2tmp*v1tmp;
u1=sxtmp3;u2=sytmp3;u3=sztmp3;
dbfx=u2*v3-u3*v2;
dbfy=u3*v1-u1*v3;
dbfz=u1*v2-u2*v1;
%cross(sss,ey) Bf FLT
u1=sxtmp3;u2=sytmp3;u3=sztmp3;
v1=psjSHEx;v2=psjSHEy;v3=psjSHEz;
fbfx=u2*v3-u3*v2;
fbfy=u3*v1-u1*v3;
fbfz=u1*v2-u2*v1;
%
dsdtx=dsdt1x+alph*dsdt2x+Bd*dbdx+alph*Bd*fbdx+alph*Bf*dbfx+Bf*fbfx;
dsdty=dsdt1y+alph*dsdt2y+Bd*dbdy+alph*Bd*fbdy+alph*Bf*dbfy+Bf*fbfy;
dsdtz=dsdt1z+alph*dsdt2z+Bd*dbdz+alph*Bd*fbdz+alph*Bf*dbfz+Bf*fbfz;
%
kk4x=scal*dsdtx;kk4y=scal*dsdty;kk4z=scal*dsdtz;
%-----------------final---------------------
snx=ssx+ts/6*(kk1x+2*kk2x+2*kk3x+kk4x);
sny=ssy+ts/6*(kk1y+2*kk2y+2*kk3y+kk4y);
snz=ssz+ts/6*(kk1z+2*kk2z+2*kk3z+kk4z);
normsn=sqrt(snx^2+sny^2+snz^2);
snx=snx/normsn;
sny=sny/normsn;
snz=snz/normsn;

sxx=snx;syy=sny;szz=snz;
end