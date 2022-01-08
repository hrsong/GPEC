% Copyright (c) 2021 Haoran Song
%
% This file is part of GPEC
%
% GPEC is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% GPEC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with GPEC.  If not, see <http://www.gnu.org/licenses/>.

clc
clear
% multi order
n=12;
itest3=30
ex=zeros(300,1000,18); 
ey=zeros(300,1000,18);
dex=zeros(300,1000,18);
dey=zeros(300,1000,18);
ddex=zeros(300,1000,18);
ddey=zeros(300,1000,18);
evx=zeros(300,1000,18);
evy=zeros(300,1000,18);
devx=zeros(300,1000,18);
devy=zeros(300,1000,18);
x=zeros(300,1000,18);
y=zeros(300,1000,18);
vx=zeros(300,1000,18);
vy=zeros(300,1000,18);

xcs=zeros(300,1000,18);
ycs=zeros(300,1000,18);
vxcs=zeros(300,1000,18);
vycs=zeros(300,1000,18);

parfor a=1:18 % order 8 to order 25
    n=7+a;
    C1=zeros(1,n-1);
    C2=zeros(1,n);
    C3=zeros(1,n-1);
    C4=zeros(1,n);
    c1=zeros(1,n);
    c2=zeros(1,n+1);
    c1(1)=1;
    c2(1)=1;
    for i=1:n-1
    C1(i)=coef1(i);
    C2(i)=coef2(i);
    C3(i)=coef3(i);
    C4(i)=coef4(i);
    end
    C2(n)=coef2(n);
    C4(n)=coef4(n);
    for i=1:n-1
    for j=0:i
        c1(j+1)=c1(j+1)+(-1)^j*nchoosek(i,j)*C1(i);
        c2(j+1)=c2(j+1)+(-1)^j*nchoosek(i,j)*C2(i);
    end
    end
    for j=0:n
        c2(j+1)=c2(j+1)+(-1)^j*nchoosek(n,j)*C2(n);
    end


    for t1=1:300 % t1 steps per loop
        t1
    T2=t1;
    w=2*pi/T2;
    
    mu=w*w;
    x0=0;
    y0=1;
    vx0=1;
    vy0=0;
    dx1=x0-sin(-w);
    dx2=dx1;
    dy1=y0-cos(-w);
    dy2=dy1;
    
    x1=x0;
    x2=x0;
    vx1=vx0;
    vx2=vx0;
    y1=y0;
    y2=y0;
    vy1=vy0;
    vy2=vy0;
    x10=x0*0;
    x20=x0*0;
    y10=y0*0;
    y20=y0*0;
    vx10=x0*0;
    vx20=x0*0;
    vy10=y0*0;
    vy20=y0*0;
    dx10=x0*0;
    dx20=x0*0;
    dy10=y0*0;
    dy20=y0*0;
    
    fx=zeros(n+1,1);
    fy=zeros(n+1,1);
    for i=1:n+1
        fx(:,i)=-w*w*sin(-(i-1)*w);
        fy(:,i)=-w*w*cos(-(i-1)*w);
    end
    dfx=fx;
    dfy=fy;
    for i=2:n+1
       dfx(:,i:n+1)=dfx(:,(i-1):n)-dfx(:,i:n+1);
       dfy(:,i:n+1)=dfy(:,(i-1):n)-dfy(:,i:n+1);
    end
    dfx0=dfx;
    dfy0=dfy;
    for t2=1:1000
        dx0=0;
        dx00=0;
        dy0=0;
        dy00=0;
        dvx0=0;
        dvx00=0;
        dvy0=0;
        dvy00=0;
    for t=1:T2 % loop count
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmpx=C1(n-1)*dfx(:,n);
        tmpy=C1(n-1)*dfy(:,n);
        for i=1:n-2
            tmpx=tmpx+C1(n-1-i)*dfx(:,n-i);
            tmpy=tmpy+C1(n-1-i)*dfy(:,n-i);
        end
        tmpx=tmpx+dfx(:,1);
        tmpy=tmpy+dfy(:,1);
        for j=1:1
            [tmpx(j),~]=add_cs(dx2(j),dx20(j),tmpx(j));
            [tmpy(j),~]=add_cs(dy2(j),dy20(j),tmpy(j));
            [tmpx(j),~]=add_cs(x2(j),x20(j),tmpx(j));
            [tmpy(j),~]=add_cs(y2(j),y20(j),tmpy(j));
        end
        g=grav(mu,tmpx,tmpy);
        dfx0(:,1)=g.*tmpx;
        dfy0(:,1)=g.*tmpy;
        for i=2:n+1
            dfx0(:,i)=dfx0(:,i-1)-dfx(:,i-1);
            dfy0(:,i)=dfy0(:,i-1)-dfy(:,i-1);
        end
        tmpx=C2(n)*dfx0(:,n+1);
        tmpy=C2(n)*dfy0(:,n+1);
        tmpvx=C4(n)*dfx0(:,n+1);
        tmpvy=C4(n)*dfy0(:,n+1);
        for i=1:n-1
            tmpx=tmpx+C2(n-i)*dfx0(:,n+1-i);
            tmpy=tmpy+C2(n-i)*dfy0(:,n+1-i);
            tmpvx=tmpvx+C4(n-i)*dfx0(:,n+1-i);
            tmpvy=tmpvy+C4(n-i)*dfy0(:,n+1-i);
        end
        tmpx=tmpx+dfx0(:,1);
        tmpy=tmpy+dfy0(:,1);
        tmpvx=tmpvx+dfx0(:,1);
        tmpvy=tmpvy+dfy0(:,1);
        for j=1:1
            [dx0(j),dx00(j)]=add_cs(dx0(j),dx00(j),tmpx(j));
            [dy0(j),dy00(j)]=add_cs(dy0(j),dy00(j),tmpy(j));
            [dx2(j),dx20(j)]=add_cs(dx2(j),dx20(j),tmpx(j));
            [dy2(j),dy20(j)]=add_cs(dy2(j),dy20(j),tmpy(j));
            [x2(j),x20(j)]=add_cs(x2(j),x20(j),dx2(j));
            [y2(j),y20(j)]=add_cs(y2(j),y20(j),dy2(j));
            [vx2(j),vx20(j)]=add_cs(vx2(j),vx20(j),tmpvx(j));
            [vy2(j),vy20(j)]=add_cs(vy2(j),vy20(j),tmpvy(j));
            [dvx0(j),dvx00(j)]=add_cs(dvx0(j),dvx00(j),tmpvx(j));
            [dvy0(j),dvy00(j)]=add_cs(dvy0(j),dvy00(j),tmpvy(j));
        end
        g=grav(mu,x2,y2);
        dfx0(:,1)=g.*x2;
        dfy0(:,1)=g.*y2;
        for i=2:n+1
            dfx0(:,i)=dfx0(:,i-1)-dfx(:,i-1);
            dfy0(:,i)=dfy0(:,i-1)-dfy(:,i-1);
        end
        dfx=dfx0;
        dfy=dfy0;
    end
%     (x0-x2)
%     (vx0-vx2)
    x(t1,t2,a)=x2;
    y(t1,t2,a)=y2;
    vx(t1,t2,a)=vx2;
    vy(t1,t2,a)=vy2; 
    xcs(t1,t2,a)=x20;
    ycs(t1,t2,a)=y20;
    vxcs(t1,t2,a)=vx20;
    vycs(t1,t2,a)=vy20; 
    ex(t1,t2,a)=(x0-x2)-x20;
    ey(t1,t2,a)=(y0-y2)-y20;
    dex(t1,t2,a)=(dx1-dx2)-dx20;
    dey(t1,t2,a)=(dy1-dy2)-dy20;
    ddex(t1,t2,a)=dx0;
    ddey(t1,t2,a)=dy0;
    evx(t1,t2,a)=(vx0-vx2)-vx20;
    evy(t1,t2,a)=(vy0-vy2)-vy20;
    devx(t1,t2,a)=dvx0;
    devy(t1,t2,a)=dvy0;
    end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%5
function out=coef1(m)
syms s
out1=1;
out2=1;
    for i=1:m
        out1=out1*(s+1-i)/(m-i+1);
        out2=out2*(-s+1-i)/(m-i+1);
    end
    out=int((-1)^m*(1-s)*(out1+out2),s,0,1);
end
function out=coef2(m)
syms s
out1=1;
out2=1;
    for i=1:m
        out1=out1*(s+2-i)/(m-i+1);
        out2=out2*(-s+2-i)/(m-i+1);
    end
    out=int((-1)^m*(1-s)*(out1+out2),s,0,1);
end
function out=coef3(m)
syms s
out2=1;
    for i=1:m
        out2=out2*(-s+1-i)/(m-i+1);
    end
    out=int((-1)^m*(out2),s,0,1);
end
function out=coef4(m)
syms s
out1=1;
out2=1;
    for i=1:m
        out2=out2*(-s+2-i)/(m-i+1);
    end
    out=int((-1)^m*(out2),s,0,1);
end

function [a,csa]=add_cs(a,csa,b)
        tmp=b+csa;
        if(b==tmp)
            if(a*b<0)
                tmp=a+b;
                csa=csa+((b-tmp)+a);
                a=tmp+csa;
                csa=csa-(a-tmp);
            else
                tmp=tmp+a;
                csa=csa+(a-(tmp-b));
                a=tmp;
            end
        else
            csa=csa-(tmp-b);
            b=tmp;
            tmp=tmp+a;
            csa=csa+(b-(tmp-a));
            a=tmp;
        end
end

function out =grav(mu,x,y)
    out=1./(x.*x+y.*y);
    out=out.*sqrt(out);
    out=-out.*mu;
end
