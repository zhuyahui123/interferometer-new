clear 
clc
close all;
lamda=0.6328*10^-3;
x=linspace(-1,1,1000);
% x=-10:0.1:10;
[vv nn]=size(x); 
[x,y]=meshgrid(x,x);
% % %  z=0.0000001*peaks(500);
% c=1/100;%顶点曲率半径
% 
% %asphere 抛物面
% k=-3;
% x1=(c*(x.^2+y.^2))./(1+(1-(1+k)*c^2*(x.^2+y.^2)).^(1/2));
% % x1=x1-min(min(x1));
% %sphere 球面
% k=0;
% x2=(c*(x.^2+y.^2))./(1+(1-(1+k)*c^2*(x.^2+y.^2)).^(1/2));
% % x2=x2-min(min(x2));
% % xxx=2*(x2-x1);
% x2=x2-x1;
% xxx=2*x2;
% 
% xxx=xxx-min(min(xxx));
ff=zeros(nn);
fn=x.^2+y.^2;
for i=1:nn
     for j=1:nn
         if fn(i,j)<=(0.9*1)^2
            ff(i,j)=1;
         end
     end
end
ffx=zeros(nn);
fn=x.^2+y.^2;
for i=1:nn
     for j=1:nn
         if fn(i,j)<=(0.9*1)^2
            ffx(i,j)=1;
         end
     end
end
figure,imshow(ffx,[]);
%向上平移30Pixel
dim=1; 
k=5;
idy=repmat({':'}, ndims(ffx), 1); % initialize subscripts
nn=size(ffx,dim);   % length along dimension dim                      
idy{dim}=[k+1:nn 1:k];%向上平移5个元素
ffx=ffx(idy{:});
%向右平移30Pixel
dim=2; 
k=31;
idx=repmat({':'}, ndims(ffx), 1); % initialize subscripts
nn=size(ffx,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向you平移40个元素
ffx=ffx(idx{:});
dim=2; 
k=98;
idx=repmat({':'}, ndims(ffx), 1); % initialize subscripts
nn=size(ffx,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向zuo平移40个元素
ffx1=ffx(idx{:});
figure,imshow(ffx,[]);
figure,imshow(ffx1,[]);
Hx=ffx&ffx1;
figure,imshow(Hx,[]);

ffy=zeros(nn);
fn=x.^2+y.^2;
for i=1:nn
     for j=1:nn
         if fn(i,j)<=(0.9*1)^2
            ffy(i,j)=1;
         end
     end
end
figure,imshow(ffy,[]);
%向右平移30Pixel
dim=2; 
k=8;
idx=repmat({':'}, ndims(ffy), 1); % initialize subscripts
nn=size(ffy,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移8个元素
ffy=ffy(idx{:});
% figure,imshow(ffy,[])

%向上平移30Pixel
dim=1; 
k=31;
idy=repmat({':'}, ndims(ffy), 1); % initialize subscripts
nn=size(ffy,dim);   % length along dimension dim                      
idy{dim}=[k+1:nn 1:k];%向上平移30个元素
ffy=ffy(idy{:});
dim=1; 
k=98;
idy=repmat({':'}, ndims(ffy), 1); % initialize subscripts
nn=size(ffy,dim);   % length along dimension dim                      
idy{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
ffy1=ffy(idy{:});
figure,imshow(ffy,[]);
figure,imshow(ffy1,[]);
Hy=ffy&ffy1;
figure,imshow(Hy,[]);
% shear=nn*0.8*0.1;
shear=98;
ff1=zeros(nn);
ff1(:,1:nn-shear)=ff(:,shear+1:nn);
ffx=ff & ff1;
ffx1=ff | ff1;
% ffxx=ffx1-ffx;
figure,imshow(ffx,[]);
% %最大偏离量
% pv1=xxx.*ff;
% PV_0=(max(max(pv1))-min(min(pv1)))/lamda;
% RMS_0=max(max(std(pv1)))/lamda;
% x1=xxx.*ff;
% %xian shi bo cha
% figure(9),mesh(x,y,x1/(2*lamda)),xlabel('x(mm)');ylabel('y(mm)');zlabel('z(lamda）');
% title('原始波前图');
% 
% %direct x shearing
% xs=zeros(nn);
% xs(:,1:nn-shear)=x1(:,shear+1:nn);
% xi=(x1-xs)*2*pi/lamda;
% xs1=cos(xi).*ffx+ffx1*1;
% % figure(1),imshow(xs1,[])
% xi=(x1-xs)*2*pi/lamda+pi/2;
% xs2=cos(xi).*ffx+ffx1*1;
% % figure(2),imshow(xs2,[])
% xi=(x1-xs)*2*pi/lamda+pi;
% xs3=cos(xi).*ffx+ffx1*1;
% % figure(3),imshow(xs3,[])
% xi=(x1-xs)*2*pi/lamda+3*pi/2;
% xs4=cos(xi).*ffx+ffx1*1;
% xsgraph=[xs1 xs2;xs3 xs4];
% figure,imshow(xsgraph,[]);
% % figure(4),imshow(xs4,[])
% % %direct y shearing
ff2=zeros(nn);
ff2(1:nn-shear,:)=ff(shear+1:nn,:);
ffy=ff & ff2;
ffy1=ff | ff2;
figure,imshow(ffy,[]);
% % ffyy=ffy1-ffy;
% ys=zeros(nn);
% ys(1:nn-shear,:)=x1(shear+1:nn,:);
% yi=(x1-ys)*2*pi/lamda;
% ys1=cos(yi).*ffy+ffy1*1;
% % figure(5),imshow(ys1,[])
% yi=(x1-ys)*2*pi/lamda+pi/2;
% ys2=cos(yi).*ffy+ffy1*1;
% % figure(6),imshow(ys2,[])
% yi=(x1-ys)*2*pi/lamda+pi;
% ys3=cos(yi).*ffy+ffy1*1;
% % figure(7),imshow(ys3,[])
% yi=(x1-ys)*2*pi/lamda+3*pi/2;
% ys4=cos(yi).*ffy+ffy1*1;
% % figure(8),imshow(ys4,[])
% ysgraph=[ys1 ys2;ys3 ys4];
% figure,imshow(ysgraph,[])
% X1=imread('E:\2021-one\仿真非球面曲面\617\4\7-71-1-0.bmp');
% X2=imread('E:\2021-one\仿真非球面曲面\617\4\7-72-2-45.bmp');
% X3=imread('E:\2021-one\仿真非球面曲面\617\4\7-73-3-90.bmp');
% X4=imread('E:\2021-one\仿真非球面曲面\617\4\7-74-4-135.bmp');
% Y1=imread('E:\2021-one\仿真非球面曲面\617\4\8-81-1-0.bmp');
% Y2=imread('E:\2021-one\仿真非球面曲面\617\4\8-82-2-45.bmp');
% Y3=imread('E:\2021-one\仿真非球面曲面\617\4\8-83-3-90.bmp');
% Y4=imread('E:\2021-one\仿真非球面曲面\617\4\8-84-4-135.bmp');
X1=imread('E:\11.实验数据\20220413\4.13chuli_data\1X0.bmp');
% X1=imcrop(X1,[40 40 473 473]);
X1=double(X1(:,:,1));
% figure,imshow(X1,[])
X1=X1.*Hx;
% figure,imshow(X1,[])
X2=imread('E:\11.实验数据\20220413\4.13chuli_data\1X45.bmp');
% X2=imcrop(X2,[40 40 473 473]);
X2=double(X2(:,:,1));
% figure,imshow(X2,[])
X2=X2.*Hx;
% figure,imshow(X2,[])
X3=imread('E:\11.实验数据\20220413\4.13chuli_data\1X90.bmp');
% X3=imcrop(X3,[40 40 473 473]);
X3=double(X3(:,:,1));
% figure,imshow(X3,[])
X3=X3.*Hx;
% figure,imshow(X3,[])
X4=imread('E:\11.实验数据\20220413\4.13chuli_data\1X135.bmp');
% X4=imcrop(X4,[40 40 473 473]);
X4=double(X4(:,:,1));
% figure,imshow(X4,[])
X4=X4.*Hx;
% figure,imshow(X4,[])
Y1=imread('E:\11.实验数据\20220413\4.13chuli_data\1Y0.bmp');
% Y1=imcrop(Y1,[65 40 473 473]);
Y1=double(Y1(:,:,1));
% figure,imshow(Y1,[])
Y1=Y1.*Hy;
% figure,imshow(Y1,[])
Y2=imread('E:\11.实验数据\20220413\4.13chuli_data\1Y45.bmp');
% Y2=imcrop(Y2,[65 40 473 473]);
Y2=double(Y2(:,:,1));
% figure,imshow(Y2,[])
Y2=Y2.*Hy;
% figure,imshow(Y2,[])
Y3=imread('E:\11.实验数据\20220413\4.13chuli_data\1Y90.bmp');
% Y3=imcrop(Y3,[65 40 473 473]);
Y3=double(Y3(:,:,1));
% figure,imshow(Y3,[])
Y3=Y3.*Hy;
% figure,imshow(Y3,[])
Y4=imread('E:\11.实验数据\20220413\4.13chuli_data\1Y135.bmp');
% Y4=imcrop(Y4,[65 40 473 473]);
Y4=double(Y4(:,:,1));
% figure,imshow(Y4,[])
Y4=Y4.*Hy;
% figure,imshow(Y4,[])
Z1=X1(:,:,1);
Z2=X2(:,:,1);
Z3=X3(:,:,1);
Z4=X4(:,:,1);
Z5=Y1(:,:,1);
Z6=Y2(:,:,1);
Z7=Y3(:,:,1);
Z8=Y4(:,:,1);
I1=double(Z1);
I2=double(Z2);
I3=double(Z3);
I4=double(Z4);
I1=medfilt2(I1,[12,12]);%中值滤波
I2=medfilt2(I2,[12,12]);
I3=medfilt2(I3,[12,12]);
I4=medfilt2(I4,[12,12]);

%解包
% I1=xs1;
% I2=xs2;
% I3=xs3;
% I4=xs4;
clear i
[width,height]=size(I1);
% phasex=atan2((I4-I2),(I3-I1));
for i=1:width
    for j=1:height
        if (I1(i,j)==I3(i,j))
            if(I2(i,j)>I4(i,j))
                phasex(i,j)=pi/2;
            elseif(I2(i,j)<I4(i,j))
                phasex(i,j)=-pi/2;
            else
                phasex(i,j)=0;
            end
        elseif (I1(i,j)>I3(i,j))
            if I2(i,j)>=I4(i,j)
                phasex(i,j)=atan((I2(i,j)-I4(i,j))/(I1(i,j)-I3(i,j)));
            else
                phasex(i,j)=-atan((I4(i,j)-I2(i,j))/(I1(i,j)-I3(i,j))); 
            end
        else
            if I2(i,j)>=I4(i,j)
                phasex(i,j)=pi-atan((I2(i,j)-I4(i,j))/(I3(i,j)-I1(i,j)));
            else
                phasex(i,j)=-pi+atan((I4(i,j)-I2(i,j))/(I3(i,j)-I1(i,j))); 
            end
        end
    end
end
figure,imshow(phasex.*Hx,[]);
phasex(~isfinite(phasex))=0;
PV_px0 = (max(max(phasex))-min(min(phasex)));
RMS_px0 =(max(max(std(phasex))));

%消除2pi突变
% phase=a;
[row,col]=size(phasex);
midr=fix(row/2);
midl=fix(col/2);
for n=1:col
for i=midr-1:-1:1
    a=phasex(i,n)-phasex(i+1,n);
    if a>=pi
        for j=i:-1:1
            phasex(j,n)=phasex(j,n)-2*pi;
        end
    elseif a<=-pi
        for j=i:-1:1
            phasex(j,n)=phasex(j,n)+2*pi;
        end
    end
end
for i=midr+1:row
    a=phasex(i,n)-phasex(i-1,n);
    if a>=pi
        for j=i:row
            phasex(j,n)=phasex(j,n)-2*pi;
        end
    elseif a<=-pi
        for j=i:row
            phasex(j,n)=phasex(j,n)+2*pi;
        end
    end
end
end
for n=1:row
for i=midl-1:-1:1
    a=phasex(n,i)-phasex(n,i+1);
    if a>=pi
        for j=i:-1:1
            phasex(n,j)=phasex(n,j)-2*pi;
        end
    elseif a<=-pi
        for j=i:-1:1
            phasex(n,j)=phasex(n,j)+2*pi;
        end
    end
end
for i=midl+1:col
    a=phasex(n,i)-phasex(n,i-1);
    if a>=pi
        for j=i:col
            phasex(n,j)=phasex(n,j)-2*pi;
        end
    elseif a<=-pi
        for j=i:col
            phasex(n,j)=phasex(n,j)+2*pi;
        end
    end
end
end

% phasex=unwrap(phasex);
phasex=phasex/(2*pi);
phasex1=phasex.*Hx; 
figure,imshow(phasex1,[])
figure,mesh(phasex1)


I5=double(Z5);
I6=double(Z6);
I7=double(Z7);
I8=double(Z8);
I5=medfilt2(I5,[12,12]);%中值滤波
I6=medfilt2(I6,[12,12]);
I7=medfilt2(I7,[12,12]);
I8=medfilt2(I8,[12,12]);
% I5=ys1;
% I6=ys2;
% I7=ys3;
% I8=ys4;

[width,height]=size(I5);
%解包
% phasey=atan2((I8-I6),(I7-I5));
for i=1:width
    for j=1:height
        if (I5(i,j)==I7(i,j))
            if(I6(i,j)>I8(i,j))
                phasey(i,j)=pi/2;
            elseif(I6(i,j)<I8(i,j))
                phasey(i,j)=-pi/2;
            else
                phasey(i,j)=0;
            end
        elseif (I5(i,j)>I7(i,j))
            if I6(i,j)>=I8(i,j)
                phasey(i,j)=atan((I6(i,j)-I8(i,j))/(I5(i,j)-I7(i,j)));
            else
                phasey(i,j)=-atan((I8(i,j)-I6(i,j))/(I5(i,j)-I7(i,j))); 
            end
        else
            if I6(i,j)>=I8(i,j)
                phasey(i,j)=pi-atan((I6(i,j)-I8(i,j))/(I7(i,j)-I5(i,j)));
            else
                phasey(i,j)=-pi+atan((I8(i,j)-I6(i,j))/(I7(i,j)-I5(i,j))); 
            end
        end
    end
end
figure,imshow(phasey.*Hy,[]);
phasey(~isfinite(phasey))=0;
PV_py0 = (max(max(phasey))-min(min(phasey)));
RMS_py0 =(max(max(std(phasey))));

%消除2pi突变
% phase=a;
[row,col]=size(phasey);
midr=fix(row/2);
midl=fix(col/2);
for n=1:col
for i=midr-1:-1:1
    a=phasey(i,n)-phasey(i+1,n);
    if a>=pi
        for j=i:-1:1
            phasey(j,n)=phasey(j,n)-2*pi;
        end
    elseif a<=-pi
        for j=i:-1:1
            phasey(j,n)=phasey(j,n)+2*pi;
        end
    end
end
for i=midr+1:row
    a=phasey(i,n)-phasey(i-1,n);
    if a>=pi
        for j=i:row
            phasey(j,n)=phasey(j,n)-2*pi;
        end
    elseif a<=-pi
        for j=i:row
            phasey(j,n)=phasey(j,n)+2*pi;
        end
    end
end
end
for n=1:row
for i=midl-1:-1:1
    a=phasey(n,i)-phasey(n,i+1);
    if a>=pi
        for j=i:-1:1
            phasey(n,j)=phasey(n,j)-2*pi;
        end
    elseif a<=-pi
        for j=i:-1:1
            phasey(n,j)=phasey(n,j)+2*pi;
        end
    end
end
for i=midl+1:col
    a=phasey(n,i)-phasey(n,i-1);
    if a>=pi
        for j=i:col
            phasey(n,j)=phasey(n,j)-2*pi;
        end
    elseif a<=-pi
        for j=i:col
            phasey(n,j)=phasey(n,j)+2*pi;
        end
    end
end
end
phasey=phasey/(2*pi);
phasey1=phasey.*Hy;
figure,imshow(phasey1,[]);
figure,mesh(phasey1)


wx=phasex1*lamda/2;
wy=-phasey1*lamda/2;
% Sx=wx(:)/(2*shear);
% Sy=wy(:)/(2*shear);
Sx=wx(:);
Sy=wy(:);%按列排一列,若按行排先加‘求转置
% S=Sy;
S=[Sx;Sy];
%% 差分zernike
% x=linspace(-1,1,1000);
% [x,y]=meshgrid(x,x);
% x=linspace(-10,10,1000);
% [x,y]=meshgrid(x,x);
% [theta,rho] = cart2pol(x,y);%将笛卡尔坐标转换为极坐标。
% index_rr = find(rho<=0.8);
% M =nan(n,n);
% M(index_rr)=1;
Z1=ones(nn);
Z2=x;
Z3=y;
Z4=2*x.^2+2*y.^2-1;
Z5=x.^2-y.^2;
Z6=2.*x.*y;
Z7=3.*x.^3+3.*x.*y.^2-2.*x;
Z8=3.*x.^2.*y+3.*y.^3-2.*y;
Z9=6.*x.^4+12.*x.^2.*y.^2+6.*y.^4-6.*x.^2-6.*y.^2+1;
Z10=x.^3-3.*x.*y.^2;
Z11=3.*x.^2.*y-y.^3;
Z12=4.*x.^4-4.*y.^4-3.*x.^2+3.*y.^2;
Z13=8.*x.^3.*y+8.*x.*y.^3-6.*x.*y;
Z14=10.*x.*y.^4+20.*x.^3.*y.^2+10.*x.^5-12.*x.*y.^2-12.*x.^3+3.*x;
Z15=10.*x.^4.*y+20.*x.^2.*y.^3+10.*y.^5-12.*x.^2.*y-12.*y.^3+3.*y;
Z16=20.*x.^6+20.*y.^6+60.*x.^4.*y.^2+60.*x.^2.*y.^4-60.*x.^2.*y.^2-30.*x.^4-30.*y.^4+12.*x.^2+12.*y.^2-1;
Z17=x.^4-6.*x.^2.*y.^2+y.^4;
Z18=4.*x.^3.*y-4.*x.*y.^3;
Z19=5.*x.^5-10.*x.^3.*y.^2-15.*x.*y.^4-4.*x.^3+12.*x.*y.^2;
Z20=-5.*y.^5+10.*x.^2.*y.^3+15.*x.^4.*y+4.*y.^3-12.*x.^2.*y;
Z21=15.*x.^6-15.*y.^6+15.*x.^4.*y.^2-15.*x.^2.*y.^4-20.*x.^4+20.*y.^4+6.*x.^2-6.*y.^2;
Z22=30.*x.^5.*y+60.*x.^3.*y.^3+30.*x.*y.^5-40.*x.^3.*y-40.*x.*y.^3+12.*x.*y;
Z23=35.*x.^7+105.*x.^5.*y.^2+105.*x.^3.*y.^4+35.*x.*y.^6-60.*x.^5-120.*x.^3.*y.^2-60.*x.*y.^4+30.*x.^3+30.*x.*y.^2-4.*x;
Z24=35.*x.^6.*y+105.*x.^4.*y.^3+105.*x.^2.*y.^5+35.*y.^7-60.*x.^4.*y-120.*x.^2.*y.^3-60.*y.^5+30.*x.^2.*y+30.*y.^3-4.*y;
Z25=70.*x.^8+280.*x.^6.*y.^2+420.*x.^4.*y.^4+280.*x.^2.*y.^6+70.*y.^8-140.*x.^6-420.*x.^4.*y.^2-420.*x.^2.*y.^4-140.*y.^6+90.*x.^4+180.*x.^2.*y.^2+90.*y.^4-20.*x.^2-20.*y.^2+1;
Z26=x.^5-10.*x.^3.*y.^2+5.*x.*y.^4;
Z27=y.^5-10.*x.^2.*y.^3+5.*x.^4.*y;
Z28=6.*x.^6-30.*x.^4.*y.^2-30.*x.^2.*y.^4+30.*x.^2.*y.^2-5.*x.^4-5.*y.^4+6.*y.^6;
Z29=24.*x.^5.*y-24.*x.*y.^5-20.*x.^3.*y+20.*x.*y.^3;
Z30=21.*x.^7-30.*x.^5+10.*x.^3-21.*x.^5.*y.^2-105.*x.^3.*y.^4-63.*x.*y.^6+60.*x.^3.*y.^2+90.*x.*y.^4-30.*x.*y.^2;
Z31=63.*x.^6.*y+105.*x.^4.*y.^3+21.*x.^2.*y.^5-90.*x.^4.*y-60.*x.^2.*y.^3+30.*x.^2.*y-21.*y.^7+30.*y.^5-10.*y.^3;
Z32=56.*x.^8+112.*x.^6.*y.^2-112.*x.^2.*y.^6-105.*x.^6-105.*x.^4.*y.^2+105.*x.^2.*y.^4+60.*x.^4-10.*x.^2-56.*y.^8+105.*y.^6-60.*y.^4+10.*y.^2;
Z33=112.*x.^7.*y+336.*x.^5.*y.^3+336.*x.^3.*y.^5+112.*x.*y.^7-210.*x.^5.*y-420.*x.^3.*y.^3-210.*x.*y.^5+120.*x.^3.*y+120.*x.*y.^3-20.*x.*y;
Z34=126.*x.^9+504.*x.^7.*y.^2+756.*x.^5.*y.^4+504.*x.^3.*y.^6+126.*x.*y.^8-280.*x.^7-840.*x.^5.*y.^2-840.*x.^3.*y.^4+210.*x.^5+420.*x.^3.*y.^2+210.*x.*y.^4-60.*x.^3-60.*x.*y.^2-280.*x.*y.^6+5.*x;
Z35=126.*x.^8.*y+504.*x.^6.*y.^3+756.*x.^4.*y.^5+504.*x.^2.*y.^7+126.*y.^9-280.*x.^6.*y-840.*x.^4.*y.^3-840.*x.^2.*y.^5+210.*y.^5+420.*x.^2.*y.^3+210.*x.^4.*y-60.*y.^3-60.*x.^2.*y-280.*y.^7+5.*y;
Z36=252.*x.^10+1260.*x.^8.*y.^2+2520.*x.^6.*y.^4+2520.*x.^4.*y.^6+1260.*x.^2.*y.^8+252.*y.^10-630.*x.^8-2520.*x.^6.*y.^2-3780.*x.^4.*y.^4-2520.*x.^2.*y.^6-630.*y.^8+560.*x.^6+1680.*x.^4.*y.^2+1680.*x.^2.*y.^4+560.*y.^6-210.*x.^4-420.*x.^2.*y.^2-210.*y.^4+30.*x.^2+30.*y.^2-1;
Z1=Z1.*ff;
Z2=Z2.*ff;
Z3=Z3.*ff;
Z4=Z4.*ff;
Z5=Z5.*ff;
Z6=Z6.*ff;
Z7=Z7.*ff;
Z8=Z8.*ff;
Z9=Z9.*ff;
Z10=Z10.*ff;
Z11=Z11.*ff;
Z12=Z12.*ff;
Z13=Z13.*ff;
Z14=Z14.*ff;
Z15=Z15.*ff;
Z16=Z16.*ff;
Z17=Z17.*ff;
Z18=Z18.*ff;
Z19=Z19.*ff;
Z20=Z20.*ff;
Z21=Z21.*ff;
Z22=Z22.*ff;
Z23=Z23.*ff;
Z24=Z24.*ff;
Z25=Z25.*ff;
Z26=Z26.*ff;
Z27=Z27.*ff;
Z28=Z28.*ff;
Z28=Z28.*ff;
Z29=Z29.*ff;
Z30=Z30.*ff;
Z31=Z31.*ff;
Z32=Z32.*ff;
Z33=Z33.*ff;
Z34=Z34.*ff;
Z35=Z35.*ff;
Z36=Z36.*ff;

dim=2; 
k=8;
idx=repmat({':'}, ndims(Z1), 1); % initialize subscripts
nn=size(Z1,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za1=Z1(idx{:});


idx=repmat({':'}, ndims(Z2), 1); % initialize subscripts
nn=size(Z2,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za2=Z2(idx{:});


idx=repmat({':'}, ndims(Z3), 1); % initialize subscripts
nn=size(Z3,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za3=Z3(idx{:});


idx=repmat({':'}, ndims(Z4), 1); % initialize subscripts
nn=size(Z4,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za4=Z4(idx{:});


idx=repmat({':'}, ndims(Z5), 1); % initialize subscripts
nn=size(Z5,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za5=Z5(idx{:});

idx=repmat({':'}, ndims(Z6), 1); % initialize subscripts
nn=size(Z6,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za6=Z6(idx{:});


idx=repmat({':'}, ndims(Z7), 1); % initialize subscripts
nn=size(Z7,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za7=Z7(idx{:});


idx=repmat({':'}, ndims(Z8), 1); % initialize subscripts
nn=size(Z8,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za8=Z8(idx{:});


idx=repmat({':'}, ndims(Z9), 1); % initialize subscripts
nn=size(Z9,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za9=Z9(idx{:});

idx=repmat({':'}, ndims(Z10), 1); % initialize subscripts
nn=size(Z10,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za10=Z10(idx{:});


idx=repmat({':'}, ndims(Z11), 1); % initialize subscripts
nn=size(Z11,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za11=Z11(idx{:});


idx=repmat({':'}, ndims(Z12), 1); % initialize subscripts
nn=size(Z12,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za12=Z12(idx{:});


idx=repmat({':'}, ndims(Z13), 1); % initialize subscripts
nn=size(Z13,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za13=Z13(idx{:});



idx=repmat({':'}, ndims(Z14), 1); % initialize subscripts
nn=size(Z14,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za14=Z14(idx{:});



idx=repmat({':'}, ndims(Z15), 1); % initialize subscripts
nn=size(Z15,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za15=Z15(idx{:});


idx=repmat({':'}, ndims(Z16), 1); % initialize subscripts
nn=size(Z16,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za16=Z16(idx{:});



idx=repmat({':'}, ndims(Z17), 1); % initialize subscripts
nn=size(Z17,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za17=Z17(idx{:});


idx=repmat({':'}, ndims(Z18), 1); % initialize subscripts
nn=size(Z18,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za18=Z18(idx{:});


idx=repmat({':'}, ndims(Z19), 1); % initialize subscripts
nn=size(Z19,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za19=Z19(idx{:});

idx=repmat({':'}, ndims(Z20), 1); % initialize subscripts
nn=size(Z20,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za20=Z20(idx{:});

idx=repmat({':'}, ndims(Z21), 1); % initialize subscripts
nn=size(Z21,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za21=Z21(idx{:});

idx=repmat({':'}, ndims(Z22), 1); % initialize subscripts
nn=size(Z22,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za22=Z22(idx{:});


idx=repmat({':'}, ndims(Z23), 1); % initialize subscripts
nn=size(Z23,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za23=Z23(idx{:});


idx=repmat({':'}, ndims(Z24), 1); % initialize subscripts
nn=size(Z24,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za24=Z24(idx{:});


idx=repmat({':'}, ndims(Z25), 1); % initialize subscripts
nn=size(Z25,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za25=Z25(idx{:});

idx=repmat({':'}, ndims(Z26), 1); % initialize subscripts
nn=size(Z26,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za26=Z26(idx{:});



idx=repmat({':'}, ndims(Z27), 1); % initialize subscripts
nn=size(Z27,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za27=Z27(idx{:});



idx=repmat({':'}, ndims(Z28), 1); % initialize subscripts
nn=size(Z28,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za28=Z28(idx{:});


idx=repmat({':'}, ndims(Z29), 1); % initialize subscripts
nn=size(Z29,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za29=Z29(idx{:});


idx=repmat({':'}, ndims(Z30), 1); % initialize subscripts
nn=size(Z30,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za30=Z30(idx{:});


idx=repmat({':'}, ndims(Z31), 1); % initialize subscripts
nn=size(Z31,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za31=Z31(idx{:});


idx=repmat({':'}, ndims(Z32), 1); % initialize subscripts
nn=size(Z32,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za32=Z32(idx{:});


idx=repmat({':'}, ndims(Z33), 1); % initialize subscripts
nn=size(Z33,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za33=Z33(idx{:});


idx=repmat({':'}, ndims(Z34), 1); % initialize subscripts
nn=size(Z34,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za34=Z34(idx{:});


idx=repmat({':'}, ndims(Z35), 1); % initialize subscripts
nn=size(Z35,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za35=Z35(idx{:});


idx=repmat({':'}, ndims(Z36), 1); % initialize subscripts
nn=size(Z36,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Za36=Z36(idx{:});
figure, imshow(Za1,[]);

%% X方向
%向上平移30Pixel
dim=1; 
k=5;
idx=repmat({':'}, ndims(Z1), 1); % initialize subscripts
nn=size(Z1,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z1=Z1(idx{:});

idx=repmat({':'}, ndims(Z2), 1); % initialize subscripts
nn=size(Z2,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z2=Z2(idx{:});

idx=repmat({':'}, ndims(Z3), 1); % initialize subscripts
nn=size(Z3,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z3=Z3(idx{:});

idx=repmat({':'}, ndims(Z4), 1); % initialize subscripts
nn=size(Z4,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z4=Z4(idx{:});

idx=repmat({':'}, ndims(Z5), 1); % initialize subscripts
nn=size(Z5,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z5=Z5(idx{:});

idx=repmat({':'}, ndims(Z6), 1); % initialize subscripts
nn=size(Z6,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z6=Z6(idx{:});

idx=repmat({':'}, ndims(Z7), 1); % initialize subscripts
nn=size(Z7,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z7=Z7(idx{:});

idx=repmat({':'}, ndims(Z8), 1); % initialize subscripts
nn=size(Z8,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z8=Z8(idx{:});

idx=repmat({':'}, ndims(Z9), 1); % initialize subscripts
nn=size(Z9,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z9=Z9(idx{:});

idx=repmat({':'}, ndims(Z10), 1); % initialize subscripts
nn=size(Z10,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z10=Z10(idx{:});

idx=repmat({':'}, ndims(Z11), 1); % initialize subscripts
nn=size(Z11,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z11=Z11(idx{:});

idx=repmat({':'}, ndims(Z12), 1); % initialize subscripts
nn=size(Z12,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z12=Z12(idx{:});

idx=repmat({':'}, ndims(Z13), 1); % initialize subscripts
nn=size(Z13,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z13=Z13(idx{:});

idx=repmat({':'}, ndims(Z14), 1); % initialize subscripts
nn=size(Z14,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z14=Z14(idx{:});

idx=repmat({':'}, ndims(Z15), 1); % initialize subscripts
nn=size(Z15,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z15=Z15(idx{:});

idx=repmat({':'}, ndims(Z16), 1); % initialize subscripts
nn=size(Z16,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z16=Z16(idx{:});

idx=repmat({':'}, ndims(Z17), 1); % initialize subscripts
nn=size(Z17,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z17=Z17(idx{:});

idx=repmat({':'}, ndims(Z18), 1); % initialize subscripts
nn=size(Z18,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z18=Z18(idx{:});

idx=repmat({':'}, ndims(Z19), 1); % initialize subscripts
nn=size(Z19,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z19=Z19(idx{:});

idx=repmat({':'}, ndims(Z20), 1); % initialize subscripts
nn=size(Z20,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z20=Z20(idx{:});

idx=repmat({':'}, ndims(Z21), 1); % initialize subscripts
nn=size(Z21,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z21=Z21(idx{:});

idx=repmat({':'}, ndims(Z22), 1); % initialize subscripts
nn=size(Z22,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z22=Z22(idx{:});

idx=repmat({':'}, ndims(Z23), 1); % initialize subscripts
nn=size(Z23,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z23=Z23(idx{:});

idx=repmat({':'}, ndims(Z24), 1); % initialize subscripts
nn=size(Z24,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z24=Z24(idx{:});

idx=repmat({':'}, ndims(Z25), 1); % initialize subscripts
nn=size(Z25,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z25=Z25(idx{:});

idx=repmat({':'}, ndims(Z26), 1); % initialize subscripts
nn=size(Z26,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z26=Z26(idx{:});

idx=repmat({':'}, ndims(Z27), 1); % initialize subscripts
nn=size(Z27,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z27=Z27(idx{:});

idx=repmat({':'}, ndims(Z28), 1); % initialize subscripts
nn=size(Z28,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z28=Z28(idx{:});

idx=repmat({':'}, ndims(Z29), 1); % initialize subscripts
nn=size(Z29,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z29=Z29(idx{:});

idx=repmat({':'}, ndims(Z30), 1); % initialize subscripts
nn=size(Z30,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z30=Z30(idx{:});

idx=repmat({':'}, ndims(Z31), 1); % initialize subscripts
nn=size(Z31,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z31=Z31(idx{:});

idx=repmat({':'}, ndims(Z32), 1); % initialize subscripts
nn=size(Z32,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z32=Z32(idx{:});

idx=repmat({':'}, ndims(Z33), 1); % initialize subscripts
nn=size(Z33,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z33=Z33(idx{:});

idx=repmat({':'}, ndims(Z34), 1); % initialize subscripts
nn=size(Z34,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z34=Z34(idx{:});

idx=repmat({':'}, ndims(Z35), 1); % initialize subscripts
nn=size(Z35,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z35=Z35(idx{:});

idx=repmat({':'}, ndims(Z36), 1); % initialize subscripts
nn=size(Z36,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移5个元素
Z36=Z36(idx{:});


%向右平移30个元素
dim=2; 
k=31;
idx=repmat({':'}, ndims(Z1), 1); % initialize subscripts
nn=size(Z1,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z1=Z1(idx{:});


idx=repmat({':'}, ndims(Z2), 1); % initialize subscripts
nn=size(Z2,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z2=Z2(idx{:});


idx=repmat({':'}, ndims(Z3), 1); % initialize subscripts
nn=size(Z3,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z3=Z3(idx{:});


idx=repmat({':'}, ndims(Z4), 1); % initialize subscripts
nn=size(Z4,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z4=Z4(idx{:});


idx=repmat({':'}, ndims(Z5), 1); % initialize subscripts
nn=size(Z5,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z5=Z5(idx{:});

idx=repmat({':'}, ndims(Z6), 1); % initialize subscripts
nn=size(Z6,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z6=Z6(idx{:});


idx=repmat({':'}, ndims(Z7), 1); % initialize subscripts
nn=size(Z7,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z7=Z7(idx{:});


idx=repmat({':'}, ndims(Z8), 1); % initialize subscripts
nn=size(Z8,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z8=Z8(idx{:});


idx=repmat({':'}, ndims(Z9), 1); % initialize subscripts
nn=size(Z9,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z9=Z9(idx{:});

idx=repmat({':'}, ndims(Z10), 1); % initialize subscripts
nn=size(Z10,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z10=Z10(idx{:});


idx=repmat({':'}, ndims(Z11), 1); % initialize subscripts
nn=size(Z11,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z11=Z11(idx{:});


idx=repmat({':'}, ndims(Z12), 1); % initialize subscripts
nn=size(Z12,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z12=Z12(idx{:});


idx=repmat({':'}, ndims(Z13), 1); % initialize subscripts
nn=size(Z13,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z13=Z13(idx{:});



idx=repmat({':'}, ndims(Z14), 1); % initialize subscripts
nn=size(Z14,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z14=Z14(idx{:});



idx=repmat({':'}, ndims(Z15), 1); % initialize subscripts
nn=size(Z15,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z15=Z15(idx{:});


idx=repmat({':'}, ndims(Z16), 1); % initialize subscripts
nn=size(Z16,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z16=Z16(idx{:});



idx=repmat({':'}, ndims(Z17), 1); % initialize subscripts
nn=size(Z17,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z17=Z17(idx{:});


idx=repmat({':'}, ndims(Z18), 1); % initialize subscripts
nn=size(Z18,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z18=Z18(idx{:});


idx=repmat({':'}, ndims(Z19), 1); % initialize subscripts
nn=size(Z19,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z19=Z19(idx{:});

idx=repmat({':'}, ndims(Z20), 1); % initialize subscripts
nn=size(Z20,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z20=Z20(idx{:});

idx=repmat({':'}, ndims(Z21), 1); % initialize subscripts
nn=size(Z21,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z21=Z21(idx{:});

idx=repmat({':'}, ndims(Z22), 1); % initialize subscripts
nn=size(Z22,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z22=Z22(idx{:});


idx=repmat({':'}, ndims(Z23), 1); % initialize subscripts
nn=size(Z23,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z23=Z23(idx{:});


idx=repmat({':'}, ndims(Z24), 1); % initialize subscripts
nn=size(Z24,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z24=Z24(idx{:});


idx=repmat({':'}, ndims(Z25), 1); % initialize subscripts
nn=size(Z25,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z25=Z25(idx{:});

idx=repmat({':'}, ndims(Z26), 1); % initialize subscripts
nn=size(Z26,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z26=Z26(idx{:});



idx=repmat({':'}, ndims(Z27), 1); % initialize subscripts
nn=size(Z27,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z27=Z27(idx{:});



idx=repmat({':'}, ndims(Z28), 1); % initialize subscripts
nn=size(Z28,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z28=Z28(idx{:});


idx=repmat({':'}, ndims(Z29), 1); % initialize subscripts
nn=size(Z29,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z29=Z29(idx{:});


idx=repmat({':'}, ndims(Z30), 1); % initialize subscripts
nn=size(Z30,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z30=Z30(idx{:});


idx=repmat({':'}, ndims(Z31), 1); % initialize subscripts
nn=size(Z31,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z31=Z31(idx{:});


idx=repmat({':'}, ndims(Z32), 1); % initialize subscripts
nn=size(Z32,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z32=Z32(idx{:});


idx=repmat({':'}, ndims(Z33), 1); % initialize subscripts
nn=size(Z33,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z33=Z33(idx{:});


idx=repmat({':'}, ndims(Z34), 1); % initialize subscripts
nn=size(Z34,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z34=Z34(idx{:});


idx=repmat({':'}, ndims(Z35), 1); % initialize subscripts
nn=size(Z35,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z35=Z35(idx{:});


idx=repmat({':'}, ndims(Z36), 1); % initialize subscripts
nn=size(Z36,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z36=Z36(idx{:});
%
dim=2; 
k=98;
idx=repmat({':'}, ndims(Z1), 1); % initialize subscripts
nn=size(Z1,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向左平移98个元素
Zs1=Z1(idx{:});


idx=repmat({':'}, ndims(Z2), 1); % initialize subscripts
nn=size(Z2,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs2=Z2(idx{:});


idx=repmat({':'}, ndims(Z3), 1); % initialize subscripts
nn=size(Z3,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs3=Z3(idx{:});


idx=repmat({':'}, ndims(Z4), 1); % initialize subscripts
nn=size(Z4,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs4=Z4(idx{:});


idx=repmat({':'}, ndims(Z5), 1); % initialize subscripts
nn=size(Z5,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Zs5=Z5(idx{:});

idx=repmat({':'}, ndims(Z6), 1); % initialize subscripts
nn=size(Z6,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs6=Z6(idx{:});


idx=repmat({':'}, ndims(Z7), 1); % initialize subscripts
nn=size(Z7,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs7=Z7(idx{:});


idx=repmat({':'}, ndims(Z8), 1); % initialize subscripts
nn=size(Z8,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs8=Z8(idx{:});


idx=repmat({':'}, ndims(Z9), 1); % initialize subscripts
nn=size(Z9,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs9=Z9(idx{:});

idx=repmat({':'}, ndims(Z10), 1); % initialize subscripts
nn=size(Z10,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs10=Z10(idx{:});


idx=repmat({':'}, ndims(Z11), 1); % initialize subscripts
nn=size(Z11,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs11=Z11(idx{:});


idx=repmat({':'}, ndims(Z12), 1); % initialize subscripts
nn=size(Z12,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs12=Z12(idx{:});


idx=repmat({':'}, ndims(Z13), 1); % initialize subscripts
nn=size(Z13,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs13=Z13(idx{:});



idx=repmat({':'}, ndims(Z14), 1); % initialize subscripts
nn=size(Z14,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs14=Z14(idx{:});



idx=repmat({':'}, ndims(Z15), 1); % initialize subscripts
nn=size(Z15,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs15=Z15(idx{:});


idx=repmat({':'}, ndims(Z16), 1); % initialize subscripts
nn=size(Z16,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs16=Z16(idx{:});



idx=repmat({':'}, ndims(Z17), 1); % initialize subscripts
nn=size(Z17,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs17=Z17(idx{:});


idx=repmat({':'}, ndims(Z18), 1); % initialize subscripts
nn=size(Z18,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs18=Z18(idx{:});


idx=repmat({':'}, ndims(Z19), 1); % initialize subscripts
nn=size(Z19,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs19=Z19(idx{:});

idx=repmat({':'}, ndims(Z20), 1); % initialize subscripts
nn=size(Z20,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs20=Z20(idx{:});

idx=repmat({':'}, ndims(Z21), 1); % initialize subscripts
nn=size(Z21,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs21=Z21(idx{:});

idx=repmat({':'}, ndims(Z22), 1); % initialize subscripts
nn=size(Z22,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs22=Z22(idx{:});


idx=repmat({':'}, ndims(Z23), 1); % initialize subscripts
nn=size(Z23,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs23=Z23(idx{:});


idx=repmat({':'}, ndims(Z24), 1); % initialize subscripts
nn=size(Z24,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs24=Z24(idx{:});


idx=repmat({':'}, ndims(Z25), 1); % initialize subscripts
nn=size(Z25,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs25=Z25(idx{:});

idx=repmat({':'}, ndims(Z26), 1); % initialize subscripts
nn=size(Z26,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs26=Z26(idx{:});



idx=repmat({':'}, ndims(Z27), 1); % initialize subscripts
nn=size(Z27,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Zs27=Z27(idx{:});



idx=repmat({':'}, ndims(Z28), 1); % initialize subscripts
nn=size(Z28,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs28=Z28(idx{:});


idx=repmat({':'}, ndims(Z29), 1); % initialize subscripts
nn=size(Z29,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs29=Z29(idx{:});


idx=repmat({':'}, ndims(Z30), 1); % initialize subscripts
nn=size(Z30,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs30=Z30(idx{:});


idx=repmat({':'}, ndims(Z31), 1); % initialize subscripts
nn=size(Z31,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs31=Z31(idx{:});


idx=repmat({':'}, ndims(Z32), 1); % initialize subscripts
nn=size(Z32,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs32=Z32(idx{:});


idx=repmat({':'}, ndims(Z33), 1); % initialize subscripts
nn=size(Z33,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs33=Z33(idx{:});


idx=repmat({':'}, ndims(Z34), 1); % initialize subscripts
nn=size(Z34,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs34=Z34(idx{:});


idx=repmat({':'}, ndims(Z35), 1); % initialize subscripts
nn=size(Z35,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs35=Z35(idx{:});


idx=repmat({':'}, ndims(Z36), 1); % initialize subscripts
nn=size(Z36,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向右平移40个元素
Zs36=Z36(idx{:});
%
Hx=Z1&Zs1;
% Hx=ff;
Z1x=(Z1-Zs1).*Hx;
Z2x=(Z2-Zs2).*Hx;
Z3x=(Z3-Zs3).*Hx;
Z4x=(Z4-Zs4).*Hx;
Z5x=(Z5-Zs5).*Hx;
Z6x=(Z6-Zs6).*Hx;
Z7x=(Z7-Zs7).*Hx;
Z8x=(Z8-Zs8).*Hx;
Z9x=(Z9-Zs9).*Hx;
Z10x=(Z10-Zs10).*Hx;
Z11x=(Z11-Zs11).*Hx;
Z12x=(Z12-Zs12).*Hx;
Z13x=(Z13-Zs13).*Hx;
Z14x=(Z14-Zs14).*Hx;
Z15x=(Z15-Zs15).*Hx;
Z16x=(Z16-Zs16).*Hx;
Z17x=(Z17-Zs17).*Hx;
Z18x=(Z18-Zs18).*Hx;
Z19x=(Z19-Zs19).*Hx;
Z20x=(Z20-Zs20).*Hx;
Z21x=(Z21-Zs21).*Hx;
Z22x=(Z22-Zs22).*Hx;
Z23x=(Z23-Zs23).*Hx;
Z24x=(Z24-Zs24).*Hx;
Z25x=(Z25-Zs25).*Hx;
Z26x=(Z26-Zs26).*Hx;
Z27x=(Z27-Zs27).*Hx;
Z28x=(Z28-Zs28).*Hx;
Z29x=(Z29-Zs29).*Hx;
Z30x=(Z30-Zs30).*Hx;
Z31x=(Z31-Zs31).*Hx;
Z32x=(Z32-Zs32).*Hx;
Z33x=(Z33-Zs33).*Hx;  
Z34x=(Z34-Zs34).*Hx;
Z35x=(Z35-Zs35).*Hx;
Z36x=(Z36-Zs36).*Hx;

%% Y方向
dim=1; 
k=31;
idx=repmat({':'}, ndims(Za1), 1); % initialize subscripts
nn=size(Za1,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z1=Za1(idx{:});


idx=repmat({':'}, ndims(Za2), 1); % initialize subscripts
nn=size(Za2,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z2=Za2(idx{:});


idx=repmat({':'}, ndims(Za3), 1); % initialize subscripts
nn=size(Za3,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z3=Za3(idx{:});


idx=repmat({':'}, ndims(Za4), 1); % initialize subscripts
nn=size(Za4,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z4=Za4(idx{:});


idx=repmat({':'}, ndims(Za5), 1); % initialize subscripts
nn=size(Za5,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z5=Za5(idx{:});

idx=repmat({':'}, ndims(Za6), 1); % initialize subscripts
nn=size(Za6,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z6=Za6(idx{:});


idx=repmat({':'}, ndims(Za7), 1); % initialize subscripts
nn=size(Za7,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z7=Za7(idx{:});


idx=repmat({':'}, ndims(Za8), 1); % initialize subscripts
nn=size(Za8,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z8=Za8(idx{:});


idx=repmat({':'}, ndims(Za9), 1); % initialize subscripts
nn=size(Za9,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z9=Za9(idx{:});

idx=repmat({':'}, ndims(Za10), 1); % initialize subscripts
nn=size(Za10,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z10=Za10(idx{:});


idx=repmat({':'}, ndims(Za11), 1); % initialize subscripts
nn=size(Za11,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z11=Za11(idx{:});


idx=repmat({':'}, ndims(Za12), 1); % initialize subscripts
nn=size(Za12,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z12=Za12(idx{:});


idx=repmat({':'}, ndims(Za13), 1); % initialize subscripts
nn=size(Za13,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z13=Za13(idx{:});



idx=repmat({':'}, ndims(Za14), 1); % initialize subscripts
nn=size(Za14,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z14=Za14(idx{:});



idx=repmat({':'}, ndims(Za15), 1); % initialize subscripts
nn=size(Za15,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z15=Za15(idx{:});


idx=repmat({':'}, ndims(Za16), 1); % initialize subscripts
nn=size(Za16,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z16=Za16(idx{:});



idx=repmat({':'}, ndims(Za17), 1); % initialize subscripts
nn=size(Za17,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z17=Za17(idx{:});


idx=repmat({':'}, ndims(Za18), 1); % initialize subscripts
nn=size(Za18,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z18=Za18(idx{:});


idx=repmat({':'}, ndims(Za19), 1); % initialize subscripts
nn=size(Za19,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z19=Za19(idx{:});

idx=repmat({':'}, ndims(Za20), 1); % initialize subscripts
nn=size(Za20,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z20=Za20(idx{:});

idx=repmat({':'}, ndims(Za21), 1); % initialize subscripts
nn=size(Za21,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z21=Za21(idx{:});

idx=repmat({':'}, ndims(Za22), 1); % initialize subscripts
nn=size(Za22,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z22=Za22(idx{:});


idx=repmat({':'}, ndims(Za23), 1); % initialize subscripts
nn=size(Za23,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z23=Za23(idx{:});


idx=repmat({':'}, ndims(Za24), 1); % initialize subscripts
nn=size(Za24,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z24=Za24(idx{:});


idx=repmat({':'}, ndims(Za25), 1); % initialize subscripts
nn=size(Za25,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z25=Za25(idx{:});

idx=repmat({':'}, ndims(Za26), 1); % initialize subscripts
nn=size(Za26,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z26=Za26(idx{:});



idx=repmat({':'}, ndims(Za27), 1); % initialize subscripts
nn=size(Za27,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z27=Za27(idx{:});



idx=repmat({':'}, ndims(Za28), 1); % initialize subscripts
nn=size(Za28,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z28=Za28(idx{:});


idx=repmat({':'}, ndims(Za29), 1); % initialize subscripts
nn=size(Za29,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z29=Za29(idx{:});


idx=repmat({':'}, ndims(Za30), 1); % initialize subscripts
nn=size(Za30,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z30=Za30(idx{:});


idx=repmat({':'}, ndims(Za31), 1); % initialize subscripts
nn=size(Za31,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z31=Za31(idx{:});


idx=repmat({':'}, ndims(Za32), 1); % initialize subscripts
nn=size(Za32,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z32=Za32(idx{:});


idx=repmat({':'}, ndims(Za33), 1); % initialize subscripts
nn=size(Za33,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z33=Za33(idx{:});


idx=repmat({':'}, ndims(Za34), 1); % initialize subscripts
nn=size(Za34,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z34=Za34(idx{:});


idx=repmat({':'}, ndims(Za35), 1); % initialize subscripts
nn=size(Za35,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z35=Za35(idx{:});


idx=repmat({':'}, ndims(Za36), 1); % initialize subscripts
nn=size(Za36,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z36=Za36(idx{:});
%
dim=1; 
k=98;
idx=repmat({':'}, ndims(Z1), 1); % initialize subscripts
nn=size(Z1,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis1=Z1(idx{:});


idx=repmat({':'}, ndims(Z2), 1); % initialize subscripts
nn=size(Z2,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis2=Z2(idx{:});


idx=repmat({':'}, ndims(Z3), 1); % initialize subscripts
nn=size(Z3,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis3=Z3(idx{:});


idx=repmat({':'}, ndims(Z4), 1); % initialize subscripts
nn=size(Z4,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis4=Z4(idx{:});


idx=repmat({':'}, ndims(Z5), 1); % initialize subscripts
nn=size(Z5,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis5=Z5(idx{:});

idx=repmat({':'}, ndims(Z6), 1); % initialize subscripts
nn=size(Z6,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis6=Z6(idx{:});


idx=repmat({':'}, ndims(Z7), 1); % initialize subscripts
nn=size(Z7,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis7=Z7(idx{:});


idx=repmat({':'}, ndims(Z8), 1); % initialize subscripts
nn=size(Z8,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis8=Z8(idx{:});


idx=repmat({':'}, ndims(Z9), 1); % initialize subscripts
nn=size(Z9,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis9=Z9(idx{:});

idx=repmat({':'}, ndims(Z10), 1); % initialize subscripts
nn=size(Z10,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis10=Z10(idx{:});


idx=repmat({':'}, ndims(Z11), 1); % initialize subscripts
nn=size(Z11,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis11=Z11(idx{:});


idx=repmat({':'}, ndims(Z12), 1); % initialize subscripts
nn=size(Z12,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis12=Z12(idx{:});


idx=repmat({':'}, ndims(Z13), 1); % initialize subscripts
nn=size(Z13,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis13=Z13(idx{:});


idx=repmat({':'}, ndims(Z14), 1); % initialize subscripts
nn=size(Z14,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis14=Z14(idx{:});



idx=repmat({':'}, ndims(Z15), 1); % initialize subscripts
nn=size(Z15,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis15=Z15(idx{:});


idx=repmat({':'}, ndims(Z16), 1); % initialize subscripts
nn=size(Z16,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis16=Z16(idx{:});



idx=repmat({':'}, ndims(Z17), 1); % initialize subscripts
nn=size(Z17,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis17=Z17(idx{:});


idx=repmat({':'}, ndims(Z18), 1); % initialize subscripts
nn=size(Z18,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis18=Z18(idx{:});


idx=repmat({':'}, ndims(Z19), 1); % initialize subscripts
nn=size(Z19,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis19=Z19(idx{:});

idx=repmat({':'}, ndims(Z20), 1); % initialize subscripts
nn=size(Z20,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis20=Z20(idx{:});

idx=repmat({':'}, ndims(Z21), 1); % initialize subscripts
nn=size(Z21,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis21=Z21(idx{:});

idx=repmat({':'}, ndims(Z22), 1); % initialize subscripts
nn=size(Z22,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis22=Z22(idx{:});


idx=repmat({':'}, ndims(Z23), 1); % initialize subscripts
nn=size(Z23,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis23=Z23(idx{:});


idx=repmat({':'}, ndims(Z24), 1); % initialize subscripts
nn=size(Z24,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis24=Z24(idx{:});


idx=repmat({':'}, ndims(Z25), 1); % initialize subscripts
nn=size(Z25,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis25=Z25(idx{:});

idx=repmat({':'}, ndims(Z26), 1); % initialize subscripts
nn=size(Z26,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis26=Z26(idx{:});



idx=repmat({':'}, ndims(Z27), 1); % initialize subscripts
nn=size(Z27,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis27=Z27(idx{:});



idx=repmat({':'}, ndims(Z28), 1); % initialize subscripts
nn=size(Z28,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis28=Z28(idx{:});


idx=repmat({':'}, ndims(Z29), 1); % initialize subscripts
nn=size(Z29,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis29=Z29(idx{:});


idx=repmat({':'}, ndims(Z30), 1); % initialize subscripts
nn=size(Z30,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis30=Z30(idx{:});


idx=repmat({':'}, ndims(Z31), 1); % initialize subscripts
nn=size(Z31,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis31=Z31(idx{:});


idx=repmat({':'}, ndims(Z32), 1); % initialize subscripts
nn=size(Z32,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis32=Z32(idx{:});


idx=repmat({':'}, ndims(Z33), 1); % initialize subscripts
nn=size(Z33,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis33=Z33(idx{:});


idx=repmat({':'}, ndims(Z34), 1); % initialize subscripts
nn=size(Z34,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis34=Z34(idx{:});


idx=repmat({':'}, ndims(Z35), 1); % initialize subscripts
nn=size(Z35,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis35=Z35(idx{:});


idx=repmat({':'}, ndims(Z36), 1); % initialize subscripts
nn=size(Z36,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Zis36=Z36(idx{:});
% figure,imshow(Z36,[]);
%
Hy=Z1&Zis1;
% Hy=ff;
figure,imshow(Z1,[]);
figure,imshow(Zis1,[]);
figure,imshow(Hy,[]);
Z1y=(Z1-Zis1).*Hy;
Z2y=(Z2-Zis2).*Hy;
Z3y=(Z3-Zis3).*Hy;
Z4y=(Z4-Zis4).*Hy;
Z5y=(Z5-Zis5).*Hy;
Z6y=(Z6-Zis6).*Hy;
Z7y=(Z7-Zis7).*Hy;
Z8y=(Z8-Zis8).*Hy;
Z9y=(Z9-Zis9).*Hy;
Z10y=(Z10-Zis10).*Hy;
Z11y=(Z11-Zis11).*Hy;
Z12y=(Z12-Zis12).*Hy;
Z13y=(Z13-Zis13).*Hy;
Z14y=(Z14-Zis14).*Hy;
Z15y=(Z15-Zis15).*Hy;
Z16y=(Z16-Zis16).*Hy;
Z17y=(Z17-Zis17).*Hy;
Z18y=(Z18-Zis18).*Hy;
Z19y=(Z19-Zis19).*Hy;
Z20y=(Z20-Zis20).*Hy;
Z21y=(Z21-Zis21).*Hy;
Z22y=(Z22-Zis22).*Hy;
Z23y=(Z23-Zis23).*Hy;
Z24y=(Z24-Zis24).*Hy;
Z25y=(Z25-Zis25).*Hy;
Z26y=(Z26-Zis26).*Hy;
Z27y=(Z27-Zis27).*Hy;
Z28y=(Z28-Zis28).*Hy;
Z29y=(Z29-Zis29).*Hy;
Z30y=(Z30-Zis30).*Hy;
Z31y=(Z31-Zis31).*Hy;
Z32y=(Z32-Zis32).*Hy;
Z33y=(Z33-Zis33).*Hy;  
Z34y=(Z34-Zis34).*Hy;
Z35y=(Z35-Zis35).*Hy;
Z36y=(Z36-Zis36).*Hy;

% h=1;
% [Z1x,Z1y]=gradient(Z1,h);
% [Z2x,Z2y]=gradient(Z2,h);
% [Z3x,Z3y]=gradient(Z3,h);
% [Z4x,Z4y]=gradient(Z4,h);
% [Z5x,Z5y]=gradient(Z5,h);
% [Z6x,Z6y]=gradient(Z6,h);
% [Z7x,Z7y]=gradient(Z7,h);
% [Z8x,Z8y]=gradient(Z8,h);
% [Z9x,Z9y]=gradient(Z9,h);
% [Z10x,Z10y]=gradient(Z10,h);
% [Z11x,Z11y]=gradient(Z11,h);
% [Z12x,Z12y]=gradient(Z12,h);
% [Z13x,Z13y]=gradient(Z13,h);
% [Z14x,Z14y]=gradient(Z14,h);
% [Z15x,Z15y]=gradient(Z15,h);
% [Z16x,Z16y]=gradient(Z16,h);
% [Z17x,Z17y]=gradient(Z17,h);
% [Z18x,Z18y]=gradient(Z18,h);
% [Z19x,Z19y]=gradient(Z19,h);
% [Z20x,Z20y]=gradient(Z20,h);
% [Z21x,Z21y]=gradient(Z21,h);
% [Z22x,Z22y]=gradient(Z22,h);
% [Z23x,Z23y]=gradient(Z23,h);
% [Z24x,Z24y]=gradient(Z24,h);
% [Z25x,Z25y]=gradient(Z25,h);
% [Z26x,Z26y]=gradient(Z26,h);
% [Z27x,Z27y]=gradient(Z27,h);
% [Z28x,Z28y]=gradient(Z28,h);
% [Z29x,Z29y]=gradient(Z29,h);
% [Z30x,Z30y]=gradient(Z30,h);
% [Z31x,Z31y]=gradient(Z31,h);
% [Z32x,Z32y]=gradient(Z32,h);
% [Z33x,Z33y]=gradient(Z33,h);
% [Z34x,Z34y]=gradient(Z34,h);
% [Z35x,Z35y]=gradient(Z35,h);
% [Z36x,Z36y]=gradient(Z36,h);



dim=1; 
k=31;
idx=repmat({':'}, ndims(Z1), 1); % initialize subscripts
nn=size(Z1,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z1=Z1(idx{:});


idx=repmat({':'}, ndims(Z2), 1); % initialize subscripts
nn=size(Z2,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z2=Z2(idx{:});


idx=repmat({':'}, ndims(Z3), 1); % initialize subscripts
nn=size(Z3,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z3=Z3(idx{:});


idx=repmat({':'}, ndims(Z4), 1); % initialize subscripts
nn=size(Z4,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z4=Z4(idx{:});


idx=repmat({':'}, ndims(Z5), 1); % initialize subscripts
nn=size(Z5,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z5=Z5(idx{:});

idx=repmat({':'}, ndims(Z6), 1); % initialize subscripts
nn=size(Z6,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z6=Z6(idx{:});


idx=repmat({':'}, ndims(Z7), 1); % initialize subscripts
nn=size(Z7,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z7=Z7(idx{:});


idx=repmat({':'}, ndims(Z8), 1); % initialize subscripts
nn=size(Z8,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z8=Z8(idx{:});


idx=repmat({':'}, ndims(Z9), 1); % initialize subscripts
nn=size(Z9,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z9=Z9(idx{:});

idx=repmat({':'}, ndims(Z10), 1); % initialize subscripts
nn=size(Z10,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z10=Z10(idx{:});


idx=repmat({':'}, ndims(Z11), 1); % initialize subscripts
nn=size(Z11,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z11=Z11(idx{:});


idx=repmat({':'}, ndims(Z12), 1); % initialize subscripts
nn=size(Z12,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z12=Z12(idx{:});


idx=repmat({':'}, ndims(Z13), 1); % initialize subscripts
nn=size(Z13,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z13=Z13(idx{:});


idx=repmat({':'}, ndims(Z14), 1); % initialize subscripts
nn=size(Z14,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z14=Z14(idx{:});



idx=repmat({':'}, ndims(Z15), 1); % initialize subscripts
nn=size(Z15,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z15=Z15(idx{:});


idx=repmat({':'}, ndims(Z16), 1); % initialize subscripts
nn=size(Z16,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z16=Z16(idx{:});



idx=repmat({':'}, ndims(Z17), 1); % initialize subscripts
nn=size(Z17,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z17=Z17(idx{:});


idx=repmat({':'}, ndims(Z18), 1); % initialize subscripts
nn=size(Z18,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z18=Z18(idx{:});


idx=repmat({':'}, ndims(Z19), 1); % initialize subscripts
nn=size(Z19,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z19=Z19(idx{:});

idx=repmat({':'}, ndims(Z20), 1); % initialize subscripts
nn=size(Z20,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z20=Z20(idx{:});

idx=repmat({':'}, ndims(Z21), 1); % initialize subscripts
nn=size(Z21,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z21=Z21(idx{:});

idx=repmat({':'}, ndims(Z22), 1); % initialize subscripts
nn=size(Z22,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z22=Z22(idx{:});


idx=repmat({':'}, ndims(Z23), 1); % initialize subscripts
nn=size(Z23,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z23=Z23(idx{:});


idx=repmat({':'}, ndims(Z24), 1); % initialize subscripts
nn=size(Z24,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z24=Z24(idx{:});


idx=repmat({':'}, ndims(Z25), 1); % initialize subscripts
nn=size(Z25,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z25=Z25(idx{:});

idx=repmat({':'}, ndims(Z26), 1); % initialize subscripts
nn=size(Z26,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z26=Z26(idx{:});



idx=repmat({':'}, ndims(Z27), 1); % initialize subscripts
nn=size(Z27,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z27=Z27(idx{:});



idx=repmat({':'}, ndims(Z28), 1); % initialize subscripts
nn=size(Z28,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z28=Z28(idx{:});


idx=repmat({':'}, ndims(Z29), 1); % initialize subscripts
nn=size(Z29,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z29=Z29(idx{:});


idx=repmat({':'}, ndims(Z30), 1); % initialize subscripts
nn=size(Z30,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z30=Z30(idx{:});


idx=repmat({':'}, ndims(Z31), 1); % initialize subscripts
nn=size(Z31,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z31=Z31(idx{:});


idx=repmat({':'}, ndims(Z32), 1); % initialize subscripts
nn=size(Z32,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z32=Z32(idx{:});


idx=repmat({':'}, ndims(Z33), 1); % initialize subscripts
nn=size(Z33,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z33=Z33(idx{:});


idx=repmat({':'}, ndims(Z34), 1); % initialize subscripts
nn=size(Z34,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z34=Z34(idx{:});


idx=repmat({':'}, ndims(Z35), 1); % initialize subscripts
nn=size(Z35,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z35=Z35(idx{:});


idx=repmat({':'}, ndims(Z36), 1); % initialize subscripts
nn=size(Z36,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向下平移90个元素
Z36=Z36(idx{:});
% figure,imshow(Z36,[]);

dim=2; 
k=8;
idx=repmat({':'}, ndims(Z1), 1); % initialize subscripts
nn=size(Z1,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z1=Z1(idx{:});


idx=repmat({':'}, ndims(Z2), 1); % initialize subscripts
nn=size(Z2,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z2=Z2(idx{:});


idx=repmat({':'}, ndims(Z3), 1); % initialize subscripts
nn=size(Z3,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z3=Z3(idx{:});


idx=repmat({':'}, ndims(Z4), 1); % initialize subscripts
nn=size(Z4,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z4=Z4(idx{:});


idx=repmat({':'}, ndims(Z5), 1); % initialize subscripts
nn=size(Z5,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z5=Z5(idx{:});

idx=repmat({':'}, ndims(Z6), 1); % initialize subscripts
nn=size(Z6,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z6=Z6(idx{:});


idx=repmat({':'}, ndims(Z7), 1); % initialize subscripts
nn=size(Z7,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z7=Z7(idx{:});


idx=repmat({':'}, ndims(Z8), 1); % initialize subscripts
nn=size(Z8,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z8=Z8(idx{:});


idx=repmat({':'}, ndims(Z9), 1); % initialize subscripts
nn=size(Z9,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z9=Z9(idx{:});

idx=repmat({':'}, ndims(Z10), 1); % initialize subscripts
nn=size(Z10,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z10=Z10(idx{:});


idx=repmat({':'}, ndims(Z11), 1); % initialize subscripts
nn=size(Z11,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z11=Z11(idx{:});


idx=repmat({':'}, ndims(Z12), 1); % initialize subscripts
nn=size(Z12,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z12=Z12(idx{:});


idx=repmat({':'}, ndims(Z13), 1); % initialize subscripts
nn=size(Z13,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z13=Z13(idx{:});



idx=repmat({':'}, ndims(Z14), 1); % initialize subscripts
nn=size(Z14,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z14=Z14(idx{:});



idx=repmat({':'}, ndims(Z15), 1); % initialize subscripts
nn=size(Z15,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z15=Z15(idx{:});


idx=repmat({':'}, ndims(Z16), 1); % initialize subscripts
nn=size(Z16,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z16=Z16(idx{:});



idx=repmat({':'}, ndims(Z17), 1); % initialize subscripts
nn=size(Z17,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z17=Z17(idx{:});


idx=repmat({':'}, ndims(Z18), 1); % initialize subscripts
nn=size(Z18,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z18=Z18(idx{:});


idx=repmat({':'}, ndims(Z19), 1); % initialize subscripts
nn=size(Z19,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z19=Z19(idx{:});

idx=repmat({':'}, ndims(Z20), 1); % initialize subscripts
nn=size(Z20,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z20=Z20(idx{:});

idx=repmat({':'}, ndims(Z21), 1); % initialize subscripts
nn=size(Z21,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z21=Z21(idx{:});

idx=repmat({':'}, ndims(Z22), 1); % initialize subscripts
nn=size(Z22,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z22=Z22(idx{:});


idx=repmat({':'}, ndims(Z23), 1); % initialize subscripts
nn=size(Z23,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z23=Z23(idx{:});


idx=repmat({':'}, ndims(Z24), 1); % initialize subscripts
nn=size(Z24,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z24=Z24(idx{:});


idx=repmat({':'}, ndims(Z25), 1); % initialize subscripts
nn=size(Z25,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z25=Z25(idx{:});

idx=repmat({':'}, ndims(Z26), 1); % initialize subscripts
nn=size(Z26,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z26=Z26(idx{:});



idx=repmat({':'}, ndims(Z27), 1); % initialize subscripts
nn=size(Z27,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z27=Z27(idx{:});



idx=repmat({':'}, ndims(Z28), 1); % initialize subscripts
nn=size(Z28,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z28=Z28(idx{:});


idx=repmat({':'}, ndims(Z29), 1); % initialize subscripts
nn=size(Z29,dim);   % length along dimension dim                      
idx{dim}=[nn-k+1:nn,1:nn-k];%向右平移40个元素
Z29=Z29(idx{:});


idx=repmat({':'}, ndims(Z30), 1); % initialize subscripts
nn=size(Z30,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z30=Z30(idx{:});


idx=repmat({':'}, ndims(Z31), 1); % initialize subscripts
nn=size(Z31,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z31=Z31(idx{:});


idx=repmat({':'}, ndims(Z32), 1); % initialize subscripts
nn=size(Z32,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z32=Z32(idx{:});


idx=repmat({':'}, ndims(Z33), 1); % initialize subscripts
nn=size(Z33,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z33=Z33(idx{:});


idx=repmat({':'}, ndims(Z34), 1); % initialize subscripts
nn=size(Z34,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z34=Z34(idx{:});


idx=repmat({':'}, ndims(Z35), 1); % initialize subscripts
nn=size(Z35,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z35=Z35(idx{:});


idx=repmat({':'}, ndims(Z36), 1); % initialize subscripts
nn=size(Z36,dim);   % length along dimension dim                      
idx{dim}=[k+1:nn 1:k];%向上平移30个元素
Z36=Z36(idx{:});
% figure,imshow(Z36,[]);
%%
Z1x=Z1x(:);
Z1y=Z1y(:);
Z1xy=[Z1x;Z1y];
Z2x=Z2x(:);
Z2y=Z2y(:);
Z2xy=[Z2x;Z2y];
Z3x=Z3x(:);
Z3y=Z3y(:);
Z3xy=[Z3x;Z3y];
Z4x=Z4x(:);
Z4y=Z4y(:);
Z4xy=[Z4x;Z4y];
Z5x=Z5x(:);
Z5y=Z5y(:);
Z5xy=[Z5x;Z5y];
Z6x=Z6x(:);
Z6y=Z6y(:);
Z6xy=[Z6x;Z6y];
Z7x=Z7x(:);
Z7y=Z7y(:);
Z7xy=[Z7x;Z7y];
Z8x=Z8x(:);
Z8y=Z8y(:);
Z8xy=[Z8x;Z8y];
Z9x=Z9x(:);
Z9y=Z9y(:);
Z9xy=[Z9x;Z9y];
Z10x=Z10x(:);
Z10y=Z10y(:);
Z10xy=[Z10x;Z10y];
Z11x=Z11x(:);
Z11y=Z11y(:);
Z11xy=[Z11x;Z11y];
Z12x=Z12x(:);
Z12y=Z12y(:);
Z12xy=[Z12x;Z12y];
Z13x=Z13x(:);
Z13y=Z13y(:);
Z13xy=[Z13x;Z13y];
Z14x=Z14x(:);
Z14y=Z14y(:);
Z14xy=[Z14x;Z14y];
Z15x=Z15x(:);
Z15y=Z15y(:);
Z15xy=[Z15x;Z15y];
Z16x=Z16x(:);
Z16y=Z16y(:);
Z16xy=[Z16x;Z16y];
Z17x=Z17x(:);
Z17y=Z17y(:);
Z17xy=[Z17x;Z17y];
Z18x=Z18x(:);
Z18y=Z18y(:);
Z18xy=[Z18x;Z18y];
Z19x=Z19x(:);
Z19y=Z19y(:);
Z19xy=[Z19x;Z19y];
Z20x=Z20x(:);
Z20y=Z20y(:);
Z20xy=[Z20x;Z20y];
Z21x=Z21x(:);
Z21y=Z21y(:);
Z21xy=[Z21x;Z21y];
Z22x=Z22x(:);
Z22y=Z22y(:);
Z22xy=[Z22x;Z22y];
Z23x=Z23x(:);
Z23y=Z23y(:);
Z23xy=[Z23x;Z23y];
Z24x=Z24x(:);
Z24y=Z24y(:);
Z24xy=[Z24x;Z24y];
Z25x=Z25x(:);
Z25y=Z25y(:);
Z25xy=[Z25x;Z25y];
Z26x=Z26x(:);
Z26y=Z26y(:);
Z26xy=[Z26x;Z26y];
Z27x=Z27x(:);
Z27y=Z27y(:);
Z27xy=[Z27x;Z27y];
Z28x=Z28x(:);
Z28y=Z28y(:);
Z28xy=[Z28x;Z28y];
Z29x=Z29x(:);
Z29y=Z29y(:);
Z29xy=[Z29x;Z29y];
Z30x=Z30x(:);
Z30y=Z30y(:);
Z30xy=[Z30x;Z30y];
Z31x=Z31x(:);
Z31y=Z31y(:);
Z31xy=[Z31x;Z31y];
Z32x=Z32x(:);
Z32y=Z32y(:);
Z32xy=[Z32x;Z32y];
Z33x=Z33x(:); 
Z33y=Z33y(:);
Z33xy=[Z33x;Z33y];
Z34x=Z34x(:);
Z34y=Z34y(:);
Z34xy=[Z34x;Z34y];
Z35x=Z35x(:);
Z35y=Z35y(:);
Z35xy=[Z35x;Z35y];
Z36x=Z36x(:);
Z36y=Z36y(:);
Z36xy=[Z36x;Z36y];

A=[Z1xy Z2xy Z3xy Z4xy Z5xy Z6xy Z7xy Z8xy Z9xy Z10xy Z11xy Z12xy Z13xy Z14xy Z15xy Z16xy Z17xy Z18xy Z19xy Z20xy Z21xy Z22xy Z23xy Z24xy Z25xy Z26xy Z27xy Z28xy Z29xy Z30xy Z31xy Z32xy Z33xy Z34xy Z35xy Z36xy];
% SS=S+S1;
tnd_31=find(isnan(S(:))==1);
S(tnd_31)=[];
A(tnd_31,:)=[];
A(~isfinite(A))=0;
aa=(pinv(A'*A))*(A'*S);

ff1=zeros(nn);
fn=x.^2+y.^2;
for i=1:nn
     for j=1:nn
         if fn(i,j)<=(0.9*1)^2
            ff1(i,j)=1;
         else
            ff1(i,j)=nan;
         end
     end
end
% dim=2; 
% k=8;
% idx=repmat({':'}, ndims(ff1), 1); % initialize subscripts
% nn=size(ff1,dim);   % length along dimension dim                      
% idx{dim}=[nn-k+1:nn,1:nn-k];%向左平移30个元素
% % idx{dim}=[k+1:nn 1:k];
% ff1=ff1(idx{:});
% dim=1; 
% k=36;
% idx=repmat({':'}, ndims(ff1), 1); % initialize subscripts
% nn=size(ff1,dim);   % length along dimension dim                      
% % idx{dim}=[nn-k+1:nn,1:nn-k];%向上平移30个元素
% idx{dim}=[k+1:nn 1:k];
% ff1=ff1(idx{:});

% 面形拟合aa(1).*Z1+aa(2).*Z2+aa(3).*Z3+aa(4).*Z4+
ZZ0=(aa(5).*Z5+aa(6).*Z6+aa(7).*Z7+aa(8).*Z8+aa(9).*Z9+aa(10).*Z10+aa(11).*Z11+aa(12).*Z12+aa(13).*Z13+aa(14).*Z14+aa(15).*Z15+aa(16).*Z16+aa(17).*Z17+aa(18).*Z18+aa(19).*Z19+aa(20).*Z20+aa(21).*Z21+aa(22).*Z22+aa(23).*Z23+aa(24).*Z24+aa(25).*Z25+aa(26).*Z26+aa(27).*Z27+aa(28).*Z28+aa(29).*Z29+aa(30).*Z30+aa(31).*Z31+aa(32).*Z32+aa(33).*Z33+aa(34).*Z34+aa(35).*Z35+aa(36).*Z36);
ZZ0=ZZ0-min(min(ZZ0));
ZZ0=ZZ0.*ff1;
%拟合面型
figure,mesh(x,y,ZZ0/(2*lamda)),title('拟合面形');
xlabel('x(cm)');ylabel('y(cm)');zlabel('z（mm）');

% xlabel('x(mm)');ylabel('y(mm)');zlabel('z（lamda）');
ZZ0(~isfinite(ZZ0))=0;
PV_Z = (max(max(ZZ0))-min(min(ZZ0)))/lamda;
RMS_Z =(max(max(std(ZZ0))))/lamda;

load('ZZ0dkj12');
figure,mesh(x,y,dkj);
title('Zygo面形')
% dkj(~isfinite(dkj))=0;
% PV_Zygo = (max(max(dkj))-min(min(dkj)));
% RMS_Zygo =(max(max(std(dkj))));

% privateCurve=dkj(500,:);
% ZZ0=ZZ0/(2*lamda);
% privateCurve1=ZZ0(500,:);
% figure,plot(privateCurve,'b.-')
% hold on
% plot(privateCurve1,'r.-')
% % legend('Zygo','BHPF(N=10)')
% xlabel('x(pixel)');
% ylabel('y(wavelength)');

% %残差
% ZZc=ZZ0-dkj;
% figure,mesh(x,y,ZZc),title('残差面形');
% xlabel('x(cm)');ylabel('y(cm)');zlabel('z（mm）');
% 
% ZZc(~isfinite(ZZc))=0;
% PV_c = (max(max(ZZc))-min(min(ZZc)));
% RMS_c =(max(max(std(ZZc))));

%direct x shearing
% dim=1; 
% k=31;
% idx=repmat({':'}, ndims(ffx), 1); % initialize subscripts
% nn=size(ffx,dim);   % length along dimension dim                      
% % idx{dim}=[nn-k+1:nn,1:nn-k];%向上平移30个元素
% idx{dim}=[k+1:nn 1:k];
% ffx=ffx(idx{:});
% figure,imshow(ffx,[]);
% 
% dim=2; 
% k=8;
% idx=repmat({':'}, ndims(ffx), 1); % initialize subscripts
% nn=size(ffx,dim);   % length along dimension dim                      
% idx{dim}=[nn-k+1:nn,1:nn-k];%向左平移30个元素
% % idx{dim}=[k+1:nn 1:k];
% ffx=ffx(idx{:});
% figure,imshow(ffx,[]);
% 
% dim=1; 
% k=31;
% idx=repmat({':'}, ndims(ffx1), 1); % initialize subscripts
% nn=size(ffx1,dim);   % length along dimension dim                      
% % idx{dim}=[nn-k+1:nn,1:nn-k];%向上平移30个元素
% idx{dim}=[k+1:nn 1:k];
% ffx1=ffx1(idx{:});
% 
% dim=2; 
% k=8;
% idx=repmat({':'}, ndims(ffx1), 1); % initialize subscripts
% nn=size(ffx1,dim);   % length along dimension dim                      
% idx{dim}=[nn-k+1:nn,1:nn-k];%向左平移30个元素
% % idx{dim}=[k+1:nn 1:k];
% ffx1=ffx1(idx{:});
% dkj=(0.001*Z1+0*Z2+0*Z3+0.003*Z4+0.034*Z5+(-0.018)*Z6+0.007*Z7+0.028*Z8+(-0.048)*Z9+0.009*Z10+(-0.009)*Z11+(-0.0037)*Z12+0.0018*Z13+(-0.003)*Z14+(-0.014).*Z15+(-0.080)*Z16+(-0.002)*Z17+(0.006)*Z18+(-0.004).*Z19+0.005*Z20+0.011*Z21+0.006*Z22+(-0.001)*Z23+0.003*Z24+0.015*Z25+(-0.004)*Z26+0.004*Z27+0.003*Z28+(-0.004)*Z29+0.002*Z30+0.002*Z31+0.003*Z32+(-0.004)*Z33+0.012*Z34+(-0.006)*Z35+(-0.031)*Z36);
% ZZ0=dkj;
figure,mesh(x,y,ZZ0);
xs=zeros(nn);
xs(:,1:nn-shear)=ZZ0(:,shear+1:nn);
figure,imshow(ZZ0-xs,[]);
xi=(ZZ0-xs)*0.03*pi/lamda;
xs1=cos(xi).*ffx;
figure,imshow(xs1,[])

xi=(ZZ0-xs)*4*pi/lamda+pi/2;
xs2=cos(xi).*ffx;
figure,imshow(xs2,[])

xi=(ZZ0-xs)*4*pi/lamda+pi;
xs3=cos(xi).*ffx;
figure,imshow(xs3,[])

xi=(ZZ0-xs)*4*pi/lamda+3*pi/2;
xs4=cos(xi).*ffx;
figure,imshow(xs4,[])


I1=xs1;
I2=xs2;
I3=xs3;
I4=xs4;

% phasex=atan2((I4-I2),(I3-I1));
for i=1:width
    for j=1:height
        if (I1(i,j)==I3(i,j))
            if(I2(i,j)>I4(i,j))
                phasex(i,j)=pi/2;
            elseif(I2(i,j)<I4(i,j))
                phasex(i,j)=-pi/2;
            else
                phasex(i,j)=0;
            end
        elseif (I1(i,j)>I3(i,j))
            if I2(i,j)>=I4(i,j)
                phasex(i,j)=atan((I2(i,j)-I4(i,j))/(I1(i,j)-I3(i,j)));
            else
                phasex(i,j)=-atan((I4(i,j)-I2(i,j))/(I1(i,j)-I3(i,j))); 
            end
        else
            if I2(i,j)>=I4(i,j)
                phasex(i,j)=pi-atan((I2(i,j)-I4(i,j))/(I3(i,j)-I1(i,j)));
            else
                phasex(i,j)=-pi+atan((I4(i,j)-I2(i,j))/(I3(i,j)-I1(i,j))); 
            end
        end
    end
end
figure,imshow(phasex,[]);

phasex(~isfinite(phasex))=0;
PV_px = (max(max(phasex))-min(min(phasex)));
RMS_px =(max(max(std(phasex))));


%direct x shearing
% dim=1; 
% k=31;
% idx=repmat({':'}, ndims(ffy), 1); % initialize subscripts
% nn=size(ffy,dim);   % length along dimension dim                      
% % idx{dim}=[nn-k+1:nn,1:nn-k];%向上平移30个元素
% idx{dim}=[k+1:nn 1:k];
% ffy=ffy(idx{:});
% 
% dim=2; 
% k=8;
% idx=repmat({':'}, ndims(ffy), 1); % initialize subscripts
% nn=size(ffy,dim);   % length along dimension dim                      
% idx{dim}=[nn-k+1:nn,1:nn-k];%向左平移30个元素
% % idx{dim}=[k+1:nn 1:k];
% ffy=ffy(idx{:});
% 
% dim=1; 
% k=31;
% idx=repmat({':'}, ndims(ffy1), 1); % initialize subscripts
% nn=size(ffy1,dim);   % length along dimension dim                      
% % idx{dim}=[nn-k+1:nn,1:nn-k];%向上平移30个元素
% idx{dim}=[k+1:nn 1:k];
% ffy1=ffy1(idx{:});
% 
% dim=2; 
% k=8;
% idx=repmat({':'}, ndims(ffy1), 1); % initialize subscripts
% nn=size(ffy1,dim);   % length along dimension dim                      
% idx{dim}=[nn-k+1:nn,1:nn-k];%向左平移30个元素
% % idx{dim}=[k+1:nn 1:k];
% ffy1=ffy1(idx{:});

ys=zeros(nn);
ys(1:nn-shear,:)=ZZ0(shear+1:nn,:);
yi=(ZZ0-ys)*4*pi/lamda;
ys1=cos(yi).*ffy;
figure,imshow(ys1,[])

yi=(ZZ0-ys)*4*pi/lamda+pi/2;
ys2=cos(yi).*ffy;
figure,imshow(ys2,[])

yi=(ZZ0-ys)*4*pi/lamda+pi;
ys3=cos(yi).*ffy;
figure,imshow(ys3,[])

yi=(ZZ0-ys)*4*pi/lamda+3*pi/2;
ys4=cos(yi).*ffy;
figure,imshow(ys4,[])


I5=ys1;
I6=ys2;
I7=ys3;
I8=ys4;

[width,height]=size(I5);
for i=1:width
    for j=1:height
        if (I5(i,j)==I7(i,j))
            if(I6(i,j)>I8(i,j))
                phasey(i,j)=pi/2;
            elseif(I6(i,j)<I8(i,j))
                phasey(i,j)=-pi/2;
            else
                phasey(i,j)=0;
            end
        elseif (I5(i,j)>I7(i,j))
            if I6(i,j)>=I8(i,j)
                phasey(i,j)=atan((I6(i,j)-I8(i,j))/(I5(i,j)-I7(i,j)));
            else
                phasey(i,j)=-atan((I8(i,j)-I6(i,j))/(I5(i,j)-I7(i,j))); 
            end
        else
            if I6(i,j)>=I8(i,j)
                phasey(i,j)=pi-atan((I6(i,j)-I8(i,j))/(I7(i,j)-I5(i,j)));
            else
                phasey(i,j)=-pi+atan((I8(i,j)-I6(i,j))/(I7(i,j)-I5(i,j))); 
            end
        end
    end
end
figure,imshow(phasey,[]);

phasey(~isfinite(phasey))=0;
PV_py = (max(max(phasey))-min(min(phasey)));
RMS_py =(max(max(std(phasey))));

