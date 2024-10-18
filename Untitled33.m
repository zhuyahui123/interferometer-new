clear 
clc
close all;
lamda=0.6328*10^-3;
x=linspace(-10,10,1000);
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
fn=(x-0.3).^2+(y-0.6).^2;
for i=1:nn
     for j=1:nn
         if fn(i,j)<=(0.82*10)^2
            ff(i,j)=1;
         end
     end
end
figure,imshow(ff,[]);
% ff=zeros(nn);
% fn=x.^2+y.^2;
% for i=1:nn
%      for j=1:nn
%          if fn(i,j)<=(0.78*10)^2
%             ff(i,j)=1;
%          end
%      end
% end
% shear=nn*0.8*0.1;
shear=90;
ff1=zeros(nn);
ff1(:,1:nn-shear)=ff(:,shear+1:nn);
ffx=ff & ff1;
ffx1=ff | ff1;
ffxx=ffx1-ffx;
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
X1=imread('H:\1.实验\4.13chuli_data\1X0.bmp');
% X1=imcrop(X1,[40 40 473 473]);
X1=double(X1(:,:,1));
X1=X1.*ffx;
figure,imshow(X1,[])
X2=imread('H:\1.实验\4.13chuli_data\1X45.bmp');
% X2=imcrop(X2,[40 40 473 473]);
X2=double(X2(:,:,1));
X2=X2.*ffx;
figure,imshow(X2,[])
X3=imread('H:\1.实验\4.13chuli_data\1X90.bmp');
% X3=imcrop(X3,[40 40 473 473]);
X3=double(X3(:,:,1));
X3=X3.*ffx;
figure,imshow(X3,[])
X4=imread('H:\1.实验\4.13chuli_data\1X135.bmp');
% X4=imcrop(X4,[40 40 473 473]);
X4=double(X4(:,:,1));
X4=X4.*ffx;
figure,imshow(X4,[])
Y1=imread('H:\1.实验\4.13chuli_data\1Y0.bmp');
% Y1=imcrop(Y1,[65 40 473 473]);
Y1=double(Y1(:,:,1));
Y1=Y1.*ffy;
figure,imshow(Y1,[])
Y2=imread('H:\1.实验\4.13chuli_data\1Y45.bmp');
% Y2=imcrop(Y2,[65 40 473 473]);
Y2=double(Y2(:,:,1));
Y2=Y2.*ffy;
figure,imshow(Y2,[])
Y3=imread('H:\1.实验\4.13chuli_data\1Y90.bmp');
% Y3=imcrop(Y3,[65 40 473 473]);
Y3=double(Y3(:,:,1));
Y3=Y3.*ffy;
figure,imshow(Y3,[])
Y4=imread('H:\1.实验\4.13chuli_data\1Y135.bmp');
% Y4=imcrop(Y4,[65 40 473 473]);
Y4=double(Y4(:,:,1));
Y4=Y4.*ffy;
figure,imshow(Y4,[])
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
I1=medfilt2(I1,[8,8]);%中值滤波
I2=medfilt2(I2,[8,8]);
I3=medfilt2(I3,[8,8]);
I4=medfilt2(I4,[8,8]);

%解包
% I1=xs1;
% I2=xs2;
% I3=xs3;
% I4=xs4;
clear i
[width,height]=size(I1);
phasex=atan2((I4-I2),(I3-I1));
% for i=1:width
%     for j=1:height
%         if (I1(i,j)==I3(i,j))
%             if(I2(i,j)>I4(i,j))
%                 phasex(i,j)=pi/2;
%             elseif(I2(i,j)<I4(i,j))
%                 phasex(i,j)=-pi/2;
%             else
%                 phasex(i,j)=0;
%             end
%         elseif (I1(i,j)>I3(i,j))
%             if I2(i,j)>=I4(i,j)
%                 phasex(i,j)=atan((I2(i,j)-I4(i,j))/(I1(i,j)-I3(i,j)));
%             else
%                 phasex(i,j)=-atan((I4(i,j)-I2(i,j))/(I1(i,j)-I3(i,j))); 
%             end
%         else
%             if I2(i,j)>=I4(i,j)
%                 phasex(i,j)=pi-atan((I2(i,j)-I4(i,j))/(I3(i,j)-I1(i,j)));
%             else
%                 phasex(i,j)=-pi+atan((I4(i,j)-I2(i,j))/(I3(i,j)-I1(i,j))); 
%             end
%         end
%     end
% end
figure,imshow(phasex,[]);

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
phasex1=phasex.*ffx; 
figure,imshow(phasex1,[])
figure,mesh(phasex1)


I5=double(Z5);
I6=double(Z6);
I7=double(Z7);
I8=double(Z8);
I5=medfilt2(I5,[8,8]);%中值滤波
I6=medfilt2(I6,[8,8]);
I7=medfilt2(I7,[8,8]);
I8=medfilt2(I8,[8,8]);
% I5=ys1;
% I6=ys2;
% I7=ys3;
% I8=ys4;

[width,height]=size(I5);
%解包
phasey=atan2((I8-I6),(I7-I5));
% for i=1:width
%     for j=1:height
%         if (I5(i,j)==I7(i,j))
%             if(I6(i,j)>I8(i,j))
%                 phasey(i,j)=pi/2;
%             elseif(I6(i,j)<I8(i,j))
%                 phasey(i,j)=-pi/2;
%             else
%                 phasey(i,j)=0;
%             end
%         elseif (I5(i,j)>I7(i,j))
%             if I6(i,j)>=I8(i,j)
%                 phasey(i,j)=atan((I6(i,j)-I8(i,j))/(I5(i,j)-I7(i,j)));
%             else
%                 phasey(i,j)=-atan((I8(i,j)-I6(i,j))/(I5(i,j)-I7(i,j))); 
%             end
%         else
%             if I6(i,j)>=I8(i,j)
%                 phasey(i,j)=pi-atan((I6(i,j)-I8(i,j))/(I7(i,j)-I5(i,j)));
%             else
%                 phasey(i,j)=-pi+atan((I8(i,j)-I6(i,j))/(I7(i,j)-I5(i,j))); 
%             end
%         end
%     end
% end
figure,imshow(phasey,[]);

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
phasey1=phasey.*ffy;
figure,imshow(phasey1,[]);
figure,mesh(phasey1)


wx=lamda*phasex1;
wy=lamda*phasey1;
Sx=wx(:)/(2*shear);
Sy=wy(:)/(2*shear);
% Sx=wx(:);
% Sy=wy(:);%按列排一列,若按行排先加‘求转置
S=[Sx;Sy];
%% 差分zernike
% x=linspace(-1,1,1000);
% [x,y]=meshgrid(x,x);
x=linspace(-10,10,1000);
[x,y]=meshgrid(x,x);
[theta,rho] = cart2pol(x,y);%将笛卡尔坐标转换为极坐标。
index_rr = find(rho<=0.8);
M =nan(n,n);
M(index_rr)=1;
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
% Z1=Z1.*ff;
% Z2=Z2.*ff;
% Z3=Z3.*ff;
% Z4=Z4.*ff;
% Z5=Z5.*ff;
% Z6=Z6.*ff;
% Z7=Z7.*ff;
% Z8=Z8.*ff;
% Z9=Z9.*ff;
% Z10=Z10.*ff;
% Z11=Z11.*ff;
% Z12=Z12.*ff;
% Z13=Z13.*ff;
% Z14=Z14.*ff;
% Z15=Z15.*ff;
% Z16=Z16.*ff;
% Z17=Z17.*ff;
% Z18=Z18.*ff;
% Z19=Z19.*ff;
% Z20=Z20.*ff;
% Z21=Z21.*ff;
% Z22=Z22.*ff;
% Z23=Z23.*ff;
% Z24=Z24.*ff;
% Z25=Z25.*ff;
% Z26=Z26.*ff;
% Z27=Z27.*ff;
% Z28=Z28.*ff;
% Z28=Z28.*ff;
% Z29=Z29.*ff;
% Z30=Z30.*ff;
% Z31=Z31.*ff;
% Z32=Z32.*ff;
% Z33=Z33.*ff;
% Z34=Z34.*ff;
% Z35=Z35.*ff;
% Z36=Z36.*ff;
Z1=Z1.*M;
Z2=Z2.*M;
Z3=Z3.*M;
Z4=Z4.*M;
Z5=Z5.*M;
Z6=Z6.*M;
Z7=Z7.*M;
Z8=Z8.*M;
Z9=Z9.*M;
Z10=Z10.*M;
Z11=Z11.*M;
Z12=Z12.*M;
Z13=Z13.*M;
Z14=Z14.*M;
Z15=Z15.*M;
Z16=Z16.*M;
Z17=Z17.*M;
Z18=Z18.*M;
Z19=Z19.*M;
Z20=Z20.*M;
Z21=Z21.*M;
Z22=Z22.*M;
Z23=Z23.*M;
Z24=Z24.*M;
Z25=Z25.*M;
Z26=Z26.*M;
Z27=Z27.*M;
Z28=Z28.*M;
Z28=Z28.*M;
Z29=Z29.*M;
Z30=Z30.*M;
Z31=Z31.*M;
Z32=Z32.*M;
Z33=Z33.*M;
Z34=Z34.*M;
Z35=Z35.*M;
Z36=Z36.*M;
%% 差分zernike1
Zs1=zeros(nn);
Zs2=zeros(nn);
Zs3=zeros(nn);
Zs4=zeros(nn);
Zs5=zeros(nn);
Zs6=zeros(nn);
Zs7=zeros(nn);
Zs8=zeros(nn);
Zs9=zeros(nn);
Zs10=zeros(nn);
Zs11=zeros(nn);
Zs12=zeros(nn);
Zs13=zeros(nn);
Zs14=zeros(nn);
Zs15=zeros(nn);
Zs16=zeros(nn);
Zs17=zeros(nn);
Zs18=zeros(nn);
Zs19=zeros(nn);
Zs20=zeros(nn);
Zs21=zeros(nn);
Zs22=zeros(nn);
Zs23=zeros(nn);
Zs24=zeros(nn);
Zs25=zeros(nn);
Zs26=zeros(nn);
Zs27=zeros(nn);
Zs28=zeros(nn);
Zs29=zeros(nn);
Zs30=zeros(nn);
Zs31=zeros(nn);
Zs32=zeros(nn);
Zs33=zeros(nn);
Zs34=zeros(nn);
Zs35=zeros(nn);
Zs36=zeros(nn);
Zs1=Zs1.*M;
Zs2=Zs2.*M;
Zs3=Zs3.*M;
Zs4=Zs4.*M;
Zs5=Zs5.*M;
Zs6=Zs6.*M;
Zs7=Zs7.*M;
Zs8=Zs8.*M;
Zs9=Zs9.*M;
Zs10=Zs10.*M;
Zs11=Zs11.*M;
Zs12=Zs12.*M;
Zs13=Zs13.*M;
Zs14=Zs14.*M;
Zs15=Zs15.*M;
Zs16=Zs16.*M;
Zs17=Zs17.*M;
Zs18=Zs18.*M;
Zs19=Zs19.*M;
Zs20=Zs20.*M;
Zs21=Zs21.*M;
Zs22=Zs22.*M;
Zs23=Zs23.*M;
Zs24=Zs24.*M;
Zs25=Zs25.*M;
Zs26=Zs26.*M;
Zs27=Zs27.*M;
Zs28=Zs28.*M;
Zs28=Zs28.*M;
Zs29=Zs29.*M;
Zs30=Zs30.*M;
Zs31=Zs31.*M;
Zs32=Zs32.*M;
Zs33=Zs33.*M;
Zs34=Zs34.*M;
Zs35=Zs35.*M;
Zs36=Zs36.*M;
%
% 差分zernike1
Zs1(:,1:nn-shear)=Z1(:,shear+1:nn);
Zs2(:,1:nn-shear)=Z2(:,shear+1:nn);
Zs3(:,1:nn-shear)=Z3(:,shear+1:nn);
Zs4(:,1:nn-shear)=Z4(:,shear+1:nn);
Zs5(:,1:nn-shear)=Z5(:,shear+1:nn);
Zs6(:,1:nn-shear)=Z6(:,shear+1:nn);
Zs7(:,1:nn-shear)=Z7(:,shear+1:nn);
Zs8(:,1:nn-shear)=Z8(:,shear+1:nn);
Zs9(:,1:nn-shear)=Z9(:,shear+1:nn);
Zs10(:,1:nn-shear)=Z10(:,shear+1:nn);
Zs11(:,1:nn-shear)=Z11(:,shear+1:nn);
Zs12(:,1:nn-shear)=Z12(:,shear+1:nn);
Zs13(:,1:nn-shear)=Z13(:,shear+1:nn);
Zs14(:,1:nn-shear)=Z14(:,shear+1:nn);
Zs15(:,1:nn-shear)=Z15(:,shear+1:nn);
Zs16(:,1:nn-shear)=Z16(:,shear+1:nn);
Zs17(:,1:nn-shear)=Z17(:,shear+1:nn);
Zs18(:,1:nn-shear)=Z18(:,shear+1:nn);
Zs19(:,1:nn-shear)=Z19(:,shear+1:nn);
Zs20(:,1:nn-shear)=Z20(:,shear+1:nn);
Zs21(:,1:nn-shear)=Z21(:,shear+1:nn);
Zs22(:,1:nn-shear)=Z22(:,shear+1:nn);
Zs23(:,1:nn-shear)=Z23(:,shear+1:nn);
Zs24(:,1:nn-shear)=Z24(:,shear+1:nn);
Zs25(:,1:nn-shear)=Z25(:,shear+1:nn);
Zs26(:,1:nn-shear)=Z26(:,shear+1:nn);
Zs27(:,1:nn-shear)=Z27(:,shear+1:nn);
Zs28(:,1:nn-shear)=Z28(:,shear+1:nn);
Zs29(:,1:nn-shear)=Z29(:,shear+1:nn);
Zs30(:,1:nn-shear)=Z30(:,shear+1:nn);
Zs31(:,1:nn-shear)=Z31(:,shear+1:nn);
Zs32(:,1:nn-shear)=Z32(:,shear+1:nn);
Zs33(:,1:nn-shear)=Z33(:,shear+1:nn);
Zs34(:,1:nn-shear)=Z34(:,shear+1:nn);
Zs45(:,1:nn-shear)=Z35(:,shear+1:nn);
Zs36(:,1:nn-shear)=Z36(:,shear+1:nn);
Z1x=(Zs1-Z1).*ffx;
Z2x=(Zs2-Z2).*ffx;
Z3x=(Zs3-Z3).*ffx;
Z4x=(Zs4-Z4).*ffx;
Z5x=(Zs5-Z5).*ffx;
Z6x=(Zs6-Z6).*ffx;
Z7x=(Zs7-Z7).*ffx;
Z8x=(Zs8-Z8).*ffx;
Z9x=(Zs9-Z9).*ffx;
Z10x=(Zs10-Z10).*ffx;
Z11x=(Zs11-Z11).*ffx;
Z12x=(Zs12-Z12).*ffx;
Z13x=(Zs13-Z13).*ffx;
Z14x=(Zs14-Z14).*ffx;
Z15x=(Zs15-Z15).*ffx;
Z16x=(Zs16-Z16).*ffx;
Z17x=(Zs17-Z17).*ffx;
Z18x=(Zs18-Z18).*ffx;
Z19x=(Zs19-Z19).*ffx;
Z20x=(Zs20-Z20).*ffx;
Z21x=(Zs21-Z21).*ffx;
Z22x=(Zs22-Z22).*ffx;
Z23x=(Zs23-Z23).*ffx;
Z24x=(Zs24-Z24).*ffx;
Z25x=(Zs25-Z25).*ffx;
Z26x=(Zs26-Z26).*ffx;
Z27x=(Zs27-Z27).*ffx;
Z28x=(Zs28-Z28).*ffx;
Z29x=(Zs29-Z29).*ffx;
Z30x=(Zs30-Z30).*ffx;
Z31x=(Zs31-Z31).*ffx;
Z32x=(Zs32-Z32).*ffx;
Z33x=(Zs33-Z33).*ffx;
Z34x=(Zs34-Z34).*ffx;
Z35x=(Zs35-Z35).*ffx;
Z36x=(Zs36-Z36).*ffx;
Zi1=zeros(nn);
Zi2=zeros(nn);
Zi3=zeros(nn);
Zi4=zeros(nn);
Zi5=zeros(nn);
Zi6=zeros(nn);
Zi7=zeros(nn);
Zi8=zeros(nn);
Zi9=zeros(nn);
Zi10=zeros(nn);
Zi11=zeros(nn);
Zi12=zeros(nn);
Zi13=zeros(nn);
Zi14=zeros(nn);
Zi15=zeros(nn);
Zi16=zeros(nn);
Zi17=zeros(nn);
Zi18=zeros(nn);
Zi19=zeros(nn);
Zi20=zeros(nn);
Zi21=zeros(nn);
Zi22=zeros(nn);
Zi23=zeros(nn);
Zi24=zeros(nn);
Zi25=zeros(nn);
Zi26=zeros(nn);
Zi27=zeros(nn);
Zi28=zeros(nn);
Zi29=zeros(nn);
Zi30=zeros(nn);
Zi31=zeros(nn);
Zi32=zeros(nn);
Zi33=zeros(nn);
Zi34=zeros(nn);
Zi35=zeros(nn);
Zi36=zeros(nn);
Zi1=Zi1.*M;
Zi2=Zi2.*M;
Zi3=Zi3.*M;
Zi4=Zi4.*M;
Zi5=Zi5.*M;
Zi6=Zi6.*M;
Zi7=Zi7.*M;
Zi8=Zi8.*M;
Zi9=Zi9.*M;
Zi10=Zi10.*M;
Zi11=Zi11.*M;
Zi12=Zi12.*M;
Zi13=Zi13.*M;
Zi14=Zi14.*M;
Zi15=Zi15.*M;
Zi16=Zi16.*M;
Zi17=Zi17.*M;
Zi18=Zi18.*M;
Zi19=Zi19.*M;
Zi20=Zi20.*M;
Zi21=Zi21.*M;
Zi22=Zi22.*M;
Zi23=Zi23.*M;
Zi24=Zi24.*M;
Zi25=Zi25.*M;
Zi26=Zi26.*M;
Zi27=Zi27.*M;
Zi28=Zi28.*M;
Zi28=Zi28.*M;
Zi29=Zi29.*M;
Zi30=Zi30.*M;
Zi31=Zi31.*M;
Zi32=Zi32.*M;
Zi33=Zi33.*M;
Zi34=Zi34.*M;
Zi35=Zi35.*M;
Zi36=Zi36.*M;
Zi1(1:nn-shear,:)=Z1(shear+1:nn,:);
Zi2(1:nn-shear,:)=Z2(shear+1:nn,:);
Zi3(1:nn-shear,:)=Z3(shear+1:nn,:);
Zi4(1:nn-shear,:)=Z4(shear+1:nn,:);
Zi5(1:nn-shear,:)=Z5(shear+1:nn,:);
Zi6(1:nn-shear,:)=Z6(shear+1:nn,:);
Zi7(1:nn-shear,:)=Z7(shear+1:nn,:);
Zi8(1:nn-shear,:)=Z8(shear+1:nn,:);
Zi9(1:nn-shear,:)=Z9(shear+1:nn,:);
Zi10(1:nn-shear,:)=Z10(shear+1:nn,:);
Zi11(1:nn-shear,:)=Z11(shear+1:nn,:);
Zi12(1:nn-shear,:)=Z12(shear+1:nn,:);
Zi13(1:nn-shear,:)=Z13(shear+1:nn,:);
Zi14(1:nn-shear,:)=Z14(shear+1:nn,:);
Zi15(1:nn-shear,:)=Z15(shear+1:nn,:);
Zi16(1:nn-shear,:)=Z16(shear+1:nn,:);
Zi17(1:nn-shear,:)=Z17(shear+1:nn,:);
Zi18(1:nn-shear,:)=Z18(shear+1:nn,:);
Zi19(1:nn-shear,:)=Z19(shear+1:nn,:);
Zi20(1:nn-shear,:)=Z20(shear+1:nn,:);
Zi21(1:nn-shear,:)=Z21(shear+1:nn,:);
Zi22(1:nn-shear,:)=Z22(shear+1:nn,:);
Zi23(1:nn-shear,:)=Z23(shear+1:nn,:);
Zi24(1:nn-shear,:)=Z24(shear+1:nn,:);
Zi25(1:nn-shear,:)=Z25(shear+1:nn,:);
Zi26(1:nn-shear,:)=Z26(shear+1:nn,:);
Zi27(1:nn-shear,:)=Z27(shear+1:nn,:);
Zi28(1:nn-shear,:)=Z28(shear+1:nn,:);
Zi29(1:nn-shear,:)=Z29(shear+1:nn,:);
Zi30(1:nn-shear,:)=Z30(shear+1:nn,:);
Zi31(1:nn-shear,:)=Z31(shear+1:nn,:);
Zi32(1:nn-shear,:)=Z32(shear+1:nn,:);
Zi33(1:nn-shear,:)=Z33(shear+1:nn,:);
Zi34(1:nn-shear,:)=Z34(shear+1:nn,:);
Zi35(1:nn-shear,:)=Z35(shear+1:nn,:);
Zi36(1:nn-shear,:)=Z36(shear+1:nn,:);
Z1y=(Zi1-Z1).*ffy;
Z2y=(Zi2-Z2).*ffy;
Z3y=(Zi3-Z3).*ffy;
Z4y=(Zi4-Z4).*ffy;
Z5y=(Zi5-Z5).*ffy;
Z6y=(Zi6-Z6).*ffy;
Z7y=(Zi7-Z7).*ffy;
Z8y=(Zi8-Z8).*ffy;
Z9y=(Zi9-Z9).*ffy;
Z10y=(Zi10-Z10).*ffy;
Z11y=(Zi11-Z11).*ffy;
Z12y=(Zi12-Z12).*ffy;
Z13y=(Zi13-Z13).*ffy;
Z14y=(Zi14-Z14).*ffy;
Z15y=(Zi15-Z15).*ffy;
Z16y=(Zi16-Z16).*ffy;
Z17y=(Zi17-Z17).*ffy;
Z18y=(Zi18-Z18).*ffy;
Z19y=(Zi19-Z19).*ffy;
Z20y=(Zi20-Z20).*ffy;
Z21y=(Zi21-Z21).*ffy;
Z22y=(Zi22-Z22).*ffy;
Z23y=(Zi23-Z23).*ffy;
Z24y=(Zi24-Z24).*ffy;
Z25y=(Zi25-Z25).*ffy;
Z26y=(Zi26-Z26).*ffy;
Z27y=(Zi27-Z27).*ffy;
Z28y=(Zi28-Z28).*ffy;
Z29y=(Zi29-Z29).*ffy;
Z30y=(Zi30-Z30).*ffy;
Z31y=(Zi31-Z31).*ffy;
Z32y=(Zi32-Z32).*ffy;
Z33y=(Zi33-Z33).*ffy;
Z34y=(Zi34-Z34).*ffy;
Z35y=(Zi35-Z35).*ffy;
Z36y=(Zi36-Z36).*ffy;
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
% 面形拟合
ZZ0=(aa(2).*Z2+aa(3).*Z3+aa(4).*Z4+aa(5).*Z5+aa(6).*Z6+aa(7).*Z7+aa(8).*Z8+aa(9).*Z9+aa(10).*Z10+aa(11).*Z11+aa(12).*Z12+aa(13).*Z13+aa(14).*Z14+aa(15).*Z15+aa(16).*Z16+aa(17).*Z17+aa(18).*Z18+aa(19).*Z19+aa(20).*Z20+aa(21).*Z21+aa(22).*Z22+aa(23).*Z23+aa(24).*Z24+aa(25).*Z25+aa(26).*Z26+aa(27).*Z27+aa(28).*Z28+aa(29).*Z29+aa(30).*Z30+aa(31).*Z31+aa(32).*Z32+aa(33).*Z33+aa(34).*Z34+aa(35).*Z35+aa(36).*Z36);
ZZ0=ZZ0-min(min(ZZ0));
ZZ0=ZZ0.*ff;
%拟合面型
figure,mesh(x,y,ZZ0),title('拟合面形');
xlabel('x(cm)');ylabel('y(cm)');zlabel('z（mm）');

% xlabel('x(mm)');ylabel('y(mm)');zlabel('z（lamda）');
ZZ0(~isfinite(ZZ0))=0;
% 
% Z=10.*ZZ0-x1;
PV_Z = (max(max(ZZ0))-min(min(ZZ0)))/lamda;
RMS_Z =(max(max(std(ZZ0))))/lamda;

% %残差
% ZZc=ZZ0-x1;
% figure,mesh(x,y,ZZc/(2*lamda)),title('残差面形');
% xlabel('x(cm)');ylabel('y(cm)');zlabel('z（mm）');
% ZZc(~isfinite(ZZc))=0;
% 
% PV_c = (max(max(ZZc))-min(min(ZZc)))/lamda;
% RMS_c =(max(max(std(ZZc))))/lamda;
