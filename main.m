close all
clear
clc
[x]=audioread('FAML_Sa1.wav');
% MFCC implement with Matlab %  
bank=melbankm(24,256,16000,0,0.5,'t'); %Mel滤波器的阶数为24，FFT变换的长度为256，采样频率为16000Hz  
%归一化Mel滤波器组系数  
bank=full(bank); %full() convert sparse matrix to full matrix  
bank=bank/max(bank(:));  
for k=1:12  
    n=0:23;  
    dctcoef(k,:)=cos((2*n+1)*k*pi/(2*24));  
end  
w=1+6*sin(pi*[1:12]./12);%归一化倒谱提升窗口  
w=w/max(w);%预加重滤波器  
xx=double(x);  
xx=filter([1-0.9375],1,xx);%语音信号分帧  
xx=enframe(xx,256,80);%对xx 256点分为一帧  
%计算每帧的MFCC参数  
for i=1:size(xx,1)  
    y=xx(i,:);  
    s=y'.*hamming(256);  
    t=abs(fft(s));%FFT快速傅里叶变换  
    t=t.^2;  
    c1=dctcoef*log(bank*t(1:129));  
    c2=c1.*w';  
    m(i,:)=c2;  
end  
%求一阶差分系数  
dtm=zeros(size(m));  
for i=3:size(m,1)-2  
    dtm(i,:)=-2*m(i-2,:)-m(i-1,:)+m(i+1,:)+2*m(i+2,:);  
end  
dtm=dtm/3;  
%求取二阶差分系数  
dtmm=zeros(size(dtm));  
for i=3:size(dtm,1)-2  
    dtmm(i,:)=-2*dtm(i-2,:)-dtm(i-1,:)+dtm(i+1,:)+2*dtm(i+2,:);  
end  
dtmm=dtmm/3;  
%合并mfcc参数和一阶差分mfcc参数  
ccc=[m dtm dtmm];  
%去除首尾两帧，以为这两帧的一阶差分参数为0  
ccc=ccc(3:size(m,1)-2,:);  
ccc; 
% subplot(2,1,1);  
% ccc_1=ccc(:,1); 
% m1=m(:,1)
figure(1)
plot(m);xlabel('Number of frames');ylabel('Amplitude'); 
figure(2)
plot(ccc);xlabel('Number of frames');ylabel('Amplitude');
  