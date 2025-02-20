 % wide band ambiguity function
clc
clear all
close all
%

fs = 4e3;
t = 0:1/fs:2-1/fs;%2s
TIME=2;

tb=-TIME;
tend=TIME;
dialation=.8;%扩张系数  <1扩张，>1收缩



signal=[zeros(1,round(TIME/4*fs)) exp(1j*2*pi*fs*t) zeros(1,TIME/4*fs)];

conj_rev_signal=conj(signal(end:-1:1));
time=linspace(tb,tend,length(signal));
delay_time=linspace(2*tb,2*tend,2*length(signal)-1);
lambda=linspace(dialation,1/dialation,512);
out=zeros(length(lambda),2*length(signal)-1);
for i=1:length(lambda)
out(i,:)=conv(signal,interp1(time,conj_rev_signal,time*lambda(i),'spline'));
end
figure(1)
imagesc(delay_time,lambda,abs(out))
xlabel('\tau(s)');ylabel('\lambda');
