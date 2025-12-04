function Z=wavelet2(x)
debug=0;
x=x(1:60000); 
t=1:60000; 
n=1:60000;
fs=1000; 
show=0;
N=length(x);
if debug
figure(1); 
plot(t,x); 
grid;
xlabel('时间(ms)','fontsize',8); 
ylabel('幅值v)','fontsize',8); 
title('原始信号','fontsize',8); 
end
%低通去除基线漂移
k = .7;
% cut-off value
fc=0.8/fs; 
alpha = (1-k*cos(2*pi*fc)-sqrt(2*k*(1-cos(2*pi*fc))-k^2*sin(2*pi*fc)^2))/(1-k);
y = zeros(size(x)); 
for i = 1:size(x,1)  %可以去掉   
    y(i,:) = filtfilt(1-alpha,[1 -alpha],x(i,:));
end
x1=x-y;
%figure(2) 
%plot(x1); 
%grid;
%xlabel('时间/ms');
%ylabel('幅值/mv'); 
%title('去除基线漂移后的信号'); 
%axis tight; 
%figure(3) 
%plot(x1-x); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%小波去噪%%%%%%%%%%%%%%%%%%%% 
%小波分解%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[c,l]=wavedec(x1,5,'db5');
%小波分解
if(show)
%a6=wrcoef('a',c,l,'db6',6);
a5=wrcoef('a',c,l,'db5',5);
%a4=wrcoef('a',c,l,'db6',5);
%d6=wrcoef('d',c,l,'db6',6); 
d5=wrcoef('d',c,l,'db5',5); 
d4=wrcoef('d',c,l,'db5',4); 
d3=wrcoef('d',c,l,'db5',3); 
d2=wrcoef('d',c,l,'db5',2);
d1=wrcoef('d',c,l,'db5',1);
%显示各层信号 
figure(4) 
subplot(3,2,1) 
plot((1:N)/fs,a5,'LineWidth',1);
%plot((1:N)/fs,a6,'LineWidth',1);
grid; 
ylabel('a5');
%ylabel('a6');
subplot(3,2,2); 
plot((1:N)/fs,d5,'LineWidth',1);
grid; 
ylabel('d5');
%plot((1:N)/fs,d6,'LineWidth',1);
%grid; 
%ylabel('d6');
subplot(3,2,3);
%plot((1:N)/fs,d5,'LineWidth',1);
%grid; ylabel('d5');
%xlabel('时间 t/s'); 
%subplot(3,2,4);
plot((1:N)/fs,d4,'LineWidth',1);
grid; ylabel('d4');
xlabel('时间 t/s'); 
subplot(3,2,4); 
plot((1:N)/fs,d3,'LineWidth',1);
grid; ylabel('d3');
subplot(3,2,5);
plot((1:N)/fs,d2,'LineWidth',1);
grid; ylabel('d2'); 
subplot(3,2,6);
plot((1:N)/fs,d1,'LineWidth',1);
grid; ylabel('d1'); 
xlabel('时间/s'); 
grid; 
%信号重构
figure(5) 
y=a5+d5; 
%y=a6+d6; 
plot((1:N)/fs,y,'LineWidth',1);
grid; xlabel('时间 t/s'); 
ylabel('电压 A/uV');
title('降噪后的信号'); 
diffence=x-y;
%重构误差
figure(6); 
plot(diffence);
grid; 
title('小波重构误差'); 
%%%%%%%%%%%%%%%%%%%多尺度分解提取低频和高频系数%%%%%%%%%%%%%%%%%%%%%%% 
%提取高频系数
ca1=appcoef(c,l,'db5',1); 
ca2=appcoef(c,l,'db5',2); 
ca3=appcoef(c,l,'db5',3); 
ca4=appcoef(c,l,'db5',4); 
ca5=appcoef(c,l,'db5',5); 
%ca6=appcoef(c,l,'db6',6); 
figure(7); title('低频系数');
subplot(5,1,1);
plot(ca1);
grid; 
title('尺度1的低频系数'); 
subplot(5,1,2); 
plot(ca2);
grid; 
title('尺度2的低频系数');
subplot(5,1,3); 
plot(ca3);
grid; 
title('尺度3的低频系数');
subplot(5,1,4);
plot(ca4);
grid; 
title('尺度4的低频系数'); 
subplot(5,1,5);
plot(ca5);
grid;
title('尺度5的低频系数');
%subplot(6,1,6);
%plot(ca6);
%grid;
%title('尺度6的低频系数');
%提取高频系数
cd1=detcoef(c,l,1); 
cd2=detcoef(c,l,2);
cd3=detcoef(c,l,3); 
cd4=detcoef(c,l,4); 
cd5=detcoef(c,l,5); 
%cd6=detcoef(c,l,6); 
figure(8);
title('高频系数'); 
subplot(5,1,1);
plot(cd1);
grid;
title('尺度1的高频系数');
subplot(5,1,2); 
plot(cd2);
grid; 
title('尺度2的高频系数'); 
subplot(5,1,3);
plot(cd3);
grid;
title('尺度3的高频系数');
subplot(5,1,4); 
plot(cd4);
grid; 
title('尺度4的高频系数'); 
subplot(5,1,5); 
plot(cd5);
grid;
title('尺度5的高频系数');
%subplot(6,1,6); 
%plot(cd6);
%grid;
%title('尺度6的高频系数');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%小波压缩%%%%%%%%%%%%%%%%%%%%%%%%% 
alpha=1.5; 
sorh='h'; 
[thr,nkeep]=wdcbm(c,l,alpha);
[xc,cxc,lxc,perf0,perfl2]=wdencmp('lvd',c,l,'db5',5,thr,sorh);
%figure(9); 
%subplot(2,1,1);
%plot(x);
%grid;
%title('原始信号') 
%subplot(2,1,2); 
% plot(xc);
%grid; 
%title('压缩后信号');  
%figure(10);
%plot(xc-x1);
%figure(11);
%plot(xc,'r');hold on,plot(x1);
%figure(12);
% plot(xc,'r');
Z=xc;
