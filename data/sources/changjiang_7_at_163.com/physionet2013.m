function [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG)
% Template algorithm for Physionet/CinC competition 2013. This function can
% be used for events 1 and 2. Participants are free to modify any
% components of the code. However the function prototype must stay the
% same:
%
% [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG) where the inputs and outputs are specified
% below.
%
% inputs:
%   ECG: 4x60000 (4 channels and 1min of signal at 1000Hz) matrix of
%   abdominal ECG channels.
%   tm : Nx1 vector of time in milliseconds
% output:
%   FQRS: FQRS markers in seconds. Each marker indicates the position of one
%   of the FQRS detected by the algorithm.
%   QT_Interval:   1x1 estimated fetal QT duration (enter NaN or 0 if you do wish to calculate)
%
%
% Author: Joachim Behar - IPMG Oxford (joachim.behar@eng.ox.ac.uk)
% Last updated: March 3, 2013 Ikaro Silva


% ---- check size of ECG ----
if size(ECG,2)>size(ECG,1)
    ECG = ECG';
end

fs      = 1000;             % sampling frequency
N       = size(ECG,2);      % number of abdominal channels
debug   = 0;  % enter debug mode?

ECG(find(isnan(ECG)==1)) = 0;
x=ECG(:,4);
x=x'; 
z=wavelet(x);
points=60000;  
level=5;    num_inter=6;   wf='db5';  fs=1000;
clr={'k-','k-','k-','k-','b-'};

[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wf);
[swa,swd] = swt(z,level,Lo_D,Hi_D);
if  debug 
   figure(1);
   subplot(level,1,1); plot(real(z)); grid on;axis tight;
   for i=1:level
       subplot(level+1,2,2*(i)+1);
       plot(swa(i,:)); axis tight;grid on;xlabel('time');
       ylabel(strcat('a   ',num2str(i)));
       subplot(level+1,2,2*(i)+2);
       plot(swd(i,:)); axis tight;grid on;
       ylabel(strcat('d   ',num2str(i)));
   end 
end

ddw=zeros(size(swd));
pddw=ddw;
nddw=ddw;
posw=swd.*(swd>0);
pdw=((posw(:,1:points-1)-posw(:,2:points))<0);
pddw(:,2:points-1)=((pdw(:,1:points-2)-pdw(:,2:points-1))>0);
negw=swd.*(swd<0);
ndw=((negw(:,1:points-1)-negw(:,2:points))>0);
nddw(:,2:points-1)=((ndw(:,1:points-2)-ndw(:,2:points-1))>0);
ddw=pddw|nddw;
ddw(:,1)=1;
ddw(:,points)=1;
wpeak=ddw.*swd;
wpeak(:,1)=wpeak(:,1)+1e-10;
wpeak(:,points)=wpeak(:,points)+1e-10;
if debug
figure(2);
subplot(level+1,1,1); plot(real(z)); grid on;axis tight;
   for i=1:level
       subplot(level+1,1,i+1);
       plot(wpeak(i,:)); axis tight;grid on;
   ylabel(strcat('j=   ',num2str(i)));
   end
end

hlevel_wpeak=wpeak(level,:);
selevel_wpeak=wpeak(4,:);
sig=z;
peaks = PeakDetection(hlevel_wpeak,1/fs);
selevel_wpeak2=selevel_wpeak;
 guilinglet=zeros(1,length(selevel_wpeak2));
 for j=1:length(peaks)
     guilinglet(peaks(j))=1;
 end
D4_p=(guilinglet~=0);
D4_p2=(guilinglet~=0);
O_d4=60;%该参数确定在上一级搜索极大值的范围，可以调整。
O_d42=60;
for P_d4=O_d4:(length(guilinglet)-O_d4);
    if D4_p(P_d4)==1; 
        for i=1:O_d4-1;
        D4_p(P_d4-i)=1;
        
        end ;
    end;     
end;
for P_d4=O_d42:(length(guilinglet)-O_d42);
    if D4_p2(P_d4)==1; 
        for i=1:O_d42-1;
        D4_p(P_d4+i)=1;
        
        end ;
    end;     
end;
selevel_wpeak2=selevel_wpeak2.*D4_p;

selevel_wpeak3=selevel_wpeak-selevel_wpeak2;

 fpeaks= PeakDetection(selevel_wpeak3,2/fs);
 if debug
 figure(10);
 plot(selevel_wpeak3);
xlim([0,10000]);
% hold on,grid on;
%  plot(QRSAnn1,sig(QRSAnn1)+(m-1)*2,'x','MarkerSize',7,...
%                          'MarkerEdgeColor','g','LineWidth',2);
      grid on;
   xlim([0,10000]);                
   hold on,plot(fpeaks,sig(fpeaks)+(m-1)*2,'o','MarkerSize',7,...
                         'MarkerEdgeColor','r','LineWidth',2);grid on;   
 end
 
                     
aa=1;
for i=1:length(fpeaks)
     if  fpeaks(i)>0
         refpeaks(aa)=fpeaks(i);
         aa=aa+1;
     end
end

sortlength=20;
sign= abs(max(z))>abs(min(z));
for j=1:length(refpeaks)
     maxorminpost=refpeaks(j);
     if sign 
        for tt=1:sortlength
            if  refpeaks(j)-tt<=0 
                break;
            else  
                if z(refpeaks(j)-tt)<z(maxorminpost)
                   maxorminpost=refpeaks(j)-tt;
                end 
            end
        end
        for tt=1:sortlength 
             if  refpeaks(j)+tt>=points
                 break;
             else
                 if   z(refpeaks(j)+tt)<z(maxorminpost)
                      maxorminpost=refpeaks(j)+tt;
                 end
             end
        end
     else
        for tt=1:sortlength 
            if  refpeaks(j)-tt<=0 
                break;
            else  
                if z(refpeaks(j)-tt)>z(maxorminpost)
                    maxorminpost=refpeaks(j)-tt;
                end 
            end
        end
        for tt=1:sortlength 
             if  refpeaks(j)+tt>=points
                 break;
             else
                 if   z(refpeaks(j)+tt)>z(maxorminpost)
                      maxorminpost=refpeaks(j)+tt;
                 end
             end
        end 
     end
    refpeaks(j) =maxorminpost;
end
if  debug
 figure(13);
 plot(sig + (m-1)*2,clr{m});
                xlim([0,10000])
%                 hold on;grid on
%      plot(QRSAnn,sig(QRSAnn)+(m-1)*2,'x','MarkerSize',7,...
%                          'MarkerEdgeColor','g','LineWidth',2);
      grid on;
   xlim([0,10000]);                
   hold on,plot(refpeaks,sig(refpeaks)+(m-1)*2,'o','MarkerSize',7,...
                         'MarkerEdgeColor','r','LineWidth',2);grid on; 
end
fetal_QRSAnn_est=refpeaks;
QT_Interval=0;

end



