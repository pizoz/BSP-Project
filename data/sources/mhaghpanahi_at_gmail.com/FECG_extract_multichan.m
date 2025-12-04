function [Mout,Fout,corrcoef]=FECG_extract_multichan(Fin,fs,peaks,L)

[nleads,nsamples]=size(Fin);
phase = PhaseCalculation(peaks);     % phase calculation

JJ = find(peaks);
fm = fs./diff(JJ);          % heart-rate
w = mean(2*pi*fm);          % average heart-rate in rads.
wsd = std(2*pi*fm,1);       % heart-rate standard deviation in rads.


%%
% Pi-ICA analysis
[indep_comp,W,A] = PiCA([Fin -Fin],[peaks peaks]);

indep_comp = indep_comp(:,1:nsamples);
num_comp=size(indep_comp,1);


%%
% MECG Extraction

bins = 250; % number of phase bins
ECGsd= zeros(L,bins);
nk= zeros(1,L); % a vector indicating number of Gaussian kernels for each independent component
Q= (wsd)^2; % the process noise covariance matrix initialized with its first component. Note: size of Q changes in the for loop.
Wmean= w; % the process noise mean vector initialized with its first component. Note: size of Wmean changes in the for loop.

for comp_ind=1:L
    [ECGmean,ECGsd(comp_ind,:),meanphase] = MeanECGExtraction(indep_comp(comp_ind,:),phase,bins,1); % mean ECG extraction
    ECGsmooth=smooth(ECGmean,6,'moving')';
     
     test1=diff(ECGsmooth);
     test1=sign(test1);
     
     test2=filter([1 1],1,test1);
     gauss_ind=find(test2==0);
     lowchange=abs(diff(ECGsmooth(gauss_ind)))<1e-3; % ECGmean replaced by ECGsmooth 
     gauss_ind=setdiff(gauss_ind,gauss_ind(find(lowchange==1)));
     
    
    theta = meanphase(gauss_ind);
    alpha = 1.2*ECGsmooth(gauss_ind);% ECGmean replaced by ECGsmooth
    b = .04*ones(size(alpha));
    
    InitParams1= [alpha b theta];

    options = optimset('TolX',1e-4,'TolFun',1e-4,'MaxIter',100);
    OptimumParams = nlinfit(meanphase,ECGmean,@ECGModel,InitParams1,options);
    
%         Model = ECGModelError(OptimumParams,ECGmean,meanphase,1);
% 
%     figure;
%     plot(Model,'k','LineWidth',2);
%     hold on
%     plot(ECGmean,'r');
%     plot(ECGsmooth);
%     plot(gauss_ind,ECGsmooth(gauss_ind),'*');
    
    
 
% KF modeling and Parameter setting
    N = length(OptimumParams)/3;% number of Gaussian kernels
    nk(comp_ind)= N;
    subQ= [(.1*OptimumParams(1:N)).^2; (.05*ones(1,N)).^2; (.05*ones(1,N)).^2];
    Q=[Q zeros(size(Q,1),N*3); zeros(N*3, size(Q,2)) diag(subQ(:))]; 
    kp= reshape(OptimumParams,[],3)'; % Kernel parameters reshaped into a matrix. 1st row: alpha, 2nd row: b, 3rd row: theta
    kp= [1 0 0;0 0 1; 0 1 0]*kp; % Kernel parameters reshaped into a matrix. 1st row: alpha, 2nd row: theta, 3rd row: b
    Wmean= [Wmean; kp(:)];    
end

Y = [phase ; indep_comp(1:L,:)]; % Observation matrix
X0 = [-pi zeros(1,L)]'; % mean of the initial state
P0 = diag([(2*pi)^2; (10*max(abs(indep_comp(1:L,:)),[],2)).^2]);
Q=[Q zeros(size(Q,1),L); zeros(L, size(Q,2)) diag((.05*mean(ECGsd(:,1:round(bins/10)),2)).^2)];
R = [(w/fs).^2/12 zeros(1,L) ;zeros(L,1) diag((mean(ECGsd(:,1:round(bins/10)),2)).^2)];
Wmean= [Wmean; zeros(L,1)];
Vmean= zeros(L+1,1);

% Reserve space for estimates.
Xminus = zeros(size(X0,1),size(Y,2));
Xekf = zeros(size(X0,1),size(Y,2));
Pminus = zeros(size(X0,1),size(X0,1),size(Y,2));
Pekf = zeros(size(X0,1),size(X0,1),size(Y,2));

X= X0;
P= P0;

params={Wmean, nk, fs};

% Estimate with EKF
for k=1:size(Y,2)

   
   [X,P] = ekf_predict1(X,P,ekf_df_dx(X,params),Q,ekf_cost(X,params), ekf_df_dw(X,params));
   
    if(abs(X(1)-Y(1,k)) > pi)
       X(1) = Y(1,k);
    end 
   
   Xminus(:,k)   = X;
   Pminus(:,:,k) = P;
   
   [X,P] = ekf_update1(X,P,Y(:,k),eye(L+1),R);
   Xekf(:,k)   = X;
   Pekf(:,:,k) = P;
     
end



[Xeks,Peks] = ekf_smooth1(Xekf,Pekf,Xminus,Pminus,@ekf_df_dx,params);

indep_comp_cleaned= indep_comp(1:L,:)-Xeks(2:L+1,:);
MECG_comp_cleaned= Xeks(2:L+1,:);



%%
% Inverse Periodic-ICA

FICA=[indep_comp_cleaned;indep_comp(L+1:4,:)];
MICA=[MECG_comp_cleaned;zeros(num_comp-L,nsamples)];
Fout=A*FICA; 
Mout=A*MICA;

corrcoef = zeros(4,1);
for i=1:4
    xx=dualSample(Fout(i,:),peaks,phase);
    corrcoef(i)=mean(xx.*Fout(i,1:length(xx)))/sqrt(mean(xx.^2)*mean(Fout(i,1:length(xx)).^2));
end

end








    
   