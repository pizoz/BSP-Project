function df = ekf_df_dw(X,params)

    % X: state vector
    % params: a cell containing all parameters required:
        % W: process noise vector
        % nk: a vector indicating the number of Gaussian Kernels used for modeling
        % each independent component
        % fs:sampling frequency

    L=size(X,1); % number of states in the state vector X
    
    W= params{1};
    nk= params{2};
    fs= params{3};
    
    Lind = length(nk); % number of independent components used for MECG cancellation
    nParam = length(W); % number of process noise parameters
    w=W(1); % average heart-rate in rads
    df= zeros(L, nParam); % Note that L=Lind+1
    
    dt= 1/fs;
    df(1,1)=dt;
    ptr=1; % pointer to the current parameter extracted from the W vector
    for i=1:Lind
        kernelParam = reshape(W(ptr+1:ptr+nk(i)*3),3,[]); % kernel parameters for the current independent component
        dfi= feval(@df_dW, X, kernelParam, w, fs);
        df(i+1,1)= dfi(1);
        df(i+1,ptr+1:ptr+nk(i)*3)=dfi(2:end);
        df(i+1,end-Lind+i)=1;
        ptr=ptr+nk(i)*3;
    end

end

function y = df_dW(X, kp, w, fs)
    % X: state vector
    % kp: kernel parameters
    % w:  average heart-rate in rads
    % fs: sampling frequency
    % np: number of process noise parameters
    dt= 1/fs;
    alpha= kp(1,:);
    theta= kp(2,:);
    b=kp(3,:);
    
    dtheta=  rem(X(1) - theta,2*pi);
    np= numel(kp);
    partDeriv= zeros(size(kp));
    
    y=zeros(1,np+1);
    y(1)= -sum(dt*alpha.*dtheta./(b.^2).*exp(-dtheta.^2./(2*b.^2)));
    partDeriv(1,:)= -dt*w./(b.^2).*dtheta .* exp(-dtheta.^2./(2*b.^2));
    partDeriv(2,:)= 2*dt.*alpha.*w.*dtheta./b.^3.*(1 - dtheta.^2./(2*b.^2)).*exp(-dtheta.^2./(2*b.^2));
    partDeriv(3,:)= dt*w*alpha./(b.^2).*exp(-dtheta.^2./(2*b.^2)) .* (1 - dtheta.^2./b.^2);
    y(2:end)= partDeriv(:)';
    
end