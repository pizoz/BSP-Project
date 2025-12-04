function df = ekf_df_dx(X,params)

    % X: state vector
    % params: a cell containing all parameters:  
        % W: process noise vector
        % nk: a vector indicating the number of Gaussian Kernels used for modeling
        % each independent component
        % fs:sampling frequency

    L=size(X,1); % number of states in the state vector X
    
    W= params{1};
    nk= params{2};
    fs= params{3};
    
    
    Lind = length(nk); % number of independent components used for MECG cancellation

    w=W(1); % average heart-rate in rads
    df= eye(L); % Note that L=Lind+1
    ptr=1; % pointer to the current parameter extracted from the W vector
    for i=1:Lind
        kernelParam = reshape(W(ptr+1:ptr+nk(i)*3),3,[]); % kernel parameters for the current independent component
        df(i+1,1)= feval(@df_dtheta, X, kernelParam, w, fs);
        ptr=ptr+nk(i)*3;
    end

end

function y = df_dtheta(X, kp, w, fs)
    % X: state vector
    % kp: kernel parameters
    % w:  average heart-rate in rads
    % fs: sampling frequency
    dt= 1/fs;
    alpha= kp(1,:);
    theta= kp(2,:);
    b=kp(3,:);
    
    dtheta=  rem(X(1) - theta,2*pi);
    y= -dt*sum( w*alpha./(b.^2).*(1 - dtheta.^2./b.^2).*exp(-dtheta.^2./(2*b.^2)) ); 
    
end