function Xk = ekf_cost(X,params)

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

    w=W(1); % average heart-rate in rads
    dt= 1/fs;
    
    
    Xk=zeros(L,1);
    Xk(1) = X(1) + w*dt; % theta state variable
    if(Xk(1) > pi)
        Xk(1) = Xk(1) - 2*pi;
    end
    
    ptr=1; % pointer to the current parameter extracted from the W vector
    for i=1:Lind
        kernelParam = reshape(W(ptr+1:ptr+nk(i)*3),3,[]); % kernel parameters for the current independent component
        alpha= kernelParam(1,:);
        theta= kernelParam(2,:);
        b= kernelParam(3,:);
        dtheta = rem(Xk(1) - theta, 2*pi);
        Xk(i+1) = X(i+1) - dt*sum(w*alpha./(b.^2).*dtheta.*exp(-dtheta.^2./(2*b.^2))); % s state variables
        ptr=ptr+nk(i)*3;
    end

end