function peakf=PsdPeak2(X,interval,leadset)

% This function finds the frequency which has peak PSD in most
% leads. X is the input PSD matrix, and each column is for a different
% lead. The function sorts the input matrix and returns the index which
% most of the leads agree to have the max power.

peakf=-1;

if any(leadset)
    tot=sum(leadset);
    lint=length(interval);
    [Xs,Xi]=sort(X);
    Y=X(interval,:);
    [m,n]=size(Y);
    [Ys,Yi]=sort(Y);
    Yi=Yi'+interval(1)-1;
    hn=zeros(4,lint);
    hn_shiftl=zeros(4,lint);
    hn_shiftr=zeros(4,lint);
    for i=1:4
        if leadset(i)
            hn(i,:)=hist(Yi(i,end-2:end),interval);
            hn_shiftl(i,:)=[hn(i,2:end),0];
            hn_shiftr(i,:)=[0 hn(i,1:end-1)];
        end
    end
    skewed_hn=hn+hn_shiftl+hn_shiftr;
    skewed_hn=skewed_hn>0;
    test=ones(1,lint);
    for i=1:4
        if leadset(i)
            test=test.*skewed_hn(i,:);
        end
    end
    windowed_test=filter([1 1 1],1,test);
    [tm,ti]=max(windowed_test);
    if tm>0
        if sum(windowed_test==tm)==1
            peakf=interval(ti-1);
        else
            k=interval(find(windowed_test==tm)-1);
            s=zeros(1,length(k));
            for i=1:length(k);
                s(i)=sum(any(Xi(end-2:end,leadset)==k(i))+any(Xi(end-2:end,leadset)==k(i)-1)+any(Xi(end-2:end,leadset)==k(i)+1));
            end
            [sm,si]=max(s);
            if sum(s==sm)==1
                peakf=k(si);
            end
            if all(diff(k)==1)
                peakf=k(floor((length(k)+1)/2));
            end
                    
        end
    end
end

