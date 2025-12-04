function Y=MissedValCheck(X)
    
    if size(X,1)>size(X,2)
        X=X';
    end

    Y=X;
    nancheck=sum(isnan(X'));
    missChannel=find(nancheck);
    
    for i=1:length(missChannel)
        Y(missChannel(i),:) = fillMissedChan(X(missChannel(i),:));
    end
    
end

function y=fillMissedChan(x)
    nsamples=length(x);
    y=x;
    nancheck=isnan(x);
    
    %finding start and end indexes of each missed data interval
    diffcheck=diff(nancheck);
    ne=find(diffcheck==-1);
    ns=find(diffcheck==1)+1;
    
    %checking special cases (first missing interval)
    N=6*(ne(1)-ns(1)+1);
    if (ns(1)==ne(1) && ns(1)>1)
        y(ns(1))=(x(ns(1)-1)+x(ns(1)+1))/2;
    else
        if ns(1)<N+1
            y(ns(1):ne(1))=0;
        else
            fwdata=x(ns(1)-N:ns(1)-1);
            bwdata=x(ne(1)+N:-1:ne(1)+1);
            estimate=fillMissed(ns(1),ne(1),fwdata,bwdata);
            y(ns(1):ne(1))=estimate;
        end
    end
    
    
    for j=2:length(ns)-1
        if ns(j)==ne(j)
            estimate=(x(ns(j)-1)+x(ns(j)+1))/2;
        else
            N=6*(ne(j)-ns(j)+1); % length of the pre and postdata sequences used for AR model
            fwdata=x(ns(j)-N:ns(j)-1);
            bwdata=x(ne(j)+N:-1:ne(j)+1);
            estimate=fillMissed(ns(j),ne(j),fwdata,bwdata);
        end
        y(ns(j):ne(j))=estimate;
    end
    
    %checking special cases (last missing interval)
    N=6*(ne(end)-ns(end)+1);
    
    if (ns(end)==ne(end)) && ne(end)<nsamples
        y(ns(end))=(x(ns(end)-1)+x(ns(end)+1))/2;
    else
        if ne(end)>nsamples-N
            y(ns(end):ne(end))=0;
        else
            fwdata=x(ns(end)-N:ns(end)-1);
            bwdata=x(ne(end)+N:-1:ne(end)+1);
            estimate=fillMissed(ns(end),ne(end),fwdata,bwdata);
            y(ns(end):ne(end))=estimate;
        end
    end
    
end