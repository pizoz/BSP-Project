function y=FillMissedPeaks3(inp_peaks,fs)
% N=4;
% p=1;
gap=0.6*fs;
test=diff(inp_peaks);
gapcheck=test>=gap;
gapind=find(gapcheck==1);
defaultval=max(test(~gapcheck));
meanval=mean(test(~gapcheck));
numgap=sum(gapcheck);
gapdepth=floor(test(gapcheck)/defaultval);
y=inp_peaks;

for i=2:numgap-1
    ns=1;
    ne=gapdepth(i);
    fwdata=inp_peaks(gapind(i)-min(3,gapind(i)-gapind(i-1))+1:gapind(i));
    bwdata=inp_peaks(gapind(i)+min(3,gapind(i+1)-gapind(i)):-1:gapind(i)+1);
    new_peak=fillMissed2(ns,ne,fwdata,bwdata,meanval);
    if ~isnan(new_peak)
        y=[y new_peak];
        clear new_peak
    end
end

if numgap>0
    i=1;
    ns=1;
    ne=gapdepth(1);
    si=0;
    fwdata=inp_peaks(gapind(i)-min(3,gapind(i)-si)+1:gapind(i));
    if numgap==1
        bwdata=inp_peaks(gapind(i)+3:-1:gapind(i)+1);
    else
        bwdata=inp_peaks(gapind(i)+min(3,gapind(i+1)-gapind(i)):-1:gapind(i)+1);
    end
    new_peak=fillMissed2(ns,ne,fwdata,bwdata,meanval);
    if ~isnan(new_peak)
        y=[y new_peak];
        clear new_peak
    end
end


if numgap>1
    i=numgap;
    ns=1;
    ne=gapdepth(numgap);
    si=length(inp_peaks)-gapind(i);
    fwdata=inp_peaks(gapind(i)-min(3,gapind(i)-gapind(i-1))+1:gapind(i));
    bwdata=inp_peaks(gapind(i)+min(3,si):-1:gapind(i)+1);
    new_peak=fillMissed2(ns,ne,fwdata,bwdata,meanval);

    if ~isnan(new_peak)
        y=[y new_peak];
        clear new_peak
    end
end


y=sort(y);


m1=mean(diff(y(end-2:end)));
last_interval=60000-y(end);
last_peak_notdetected=last_interval>m1;
if last_peak_notdetected
    nl=floor(last_interval/m1);
    y=[y floor(y(end)+m1*[1:nl])];
end


m2=mean(diff(y(1:3)));
first_peak_notdetected=y(1)>m2;
if first_peak_notdetected
    ns=floor(y(1)/m2);
    y=[floor(y(1)-m2*[1:ns]) y];
end




    