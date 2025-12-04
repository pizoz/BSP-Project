function y=movMedian(x,w)
%x: input vector (signal)
%w: desired window size (odd number)
%y: smoothed output signal

m=length(x);
y=zeros(1,m);
halfw=floor(w/2); 
idx1=halfw+1;%starting index of input signal with moving window of size w
winsize=1:2:w-2;% window size for indexes before idx1

for i=1:idx1-1
    halfvarw=floor(winsize(i)/2);
    temp=sort(x(i-halfvarw:i+halfvarw));
    y(i)=temp(halfvarw+1);    
end

for i=idx1:m-idx1+1
   temp=sort(x(i-halfw:i+halfw));
   y(i)=temp(halfw+1);
end

for i=m-idx1+2:m
    halfvarw=floor(winsize(m+1-i)/2);
    temp=sort(x(i-halfvarw:i+halfvarw));
    y(i)=temp(halfvarw+1);
end


