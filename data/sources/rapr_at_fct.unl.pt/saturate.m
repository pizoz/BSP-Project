function data=saturate(data)

  numchannels=size(data,1);

  maximums=max(abs(data'));

  factor=1.2;

  limit=factor*median(maximums);

  data(data>limit)=limit;

  data(data<-limit)=-limit;

endfunction
