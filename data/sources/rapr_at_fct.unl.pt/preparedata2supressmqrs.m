function  olsdata=preparedata2supressmqrs(allchannelsecg,mqrs,
					  qrsL, modelorder,samplingrate)

#olsdata consists on an (numqrs*(2*halfqrsL+1),modelorder,numchannels)
#matrix with input data for arx model,
#numqrs,numchannels are defined below

  L=length(allchannelsecg);
  numchannels=size(allchannelsecg,1);
  halfqrsL=floor(qrsL/2);
  halfmodelorder=(modelorder-1)/2;

  numqrs=length(mqrs);
  if(mqrs(1)<=halfqrsL+halfmodelorder)
    mqrs=mqrs(1,2:end);
    numqrs=numqrs-1;
  endif

  if(mqrs(end)>L-halfqrsL-halfmodelorder)
    mqrs=mqrs(1,1:end-1);
    numqrs=numqrs-1;
  endif

  patches=zeros(numqrs,2*halfqrsL+modelorder,numchannels);
  for c=1:numchannels
    for i=1:numqrs
      patches(i,:,c)=allchannelsecg(c,mqrs(i)-halfqrsL-halfmodelorder:\
				    mqrs(i)+halfqrsL+halfmodelorder);
    endfor
  endfor

  #numsamplesperqrs=2*halfqrsL+1-modelorder+1;

  olsdata=zeros(numqrs*(2*halfqrsL+1),modelorder,numchannels);
  for c=1:numchannels
    for i=1:numqrs
      for j=1:2*halfqrsL+1
	olsdata((i-1)*(2*halfqrsL+1)+j,:,c)=patches(i,j:j+modelorder-1,c);
      endfor
    endfor
  endfor

endfunction