function [fqrs ecgsupressmqrs auxvec]=detectfqrsupressnew(ecgsupressmqrs,samplingfreq)

# ecgsupressmqrs-> ecg where mqrs was suppressed using 'supressmqrsuseall'

  L=size(ecgsupressmqrs,2);
  assert(size(ecgsupressmqrs,1)==1);


  ### filtering again------------------------

  lengthmovingaverage=8;
  ecgsupressmqrs=medfilt1(ecgsupressmqrs,12);
  ecgsupressmqrs=filtfilt(ones(1,lengthmovingaverage)/lengthmovingaverage,1,
			  ecgsupressmqrs);

  ###--------------------------------------------
 

  #search for peaks in diff signal

  #diffsignal=[0,diff(ecgnmqrs)];
  diffsignalsupress=[0,diff(ecgsupressmqrs)];

  auxvec=zeros(1,L);
  windowtime=10;#ms
  windowL=floor(windowtime*samplingfreq/1000);
  halfwindow=floor(windowL/2);

  #sensible for peaks
  auxsignal=[-ones(1,halfwindow),ones(1,halfwindow+1)];
  for i=halfwindow+1:L-halfwindow
    auxvec(i)=abs(diffsignalsupress(i-halfwindow:i+halfwindow)*auxsignal');
  endfor
 

  ############### parameters that must be better tuned 
  minimumfqrsinterval=floor(320*samplingfreq/1000);#320 ms
  threshold=5000;
  maximumfqrsinterval=floor(550*samplingfreq/1000);#550ms
  minimumnumberfqrs=110;90;#after that use median fqrs period to search fqrs on mother qrs
  ############### 

  fqrs=[];
  emptysegments=[1,length(auxvec)];

  while(!isempty(emptysegments)&&(threshold>0.1))

    printf("threshold=%f length(fqrs)=%d  \r",threshold,length(fqrs));
    fflush(stdout);

    for i=1:length(emptysegments)/2

      fqrs=searchlocalmaxnew2(threshold, auxvec, fqrs,
			 emptysegments(2*i-1), emptysegments(2*i),
			 minimumfqrsinterval);

      #debug
      #assert(min(diff(fqrs))>=minimumfqrsinterval);
    endfor

    threshold=3/4*threshold;

    #recompute empty intervals
    if(!isempty(fqrs))
      aux=zeros(1,length(auxvec));
      aux(fqrs(:))=1;
      
      emptysegments=findemptysegments(aux, maximumfqrsinterval);
    endif

  endwhile
 
  if(length(fqrs)>minimumnumberfqrs)

    quantile1fqrsperiod=quantile(diff(fqrs),0.30,2);
    
    #start from beginning
    fqrs=[];
    emptysegments=[1,length(auxvec)];
    threshold=5000;


    ##thats why we start from beginning !!!!!!!!
    #from now on maximumfqrsinterval depends on quantile1fqrsperiod
    maximumfqrsinterval=9/8*quantile1fqrsperiod;#8/7;
    #maximumfqrsinterval=floor(7/6*quantilebig);
    minimumfqrsinterval=floor(max([minimumfqrsinterval,6/7*quantile1fqrsperiod]));

    while(!isempty(emptysegments)&&(threshold>0.1))

      for i=1:length(emptysegments)/2

	fqrs=searchlocalmaxnew2(threshold, auxvec, fqrs,
				emptysegments(2*i-1), emptysegments(2*i),
				minimumfqrsinterval);

        #debug
	#assert(min(diff(fqrs))>=minimumfqrsinterval);
      endfor
  
      threshold=3/4*threshold;

      #recompute empty intervals
      if(!isempty(fqrs))
	aux=zeros(1,length(auxvec));
	aux(fqrs(:))=1;
      
	emptysegments=findemptysegments(aux, maximumfqrsinterval);
      endif

    endwhile

    #debug
    #assert(min(diff(fqrs))>=minimumfqrsinterval);

  endif

  #debug
  #assert(min(diff(fqrs))>=minimumfqrsinterval);
    
endfunction




