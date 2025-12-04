function [qrs threshold auxvec ]=detectqrsmultichannel2(ecg,samplingfreq)

  L=size(ecg,2);
  numchannels=size(ecg,1);

  #extra filtering->inovation for detectqrsmultichannel_2
  newmovageragefiltsize=40;#ms
  movaveragekernel=ones(1,newmovageragefiltsize)/newmovageragefiltsize;
  for i=1:numchannels
    newecg(i,:)=filtfilt(movaveragekernel,1,ecg(i,:));
  endfor



  diffsignal=[zeros(numchannels,1),diff(newecg,1,2)];

  auxvec=zeros(1,L);
  windowtime=40;#ms related to qrs duration
  windowL=floor(windowtime*samplingfreq/1000);
  halfwindow=floor(windowL/2);


  peak=[ones(1,halfwindow+1),-ones(1,halfwindow)];

  for i=halfwindow+1:L-halfwindow
    for c=1:numchannels
      auxvec(i)+=abs(diffsignal(c,i-halfwindow:i+halfwindow)*peak');
    endfor
  endfor

  basicvec=auxvec;

  #replace peaks by area on a window time

   # auxvec=zeros(1,L);
   # newl=10;
   # for i=newl+1:L-newl
   #   auxvec(i)=sum(basicvec(i-newl:i+newl));
   # endfor


 

 ##################################################
 ## FIRST APPROACH ################################
 ###################################################


  ############### parameters that must be better tuned 

  threshold=50000;

  maximumqrsinterval=1400*samplingfreq/1000;
  ##reduction of threshold value is done in a local base and not globally as in
  ## the first version

  smallwindow=floor(400*samplingfreq/1000);##minimum QRS interval
  ############### 


  qrs=[];

  emptysegments=[1,length(auxvec)];

  


  while(!isempty(emptysegments)&&(threshold>10))

    

    for i=1:length(emptysegments)/2

       qrs=searchlocalmax(threshold, auxvec, qrs,
       			 emptysegments(2*i-1), emptysegments(2*i),
       			 maximumqrsinterval, smallwindow);

      # qrs=searchlocalmaxnew(threshold, auxvec, qrs,
      # 			    emptysegments(2*i-1), emptysegments(2*i),
      # 			    smallwindow);

    endfor

    if(!isempty(qrs))
      aux=zeros(1,length(auxvec));
      aux(qrs(:))=1;

      emptysegments=findemptysegments(aux, maximumqrsinterval);
    endif

    threshold=3/4*threshold;

    #debug
    #printf("size(qrs)=%d\n",size(qrs,2));fflush(stdout);

  endwhile

  

  minimumnumberqrsaux=50;

  if(length(qrs)>minimumnumberqrsaux)

  #repeat process to
  #locate possible missing qrs using new maximumqrsinterval
  
    medianqrsperiod=median(diff(qrs));
    smallwindow=floor(min(3/4*medianqrsperiod,smallwindow));#new
    maximumqrsinterval=max(floor(5/4*medianqrsperiod),2*smallwindow+10);

    refthreshold=threshold;

    aux=zeros(1,length(auxvec));
    aux(qrs(:))=1;
    
    emptysegments=findemptysegments(aux, maximumqrsinterval);

    while(!isempty(emptysegments)&&(threshold/refthreshold>0.1))

    

      for i=1:length(emptysegments)/2

	qrs=searchlocalmax(threshold, auxvec, qrs,
			   emptysegments(2*i-1), emptysegments(2*i),
			   maximumqrsinterval, smallwindow);
      endfor

      if(!isempty(qrs))
	aux=zeros(1,length(auxvec));
	aux(qrs(:))=1;

	emptysegments=findemptysegments(aux, maximumqrsinterval);
      endif

      threshold=3/4*threshold;

      
      #printf("size(qrs)=%d\n",size(qrs,2));fflush(stdout);

    endwhile

  endif

  if(length(qrs)>1)
    assert(min(diff(qrs))>=smallwindow);
  endif


 #debug
 #printf("after first approach siz qrs is %d\n",size(qrs,2));fflush(stdout);


 ##################################################
 ## SECOND APPROACH  redefining maximumqrsinterval and smallwindow
 ###################################################

 # threshold=50000;

#  ##redefining thresholds

#  differences=diff(qrs);

#  avgdiff=mean(differences);


#  maximumqrsinterval=floor(avgdiff+120*samplingfreq/1000);

#  smallwindow=floor(avgdiff-160*samplingfreq/1000);

#  ###############################################

#  #repeating the cycle

#  qrs=[];

#  emptysegments=[1,length(auxvec)];

# #printf("length emptysegments is %d and threshold is %d\n",length(emptysegments),threshold);fflush(stdout);

  

#   while(!isempty(emptysegments)&&(threshold>10))

    

#     for i=1:length(emptysegments)/2

#       qrs=searchlocalmax(threshold, auxvec, qrs,
# 			 emptysegments(2*i-1), emptysegments(2*i),
# 			 maximumqrsinterval, smallwindow);

#     endfor

#     if(!isempty(qrs))
#       aux=zeros(1,length(auxvec));
#       aux(qrs(:))=1;

#       emptysegments=findemptysegments(aux, maximumqrsinterval);
#     endif

    

#     threshold=3/4*threshold;

#     #debug
#     #printf("size(qrs)=%d\n",size(qrs,2));fflush(stdout);

#   endwhile

  #####################################################


  if(length(qrs)>1)
    assert(min(diff(qrs))>=smallwindow);
  endif

  qrs=replacebylocmax(qrs,halfwindow,sum(ecg.^2));
    
endfunction



function res=replacebylocmax(setofpoints,halfwindow,abssignal)

  numpoints=length(setofpoints);

  res=zeros(1,numpoints);

  L=length(abssignal);


  for i=1:numpoints

    if(setofpoints(i)>halfwindow)

      if(setofpoints(i)<=L-halfwindow)

	[a order]=max(abssignal(setofpoints(i)-halfwindow:\
				setofpoints(i)+halfwindow));

	res(i)=setofpoints(i)+order-halfwindow;

      else

	[a order]=max(abssignal(setofpoints(i)-halfwindow:L));

	res(i)=setofpoints(i)+order-halfwindow;

      endif

    else

      [a order]=max(abssignal(1:setofpoints(i)+halfwindow));

      res(i)=order;
    endif

  endfor

endfunction