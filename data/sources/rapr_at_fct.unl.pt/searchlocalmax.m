
function maximumsauxsignal=searchlocalmax(threshold, auxsignal, maximumsauxsignal,\
				   beginingsegment, endsegment,  windowbig, windowsmall)

  #this function 'inserts' local maximums in the interval [beginingsegment, endsegment]
  #to the set maximumsauxsignal
  #auxsignal is the signal where we search local maximums
  #maximumsauxsignal contains the indices of local maximums
  #windowbig is the maximum distance between two local maximums
  #windowsmall is the minimum distance between two local maximums
  #beginingsegment should be the next indice after the local maximum
  #endsegment should be the last indice before the local maximum
  #endsegment-beginingsegment must be >=windowbig

  #printf("searchlocmax");fflush(stdout);

  #extract subsignal
  signalsegment=auxsignal(1,beginingsegment:endsegment);

  #debug
  #assert(length(signalsegment)>=windowbig);#verify this together with findlocmaxpeaks

  #aux=findlocmaxpeaks(signalsegment, threshold,windowbig,  windowsmall);
  aux=findlocmaxpeaksnew(signalsegment, threshold,windowsmall);

  if(!isempty(aux))

    positions=find(aux>0);#positions of local maximums in the sub segment

    #insert new local maximums in the long signal: if it is near a local maximum
    # keep the biggest
    #first elemente from positions
    if(beginingsegment>1)
      if(positions(1)<windowsmall)
	if(signalsegment(positions(1))>auxsignal(beginingsegment-1))
	  p=find(maximumsauxsignal==beginingsegment-1);
	  if(p==1)
	    maximumsauxsignal=maximumsauxsignal(1,2:end);
	  elseif(p==length(maximumsauxsignal))
	    maximumsauxsignal=maximumsauxsignal(1,1:end-1);
	  else
	    maximumsauxsignal=[maximumsauxsignal(1,1:p-1),maximumsauxsignal(1,p+1:end)];
	  endif
	else
	  positions=positions(1,2:end);      
	endif
      endif
    endif

    #last elemente from positions
    if(endsegment<length(auxsignal))
      if(length(signalsegment)-positions(end)<windowsmall)
	if(signalsegment(positions(end))>auxsignal(endsegment+1))
	  p=find(maximumsauxsignal==endsegment+1);
	  if(p==length(maximumsauxsignal))
	    maximumsauxsignal=maximumsauxsignal(1,1:end-1);
	  elseif(p==1)
	    maximumsauxsignal=maximumsauxsignal(1,2:end);
	  else
	    maximumsauxsignal=[maximumsauxsignal(1,1:p-1),maximumsauxsignal(1,p+1:end)];
	  endif
	else
	  positions=positions(1,1:end-1);
	endif
      endif
    endif


    #debug
    #assert(min(diff(positions))>windowsmall);

    # insert positions into maximumsauxsignal
    positions=beginingsegment-1+positions;
    if(beginingsegment==1)
      maximumsauxsignal=[positions, maximumsauxsignal];
    elseif(endsegment==length(auxsignal))
      maximumsauxsignal=[maximumsauxsignal, positions];
    else
      p=find(maximumsauxsignal==beginingsegment-1);
      maximumsauxsignal=[maximumsauxsignal(1,1:p),positions,maximumsauxsignal(1,p+1:end)];
    endif

    #assert(diff(maximumsauxsignal)>windowsmall);

  # else
  #   maximumsauxsignal=maximumsauxsignal;
  endif  

endfunction



