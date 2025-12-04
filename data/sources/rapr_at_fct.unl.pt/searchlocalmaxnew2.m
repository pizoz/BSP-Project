
function maximumsauxsignal=searchlocalmaxnew2(threshold, auxsignal, maximumsauxsignal,\
				   beginingsegment, endsegment, \
					     windowsmall)

#in   searchlocalmaxnew2 the new maximums are inserted between
#beginingsegment+windowsmall and endsegment-windowsmall when they
#are no extremes of original signal 

#in   searchlocalmaxnew there is no 'windowbig'
  #this function 'inserts' local maximums in the interval [beginingsegment, endsegment]
  #to the set maximumsauxsignal
  #auxsignal is the signal where we search local maximums
  #maximumsauxsignal contains the indices of local maximums
  #windowsmall is the minimum distance between two local maximums
  #beginingsegment should be the next indice after the local maximum
  #endsegment should be the last indice before the local maximum


  signalsegment=auxsignal(1,beginingsegment:endsegment);

  if (beginingsegment>1)
    signalsegment(1:windowsmall)=zeros(1,windowsmall);
  endif


  if(endsegment<length(auxsignal))
     signalsegment(end-windowsmall+1:end)=zeros(1,windowsmall);
   endif

  aux=findlocmaxpeaksnew(signalsegment, threshold,windowsmall);
  


  if(any(aux))

    positions=find(aux>0);#positions of local maximums in the sub segment

 

    #debug
    #assert(min(diff(positions))>windowsmall);



    if(length(positions)>0)

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
  
    endif

    assert(diff(maximumsauxsignal)>=windowsmall);

  else
    maximumsauxsignal=maximumsauxsignal;
  endif  

endfunction



