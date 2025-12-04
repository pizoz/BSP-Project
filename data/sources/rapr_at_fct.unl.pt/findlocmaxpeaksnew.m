function res=findlocmaxpeaksnew(vector,threshold, smallwindowsize)
  #return a vector of 0s and 1s with 1s in the place of local maximums
  #smallwindowsize is the smallest distance between two maximums

  #in this version we ?????????????????????????????????????????????????????????????? 

  #printf("findlocmaxpeaks threshold=%f size(vector)=%d\n",threshold,size(vector,2));fflush(stdout);


  #debug
    if(threshold<0.2)
      2+1;
    endif

  auxdata=vector;

  auxdata(vector<threshold)=0;

  points=find(auxdata>0);

  if(!isempty(points))

    start=0;

    aaux=points(find(points>start));

    while(!isempty(aaux))

      this=aaux(1);#indice from auxdata first element >0 not yet visited 

      #set to zero all neighborous from auxdata(this) that are smaller

      value=auxdata(this);

      if (this>smallwindowsize)
	if(this+smallwindowsize<=size(auxdata,2))
	  C=auxdata(this-smallwindowsize:this+smallwindowsize);
	  C(C<=value)=0;
	  auxdata(this-smallwindowsize:this+smallwindowsize)=C;
	else
	  C=auxdata(this-smallwindowsize:end);
	  C(C<=value)=0;
	  auxdata(this-smallwindowsize:end)=C;
	endif
      elseif(this+smallwindowsize<=length(auxdata))
	C=auxdata(1:this+smallwindowsize);
	C(C<=value)=0;
	auxdata(1:this+smallwindowsize)=C;
      else
	auxdata(auxdata<=value)=0;
      endif

      auxdata(this)=value;



    
      points=find(auxdata>0);
      aaux=points(find(points>this));

      #debug
      #previouspositions=points(find(points<=this));
      #assert(diff(previouspositions)>smallwindowsize);

    endwhile

    res=auxdata;
    #res=removerepetitions(auxdata);

  else
    res=zeros(size(vector));
  endif



endfunction


