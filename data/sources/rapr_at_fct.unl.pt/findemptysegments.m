function emptysegments=findemptysegments(signal01,windowbig)

  #printf("findemptysegments");fflush(stdout);

  aux=find(signal01>0);

  aux=[0, aux, size(signal01,2)+1];

  diffs=diff(aux);

  positions=find(diffs>windowbig);

  aux=find(signal01>0);

  aux=[aux,size(signal01,2)+1];

  endings=aux(positions);#endings of big segments of zeros

  aux=[0, aux];

  beginnings=aux(positions);#beginings of big segments of zeros

  beginnings=beginnings+1;

  endings=endings-1;

  emptysegments=[];

  for i=1:size(positions,2)

    emptysegments=[emptysegments, beginnings(i), endings(i)];

  endfor
