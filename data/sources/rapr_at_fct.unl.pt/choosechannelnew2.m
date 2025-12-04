function channel=choosechannelnew2(stats)
#stats is a matrix as in 'runlearn2rankfqrsonly' but omiting the last row
#it consists (for one record) on 4 rows each one containing:
# numfqrs statistics(diff(fqrs)) (9 elements) nummqrs  statistics(diff(mqrs))
#w is the vector of weights obtained in 'run2learnrankfqrs' based on training set-a

w=zeros(1,2);
w(1)=-1;#num fqrs is the most importante factor
w(2)=0.5;#low std is the other factor 

aux=zeros(1,4);
for i=1:4
  aux(i)=stats(i,:)*w';# this is the had hoc score for channel i
endfor

[a channel]=min(aux(:));

#if the stats for these channel are all zero replace it
if(!any(stats(channel,:)))
  b=max(aux(:));
  aux(channel)=b;

  [a channel]=min(aux(:));
endif

endfunction
