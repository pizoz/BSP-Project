function newecgchannel=supressmqrsuseall(allchannelsecg,thischannel,mqrs,
				    olsdata,modelorder)
#ordermodel must be an odd integer
# olsdata is created by function 'preparedata2supressmqrs'
#use all the other channels to reconstruct this channel mqrs

  #------------------------------------------------------------
  ### this piece of code is common to all channels
  samplingrate=1000;
  qrsL=floor(130*samplingrate/1000);#130ms used in 'set2zeromqrs'
  #ordermodel=3;

  ##verify olsdata is correct ----------------------------
  numchannels=size(allchannelsecg,1);
  L=length(allchannelsecg);
  halfqrsL=floor(qrsL/2);
  halfmodelorder=(modelorder-1)/2;
  numsamplesperqrs=2*halfqrsL+1;

  numqrs=length(mqrs);
  if(mqrs(1)<=halfqrsL+halfmodelorder)
    mqrs=mqrs(1,2:end);
    numqrs=numqrs-1;
  endif
  if(mqrs(end)>L-halfqrsL-halfmodelorder)
    mqrs=mqrs(1,1:end-1);
    numqrs=numqrs-1;
  endif

  [a b c]=size(olsdata);

  assert(a==numqrs*(2*halfqrsL+1));
  assert(b==modelorder);
  assert(c==numchannels);

  ##end verify -----------------------------------------


  # build arx model --------------------------------

  otherchannels=[];
  for i=1:numchannels
    if(i!=thischannel)
      otherchannels=[otherchannels,i];
    endif
  endfor


#arx model generates middle point from  'order size' segment of another
#particular channel

  middlepoint=(modelorder-1)/2+1;

  input=[];

  for i=otherchannels #use i channel to generate the other channels in
		      #the neighbourood of mqrs


#exper###############################################
####chosen the reconstructing channel based only on it ability to
####recontstruct 'thischannel'

    input=[input,olsdata(:,:,i)];

  endfor
  output=olsdata(:,middlepoint,thischannel);

  [beta, sigma, r] = ols (output, input);
  recqrsthischannel=output-r;

    # if(i==otherchannels(1))
    #   best=i;
    #   error=sum(r.^2);
    #   fbeta=beta;
    #   recqrsthischannel=output-r;
    # elseif(sum(r.^2)< error)
    #   error=sum(r.^2);
    #   best=i;
    #   fbeta=beta;
    #   recqrsthischannel=output-r;
    # endif

  #endfor

##########################################



  #   complementary_i=[];#other channels than i
  #   for j=1:numchannels
  #     if(j!=i)
  # 	complementary_i=[complementary_i,j];
  #     endif
  #   endfor

  #   input=olsdata(:,:,i);
  #   output=zeros(size(olsdata,1),numchannels-1);

  #   k=0;
  #   for j=complementary_i
  #     k+=1;
  #     output(:,k)=olsdata(:,middlepoint,j);
  #   endfor

  #   ###chose the channel that best reconstruct thischannel and the other
  #   ###two channels

  #   [beta, sigma, r] = ols (output, input);

  #   if(i==otherchannels(1))
  #     best=i;
  #     error=sum(sum(r.^2));
  #     this=find(complementary_i==thischannel);
  #     fbeta=beta(:,this);
  #     recqrsthischannel=output(:, this)-r(:,this);
  #   elseif(sum(sum(r.^2))< error)
  #     error=sum(sum(r.^2));
  #     best=i;
  #     this=find(complementary_i==thischannel);
  #     fbeta=beta(:,this);
  #     recqrsthischannel=output(:, this)-r(:,this);
  #   endif

  # endfor

  newecgchannel=allchannelsecg(thischannel,:);

  ##use  allecgsegmentsthischannel to replace mqrs segments in this channel ecg 

  for i=1:numqrs
    for j=-halfqrsL:halfqrsL

      newecgchannel(mqrs(i)+j)=allchannelsecg(thischannel,mqrs(i)+j)-\
	  recqrsthischannel((i-1)*(2*halfqrsL+1)+j+1+halfqrsL);
    endfor
  endfor

#printf("reconstruction channel is %d\n",best);


endfunction