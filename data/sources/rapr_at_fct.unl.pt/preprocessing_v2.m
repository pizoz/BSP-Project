function [FilteredECG] = preprocessing_v2(ECG,fs)
#uses notch filter
#this version is a testing version for eliminating 50 Hz power line
#noise


% ---- preprocess the data ----

  FilteredECG=zeros(size(ECG));

  lmedfilt=floor(12*fs/1000);
  lmovavgkernel=8;
  movavgkernel=ones(1,floor(lmovavgkernel*fs/1000))/lmovavgkernel;
  %should it really be done before mqrs detection????


  #noch filter 50Hz
  #fs=samplingrate;
  fs2=fs/2;
  #[b,a]=pei_tseng_notch ( 50/fs2, 3/fs2 ); 
  load notchfilter50Hz.txt

  for i=1:4
    FilteredECG(i,:)=removebaseline(ECG(i,:), fs);
    FilteredECG(i,:)=medfilt1(FilteredECG(i,:),lmedfilt);
    #FilteredECG(i,:)=filtfilt(movavgkernel,1,FilteredECG(i,:));
    FilteredECG(i,:)=filtfilt(b,a,FilteredECG(i,:));
  endfor

  FilteredECG=saturate(FilteredECG);

end
