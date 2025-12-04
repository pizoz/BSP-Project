function [newecg baseline]=removebaseline(ecg,samplingfreq)

assert(size(ecg,1)==1);

#eliminate effect of qrs and p wave
window1=floor(samplingfreq/10);#window=100ms alternative window=200ms
baseline=medfilt1(ecg,window1);

#eliminate effect of t wave
window2=floor(samplingfreq/5);#window=200ms alternative window=600ms
baseline=medfilt1(baseline,window2);

newecg=ecg-baseline;