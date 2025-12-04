function [out, MFRh] = QRSdetection(in, fs, mode)
% QRS detection based on Christov algorithm published in [1]
%   (c) Jakub Kuzilek, Lenka Lhotska
% 	http://bio.felk.cvut.cz/~kuziljak/     E-mail: jakub.kuzilek@gmail.com
%   Version: 1.0                    Last update:   07/01/2013.
% 					                (Version: 1.0, 07/01/2013)
%
%======================================================
%
% PURPOSE:    This is main function of QRS detection algorithm based 
%             on Christovs detection algoritm and Independent Component
%             Analysis.
%
% MANDATORY INPUT ARGUMENTS
%   in ..... input data MxN, M - length of data, N - number of leads
%   fs ..... sampling frequency in Hz
%   mode ... 0 - Christov, 1 - Christov + ICA
% OPTIONAL INPUT ARGUMENTS
%   none
% OUTPUT ARGUMENTS
%   out .... time in samples of detected peaks
%=========================================================
%
% Related Bibliography:
%
% [1] Kuzilek J, Lhotska L. Electrocardiogram beat detection 
%      enhancement using independent component analysis. Medical
%      Engineering and Physics. 2012.
%
%=========================================================
oin = in;

in = preprocessingQRSdetection(in, fs, mode); % preprocessing

N = length(in); % length of data

% help variables
MS200 = round(0.2*fs);
MS1000 = fs;
MS50 = round(0.05*fs);
MS350 = round(0.35*fs);

% constants
K1 = 0.36/MS1000; 
K2 = -0.47/1.4/MS1000; 

% init
M = 0.6*max(in(1:min([length(in),5*fs])));
MM = ones(1,5)*M;
mMM = M;

% F threshold, it does not depend on detected beats
F = zeros(1,N);
F(1:MS350) = mean(in(1:min([length(in),MS350])));

pom = zeros(MS50+1,N-MS50);
for m = 1:(N-MS50)
    pom(:,m) = in(m:m+MS50);
end
pom = max(pom);
if(mode)
    F(MS350:N-MS50) = cumsum([F(MS350),(pom(MS350+1:end)-pom(1:end-MS350))/290]); %500
else
    F(MS350:N-MS50) = cumsum([F(MS350),(pom(MS350+1:end)-pom(1:end-MS350))/500]);
end
F(end-MS50:end) = F(end-MS50);

R = 0;
RR = [];
Rm = 0;
Rm23 = 0;

C = 1;

% detected QRS and help vars
QRS = zeros(1,1000);
iQRS = 1;

RRl = 0;

% last detected QRS time
tQRS = -200;

m = 0;

% history help
MFRh = zeros(1,N);
% Mh = zeros(1,N);
% Fh = zeros(1,N);
% Rh = zeros(1,N);

while m <= N-1

    m = m+1;

    % update weights
    
    tdiff = m - tQRS;
    
    if (tQRS > 0)
        
       if(tdiff > MS200 && tdiff < MS1000+MS200) 
           M = mMM - K1*mMM*(tdiff-MS200); 
        elseif(tdiff > MS1000+MS200)
           M = 0.6*mMM;
        end
        
    end
  
    if(RRl && tdiff > Rm23 && tdiff < Rm)
        R = K2*mMM*(tdiff-Rm23);
    end
    
    MFR = (M+F(m)+R)/C;

    MFRh(m) = MFR;
%     Mh(m) = M;
%     Fh(m) = F(m);
%     Rh(m) = R;
    
    % decide if QRS found + update weights
    
    if(m > tQRS + MS200 && in(m) > ((M+F(m)+R)/C)) % if sample is larger than threshold - QRS detected
        
        [ampl, time] = max(in(m:min([m+MS200, N]))); % find peak amplitude
        dQRS = m+time; 
        
        if(0.6*ampl > 1.5*MM(end)) % update MM and reset M threshold
            MM = [MM, 1.1*MM(end)];
        else
            MM = [MM, 0.6*ampl]; 
        end
        MM = MM(2:end);
        mMM = mean(MM);
        M = mMM;
        
        if(iQRS>=2) % if first RR interval can be calculated, start R calculation 
            if(RRl && (dQRS-tQRS)<0.4*RR(end) && (dQRS-tQRS)>1.2*RR(end)) % if detected RR is small or large than previous one
                RR(end+1) = RR(end);
                RR = RR(2:end);
            elseif(RRl) % if RR is ok
                RR(end+1) = dQRS-tQRS;
                RR = RR(2:end);
            else % add RR into storage
                RR(end+1) = dQRS-tQRS;
                if length(RR)==7
                    RRl=1;
                end
            end
            R = 0;
            Rm = sum(RR)/7;
            Rm23 = Rm*0.6667;
        end
        
        % store QRS position
        QRS(iQRS) = dQRS;
        iQRS = iQRS + 1;
        tQRS = dQRS;
        
        C = 1;

    end
    
    % if QRS not detected for long period change coefficient C
    if(C == 1 && RRl && tdiff > Rm*1.5)
        m = max([tQRS 1]);
        if mode
            C = 1.3;
        else
            C = 2;
        end
    elseif(tdiff > MS1000 && C == 1)
        m = max([tQRS 1]);
        if mode
            C = 1.3;
        else
            C = 2;
        end
    end
    
end

out = QRS(1:iQRS-1); % detected QRS peaks

for m = 1:length(out)
    [~,t] = max(abs(oin(max([1,out(m)-60]):min([length(oin),out(m)+60]))));
    out(m) = out(m)+ t -61;
end

d = diff(out);
% z = find(d<rng);
out(d<0.3*fs)=[];

out(out<1)=[];