% function% [A,FilteredECG,chann,Rmaterni,delay,MQRS,FQRS,fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG)
function [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG)
% fetal_QRSAnn_est=ones(1,50);
QT_Interval = 0;

%%c'è correttore periodicità picchi fetali e correttore posizione sul picco
%% SOLO DUE TEMPLATE, LAVORO A 1000HZ, allineo sulla massima correlazione,
%% nuovo correttore periodicità

%%CORRETTA CORREZIONE PERIODICITà

%facciamo riferimento a occorrenze_finali (quelle effettive post fusione)
%invece che a quelle immaginate come somma dei due template fusi

%rispetto alla prima sottomissione variazione del metodo di identificazione
%del picco nell'intorno trovato. Al posto della semplice QRS detection
%nuovo tentativo di wavelet come per la rimozione della madre: per il feto
%si potrebbe provare ad abbassare la soglia di correlazione di fusione dei
%template?


% Template algorithm for Physionet/CinC competition 2013. This function can
% be used for events 1 and 2. Participants are free to modify any
% components of the code. However the function prototype must stay the
% same:
%
% [fetal_QRSAnn_est,QT_Interval] = physionet2013(tm,ECG) where the inputs and outputs are specified
% below.
%
% inputs:
%   ECG: 4x60000 (4 channels and 1min of signal at 1000Hz) matrix of
%   abdominal ECG channels.
%   tm : Nx1 vector of time in milliseconds
% output:
%   FQRS: FQRS markers in seconds. Each marker indicates the position of one
%   of the FQRS detected by the algorithm.
%   QT_Interval:   1x1 estimated fetal QT duration (enter NaN or 0 if you do wish to calculate)
%
%
% Author: Joachim Behar - IPMG Oxford (joachim.behar@eng.ox.ac.uk)
% Last updated: March 3, 2013 Ikaro Silva


% ---- check size of ECG ----
if size(ECG,1)>size(ECG,2)
    ECG = ECG';
end

fs      = 1000;             % sampling frequency
N       = size(ECG,2);      % number of abdominal channels
debug   = 0;                % enter debug mode?

% ---- preprocessing ----
[FilteredECG, original_d] = preprocessing(ECG,fs);
% save(['ECG250_' num2str(segnale)],'FilteredECG')
%ALGORITMO DI QRS DETECTION SU TUTTI I CANALI

start_delay=1;
numero_intervalli=3;   %numero che dia intero dividendo 60000
dim_intervallo=floor(length(FilteredECG(1,:))/numero_intervalli);
MQRS1=[];
MQRS2=[];
MQRS3=[];
MQRS4=[];
Apost1=zeros(2*numero_intervalli, 50*4);
Apost2=zeros(2*numero_intervalli, 50*4);
Apost3=zeros(2*numero_intervalli, 50*4);
Apost4=zeros(2*numero_intervalli, 50*4);
occorrenze_finali1=zeros(2,1);
occorrenze_finali2=zeros(2,1);
occorrenze_finali3=zeros(2,1);
occorrenze_finali4=zeros(2,1);
noise=[0,0,0,0];
% soglia_incl=(dim_intervallo/250)/2;
% soglia_incl=round(soglia_incl);
for di=1:numero_intervalli
    [Ttotpost1,Ttot1,MQRS,Apost,Nt1,occorrenze_finali]=PeakDetection...
        (FilteredECG(1,(di-1)*dim_intervallo+1:di*dim_intervallo), 'bior6.8', 3, 1, 2500, 1, 'h', 1, 0, 0, 0, 1, [], 10*length(FilteredECG(1,:)), 100, 50*4, 0,0.85,0.8,0.8,original_d(1,4*((di-1)*dim_intervallo)+1:4*(di*dim_intervallo)));
    if(occorrenze_finali(1)==0)
        noise(1)=1;
    end
    if(noise(1)==0)
        dimx=0;
        if(length(occorrenze_finali)==1)
            dimx=1;
        else
            r=pearsoncoeff(Apost(1,:),Apost(2,:));
            if(r>0.7) %sono due template della stessa classe
                occorrenze_finali(2)=0;
                MQRS(2,:)=zeros(1,length(MQRS(1,:)));
                Apost(2,:)=zeros(1,length(Apost(1,:)));
            end
        end
        if(di==1)
            if(dimx==1)
                occorrenze_finali=[occorrenze_finali; 0];
                MQRS=[MQRS; zeros(1,length(MQRS))];
            end
            MQRS1=MQRS;
            Apost1=Apost;
            occorrenze_finali1=occorrenze_finali;
        else
            if(dimx==1)
                %             if(sum(Apost(1,:))==0)
                %                 Apost(1,:)=Apost(2,:);
                %                 Apost(2,:)=zeros(1,length(Apost(1,:)));
                %             end
                r1=pearsoncoeff(Apost(1,:),Apost1(end-1,:));
                r2=pearsoncoeff(Apost(1,:),Apost1(end,:));
                if(isnan(r1))
                    r1=-2;
                end
                if(isnan(r2))
                    r2=-2;
                end
                if(r1>r2)
                    MQRS1=[MQRS1, [(MQRS+(di-1)*4*dim_intervallo); zeros(1,length(MQRS))]];
                    Apost1=[Apost1; Apost];
                    occorrenze_finali1=[occorrenze_finali1, [occorrenze_finali; 0]];
                else
                    MQRS1=[MQRS1, [zeros(1,length(MQRS)); (MQRS+(di-1)*4*dim_intervallo)]];
                    Apost1=[Apost1; Apost(2,:);Apost(1,:)];
                    occorrenze_finali1=[occorrenze_finali1, [0; occorrenze_finali]];
                end
            else
                %correlazione template 1 con il precedente 1 e con il
                %precedente 2
                r(1)=pearsoncoeff(Apost(1,:),Apost1(end-1,:));
                r(2)=pearsoncoeff(Apost(1,:),Apost1(end,:));
                %correlazione del template 2 con il precedente 1 e con il
                %precedente 2
                r(3)=pearsoncoeff(Apost(2,:),Apost1(end-1,:));
                r(4)=pearsoncoeff(Apost(2,:),Apost1(end,:));
                for uh=1:4
                    if(isnan(r(uh)))
                        r(uh)=-2;
                    end
                end
                [masr,indr]=max(r);
                if((indr==2)||(indr==3))
                    tempo=MQRS(1,:);
                    MQRS(1,:)=MQRS(2,:);
                    MQRS(2,:)=tempo;
                    tempo=Apost(1,:);
                    Apost(1,:)=Apost(2,:);
                    Apost(2,:)=tempo;
                    tempo=occorrenze_finali(1);
                    occorrenze_finali(1)=occorrenze_finali(2);
                    occorrenze_finali(2)=tempo;
                end
                MQRS1=[MQRS1, (MQRS+(di-1)*4*dim_intervallo)];
                Apost1=[Apost1; Apost];
                occorrenze_finali1=[occorrenze_finali1, occorrenze_finali];
            end
        end
    end
%     figure
%     plot(Apost(1,:))
%     if(length(occorrenze_finali)>1)
%     hold on
%     plot(Apost(2,:),'r')
%     title(['1 template 1 blu, occ: ' num2str(occorrenze_finali1(1)) 'battiti media: '  num2str(Ttotpost1(1)) 'template 2 rosso,occ:' num2str(occorrenze_finali1(2)) 'battiti media: ' num2str(Ttotpost1(2))])
%     end
%     
    [Ttotpost2,Ttot2,MQRS,Apost,Nt2,occorrenze_finali]=PeakDetection...
        (FilteredECG(2,(di-1)*dim_intervallo+1:di*dim_intervallo), 'bior6.8', 3, 1, 2500, 1, 'h', 1, 0, 0, 0, 1, [], 10*length(FilteredECG(2,:)), 100, 50*4, 0,0.85,0.8,0.8,original_d(2,4*((di-1)*dim_intervallo)+1:4*(di*dim_intervallo)));
    if(occorrenze_finali(1)==0)
        noise(2)=1;
    end
    if(noise(2)==0)
        dimx=0;
        if(length(occorrenze_finali)==1)
            dimx=1;
        else
            r=pearsoncoeff(Apost(1,:),Apost(2,:));
            if(r>0.7) %sono due template della stessa classe
                occorrenze_finali(2)=0;
                MQRS(2,:)=zeros(1,length(MQRS(1,:)));
                Apost(2,:)=zeros(1,length(Apost(1,:)));
            end
        end
        if(di==1)
            if(dimx==1)
                occorrenze_finali=[occorrenze_finali; 0];
                MQRS=[MQRS; zeros(1,length(MQRS))];
            end
            MQRS2=MQRS;
            Apost2=Apost;
            occorrenze_finali2=occorrenze_finali;
        else
            if(dimx==1)
                %             if(sum(Apost(1,:))==0)
                %                 Apost(1,:)=Apost(2,:);
                %                 Apost(2,:)=zeros(1,length(Apost(1,:)));
                %             end
                r1=pearsoncoeff(Apost(1,:),Apost2(end-1,:));
                r2=pearsoncoeff(Apost(1,:),Apost2(end,:));
                if(isnan(r1))
                    r1=-2;
                end
                if(isnan(r2))
                    r2=-2;
                end
                if(r1>r2)
                    MQRS2=[MQRS2, [(MQRS+(di-1)*4*dim_intervallo); zeros(1,length(MQRS))]];
                    Apost2=[Apost2; Apost];
                    occorrenze_finali2=[occorrenze_finali2, [occorrenze_finali; 0]];
                else
                    MQRS2=[MQRS2, [zeros(1,length(MQRS)); (MQRS+(di-1)*4*dim_intervallo)]];
                    Apost2=[Apost2; Apost(2,:);Apost(1,:)];
                    occorrenze_finali2=[occorrenze_finali2, [0; occorrenze_finali]];
                end
            else
                %correlazione template 1 con il precedente 1 e con il
                %precedente 2
                r(1)=pearsoncoeff(Apost(1,:),Apost2(end-1,:));
                r(2)=pearsoncoeff(Apost(1,:),Apost2(end,:));
                %correlazione del template 2 con il precedente 1 e con il
                %precedente 2
                r(3)=pearsoncoeff(Apost(2,:),Apost2(end-1,:));
                r(4)=pearsoncoeff(Apost(2,:),Apost2(end,:));
                for uh=1:4
                    if(isnan(r(uh)))
                        r(uh)=-2;
                    end
                end
                [masr,indr]=max(r);
                if((indr==2)||(indr==3))
                    tempo=MQRS(1,:);
                    MQRS(1,:)=MQRS(2,:);
                    MQRS(2,:)=tempo;
                    tempo=Apost(1,:);
                    Apost(1,:)=Apost(2,:);
                    Apost(2,:)=tempo;
                    tempo=occorrenze_finali(1);
                    occorrenze_finali(1)=occorrenze_finali(2);
                    occorrenze_finali(2)=tempo;
                end
                MQRS2=[MQRS2, (MQRS+(di-1)*4*dim_intervallo)];
                Apost2=[Apost2; Apost];
                occorrenze_finali2=[occorrenze_finali2, occorrenze_finali];
            end
        end
    end
%     figure
%     plot(Apost(1,:))
%     if(length(occorrenze_finali)>1)
%     hold on
%     plot(Apost(2,:),'r')
%     title(['2 template 1 blu, occ: ' num2str(occorrenze_finali2(1)) 'battiti media: '  num2str(Ttotpost2(1)) 'template 2 rosso,occ:' num2str(occorrenze_finali2(2)) 'battiti media: ' num2str(Ttotpost2(2))])
%     end
    
    [Ttotpost3,Ttot3,MQRS,Apost,Nt3,occorrenze_finali]=PeakDetection...
        (FilteredECG(3,(di-1)*dim_intervallo+1:di*dim_intervallo), 'bior6.8', 3, 1, 2500, 1, 'h', 1, 0, 0, 0, 1, [], 10*length(FilteredECG(3,:)), 100, 50*4, 0,0.85,0.8,0.8,original_d(3,4*((di-1)*dim_intervallo)+1:4*(di*dim_intervallo)));
    if(occorrenze_finali(1)==0)
        noise(3)=1;
    end
    if(noise(3)==0)
        dimx=0;
        if(length(occorrenze_finali)==1)
            dimx=1;
        else
            r=pearsoncoeff(Apost(1,:),Apost(2,:));
            if(r>0.7) %sono due template della stessa classe
                occorrenze_finali(2)=0;
                MQRS(2,:)=zeros(1,length(MQRS(1,:)));
                Apost(2,:)=zeros(1,length(Apost(1,:)));
            end
        end
        if(di==1)
            if(dimx==1)
                occorrenze_finali=[occorrenze_finali; 0];
                MQRS=[MQRS; zeros(1,length(MQRS))];
            end
            MQRS3=MQRS;
            Apost3=Apost;
            occorrenze_finali3=occorrenze_finali;
        else
            if(dimx==1)
                %             if(sum(Apost(1,:))==0)
                %                 Apost(1,:)=Apost(2,:);
                %                 Apost(2,:)=zeros(1,length(Apost(1,:)));
                %             end
                r1=pearsoncoeff(Apost(1,:),Apost3(end-1,:));
                r2=pearsoncoeff(Apost(1,:),Apost3(end,:));
                if(isnan(r1))
                    r1=-2;
                end
                if(isnan(r2))
                    r2=-2;
                end
                if(r1>r2)
                    MQRS3=[MQRS3, [(MQRS+(di-1)*4*dim_intervallo); zeros(1,length(MQRS))]];
                    Apost3=[Apost3; Apost];
                    occorrenze_finali3=[occorrenze_finali3, [occorrenze_finali; 0]];
                else
                    MQRS3=[MQRS3, [zeros(1,length(MQRS)); (MQRS+(di-1)*4*dim_intervallo)]];
                    Apost3=[Apost3; Apost(2,:);Apost(1,:)];
                    occorrenze_finali3=[occorrenze_finali3, [0; occorrenze_finali]];
                end
            else
                %correlazione template 1 con il precedente 1 e con il
                %precedente 2
                r(1)=pearsoncoeff(Apost(1,:),Apost3(end-1,:));
                r(2)=pearsoncoeff(Apost(1,:),Apost3(end,:));
                %correlazione del template 2 con il precedente 1 e con il
                %precedente 2
                r(3)=pearsoncoeff(Apost(2,:),Apost3(end-1,:));
                r(4)=pearsoncoeff(Apost(2,:),Apost3(end,:));
                for uh=1:4
                    if(isnan(r(uh)))
                        r(uh)=-2;
                    end
                end
                [masr,indr]=max(r);
                if((indr==2)||(indr==3))
                    tempo=MQRS(1,:);
                    MQRS(1,:)=MQRS(2,:);
                    MQRS(2,:)=tempo;
                    tempo=Apost(1,:);
                    Apost(1,:)=Apost(2,:);
                    Apost(2,:)=tempo;
                    tempo=occorrenze_finali(1);
                    occorrenze_finali(1)=occorrenze_finali(2);
                    occorrenze_finali(2)=tempo;
                end
                MQRS3=[MQRS3, (MQRS+(di-1)*4*dim_intervallo)];
                Apost3=[Apost3; Apost];
                occorrenze_finali3=[occorrenze_finali3, occorrenze_finali];
            end
        end
    end
%     figure
%     plot(Apost(1,:))
%     if(length(occorrenze_finali)>1)
%     hold on
%     plot(Apost(2,:),'r')
%     title(['3 template 1 blu, occ: ' num2str(occorrenze_finali3(1)) 'battiti media: '  num2str(Ttotpost3(1)) 'template 2 rosso,occ:' num2str(occorrenze_finali3(2)) 'battiti media: ' num2str(Ttotpost3(2))])
% 
%     end
    
    [Ttotpost4,Ttot4,MQRS,Apost,Nt4,occorrenze_finali]=PeakDetection...
        (FilteredECG(4,(di-1)*dim_intervallo+1:di*dim_intervallo), 'bior6.8', 3, 1, 2500, 1, 'h', 1, 0, 0, 0, 1, [], 10*length(FilteredECG(4,:)), 100, 50*4, 0,0.85,0.8,0.8,original_d(4,4*((di-1)*dim_intervallo)+1:4*(di*dim_intervallo)));
    if(occorrenze_finali(1)==0)
        noise(4)=1;
    end
    if(noise(4)==0)
        dimx=0;
        if(length(occorrenze_finali)==1)
            dimx=1;
        else
            r=pearsoncoeff(Apost(1,:),Apost(2,:));
            if(r>0.7) %sono due template della stessa classe
                occorrenze_finali(2)=0;
                MQRS(2,:)=zeros(1,length(MQRS(1,:)));
                Apost(2,:)=zeros(1,length(Apost(1,:)));
            end
        end
        if(di==1)
            if(dimx==1)
                occorrenze_finali=[occorrenze_finali; 0];
                MQRS=[MQRS; zeros(1,length(MQRS))];
            end
            MQRS4=MQRS;
            Apost4=Apost;
            occorrenze_finali4=occorrenze_finali;
        else
            if(dimx==1)
                %             if(sum(Apost(1,:))==0)
                %                 Apost(1,:)=Apost(2,:);
                %                 Apost(2,:)=zeros(1,length(Apost(1,:)));
                %             end
                r1=pearsoncoeff(Apost(1,:),Apost4(end-1,:));
                r2=pearsoncoeff(Apost(1,:),Apost4(end,:));
                if(isnan(r1))
                    r1=-2;
                end
                if(isnan(r2))
                    r2=-2;
                end
                if(r1>r2)
                    MQRS4=[MQRS4, [(MQRS+(di-1)*4*dim_intervallo); zeros(1,length(MQRS))]];
                    Apost4=[Apost4; Apost];
                    occorrenze_finali4=[occorrenze_finali4, [occorrenze_finali; 0]];
                else
                    MQRS4=[MQRS4, [zeros(1,length(MQRS)); (MQRS+(di-1)*4*dim_intervallo)]];
                    Apost4=[Apost4; Apost(2,:);Apost(1,:)];
                    occorrenze_finali4=[occorrenze_finali4, [0; occorrenze_finali]];
                end
            else
                %correlazione template 1 con il precedente 1 e con il
                %precedente 2
                r(1)=pearsoncoeff(Apost(1,:),Apost4(end-1,:));
                r(2)=pearsoncoeff(Apost(1,:),Apost4(end,:));
                %correlazione del template 2 con il precedente 1 e con il
                %precedente 2
                r(3)=pearsoncoeff(Apost(2,:),Apost4(end-1,:));
                r(4)=pearsoncoeff(Apost(2,:),Apost4(end,:));
                for uh=1:4
                    if(isnan(r(uh)))
                        r(uh)=-2;
                    end
                end
                [masr,indr]=max(r);
                if((indr==2)||(indr==3))
                    tempo=MQRS(1,:);
                    MQRS(1,:)=MQRS(2,:);
                    MQRS(2,:)=tempo;
                    tempo=Apost(1,:);
                    Apost(1,:)=Apost(2,:);
                    Apost(2,:)=tempo;
                    tempo=occorrenze_finali(1);
                    occorrenze_finali(1)=occorrenze_finali(2);
                    occorrenze_finali(2)=tempo;
                end
                MQRS4=[MQRS4, (MQRS+(di-1)*4*dim_intervallo)];
                Apost4=[Apost4; Apost];
                occorrenze_finali4=[occorrenze_finali4, occorrenze_finali];
            end
        end
    end
%     figure
%     plot(Apost(1,:))
%     if(length(occorrenze_finali)>1)
%     hold on
%     plot(Apost(2,:),'r')
%     title(['4 template 1 blu, occ: ' num2str(occorrenze_finali4(1)) 'battiti media: '  num2str(Ttotpost4(1)) 'template 2 rosso,occ:' num2str(occorrenze_finali4(2)) 'battiti media: ' num2str(Ttotpost4(2))])
%     end
end

occorrenze_finali1=sum(occorrenze_finali1,2);
occorrenze_finali2=sum(occorrenze_finali2,2);
occorrenze_finali3=sum(occorrenze_finali3,2);
occorrenze_finali4=sum(occorrenze_finali4,2);
% figure
% h=original_d(1,:);
% plot(original_d(1,:))
% hold on
% plot(MQRS1(1,:),h(MQRS1(1,:)),'*r')


%%


%%%valutazione canali rumorosi 
th_rum=20;   %soglia per individuare canali troppo rumorosi con il template matching


fin_occ=[max(occorrenze_finali1), max(occorrenze_finali2), max(occorrenze_finali3), max(occorrenze_finali4)];
% for ly=1:4
%     if fin_occ(ly)<th_rum;
%         fin_occ(ly)=0;
%     end
% end

canali_no_rumore=find(fin_occ>0);
canali_rumore=find(fin_occ==0);
if(isempty(canali_rumore))
    canali_rumore=[0];
end
if(noise(1)==0)
    temp1=[Apost1(1,:);Apost1(3,:); Apost1(5,:)];
    temp2=[Apost1(2,:);Apost1(4,:); Apost1(6,:)];
    Apost1(1:3,:)=temp1;
    Apost1(4:6,:)=temp2;
end


if(noise(2)==0)
temp1=[Apost2(1,:);Apost2(3,:); Apost2(5,:)];
temp2=[Apost2(2,:);Apost2(4,:); Apost2(6,:)];
Apost2(1:3,:)=temp1;
Apost2(4:6,:)=temp2;
end

if(noise(3)==0)
temp1=[Apost3(1,:);Apost3(3,:); Apost3(5,:)];
temp2=[Apost3(2,:);Apost3(4,:); Apost3(6,:)];
Apost3(1:3,:)=temp1;
Apost3(4:6,:)=temp2;
end

if(noise(4)==0)
temp1=[Apost4(1,:);Apost4(3,:); Apost4(5,:)];
temp2=[Apost4(2,:);Apost4(4,:); Apost4(6,:)];
Apost4(1:3,:)=temp1;
Apost4(4:6,:)=temp2;
end
% figure
% plot(Apost4(1,:))
% hold on
% plot(Apost4(2,:))
% plot(Apost4(3,:))
% figure
% plot(Apost4(4,:))
% hold on
% plot(Apost4(5,:))
% plot(Apost4(6,:))


%% QRS DETECTOR PER TROVARE POSIZIONE BATTITI MATERNI %250Hz

% [qrs1]= balda(FilteredECG(1,:),1);   
[qrs1]= balda(FilteredECG(1,:),0);  
[qrs2]= balda(FilteredECG(2,:),0);
[qrs3]= balda(FilteredECG(3,:),0);
[qrs4]= balda(FilteredECG(4,:),0);
y=FilteredECG(1,:);
   
%% creazione sequenze di classi
[mark1] = match(qrs1*4, MQRS1, 15*4);
[mark2] = match(qrs2*4, MQRS2, 15*4);
[mark3] = match(qrs3*4, MQRS3, 15*4);
[mark4] = match(qrs4*4, MQRS4, 15*4);

%% ANALISI DELLA PERIODICITà PER SCEGLIERE IL CANALE MATERNO
th_periodicity=40;
forget=1/4;
delay=0;
% [occtmp1,occchk1,RRm1,count_removed1,count_added1,noise1]=periodicity_correction_v2(qrs1,length(y),th_periodicity,forget, delay);
% [occtmp2,occchk2,RRm2,count_removed2,count_added2,noise2]=periodicity_correction_v2(qrs2,length(y),th_periodicity,forget, delay);
% [occtmp3,occchk3,RRm3,count_removed3,count_added3,noise3]=periodicity_correction_v2(qrs3,length(y),th_periodicity,forget, delay);
% [occtmp4,occchk4,RRm4,count_removed4,count_added4,noise4]=periodicity_correction_v2(qrs4,length(y),th_periodicity,forget, delay);
[occtmp1,occchk1,RRm1,count_removed1,count_added1,noise1]=periodicity_cor_caller...
    (FilteredECG(1,:),mark1,qrs1,length(y),th_periodicity,forget, delay);
[occtmp2,occchk2,RRm2,count_removed2,count_added2,noise2]=periodicity_cor_caller...
    (FilteredECG(2,:),mark2,qrs2,length(y),th_periodicity,forget, delay);
[occtmp3,occchk3,RRm3,count_removed3,count_added3,noise3]=periodicity_cor_caller...
    (FilteredECG(3,:),mark3,qrs3,length(y),th_periodicity,forget, delay);
[occtmp4,occchk4,RRm4,count_removed4,count_added4,noise4]=periodicity_cor_caller...
    (FilteredECG(4,:),mark4,qrs4,length(y),th_periodicity,forget, delay);

RRm=[RRm1, RRm2, RRm3, RRm4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTA
% per la scelta del canale io farei così: per quelli che non sono noise,
% cercherei quelli con RRm simile e alto (madre), poi sceglierei sulla base
% del numero di correzioni fatte. Le soglie che metti sotto sono troppo
% basse. Inoltre metterei un miglior QRS detector (magari con un matched
% filtro un po' più largo, perché il balda che fai mi pare un po' una
% porcheria (gli piace molto più il feto che la mamma. Il problema del
% template sul segnale 3 io non riesco a vederlo proprio. Inoltre ho notato
% che su un segnale i due template sono entrambi materni e anche molto
% simili. Potremmo provare ad abbassare un po' la soglia di correlazione
% (non molto) perché è strano che non si fondano. Inoltre per quello che
% vedo io non dovrebbe essere molto difficile capire chi è feto e chi
% madre...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCELTA CANALE MATERNO
cont_tot=zeros(4,1);
cont_tot(1)=count_removed1+1.5*count_added1;  %peso di più gli aggiunti
cont_tot(2)=count_removed2+1.5*count_added2;
cont_tot(3)=count_removed3+1.5*count_added3;
cont_tot(4)=count_removed4+1.5*count_added4;

numero_picchi=[length(occtmp1);length(occtmp2);length(occtmp3);length(occtmp4)];
etichetta_picchi=[1,2,3,4];
u=1;
while u<4
    if((numero_picchi(u)<45)||cont_tot(u)>50)
        numero_picchi(u)=[];
        etichetta_picchi(u)=[];
        cont_tot(u)=[];
    else
        u=u+1;
    end
    if(u>length(numero_picchi))
        break;
    end

end

if(isempty(numero_picchi))
    madre=1;    
elseif(length(numero_picchi)==1)
    madre=etichetta_picchi(1);
else
    media_picchi=mean(numero_picchi);
    test=abs(numero_picchi-media_picchi*ones(length(numero_picchi),1));
    massimo_picchi=max(test);
    due_classi=0;
    if(massimo_picchi>5)
        due_classi=1;
    end
    if(due_classi==0)
        [minimo,ind_min]=min(cont_tot);
        madre=etichetta_picchi(ind_min);
    else
        classe1=etichetta_picchi(1);
        num_classe1=numero_picchi(1);
        classe2=[];
        num_classe2=[];
        for zj=2:length(numero_picchi)
            if(abs(num_classe1(1)-numero_picchi(zj))<5)
                classe1=[classe1, etichetta_picchi(zj)];
                num_classe1=[num_classe1, numero_picchi(zj)];
            else
                classe2=[classe2, etichetta_picchi(zj)];
                num_classe2=[num_classe2, numero_picchi(zj)];
            end
        end
        media1=mean(num_classe1);
        media2=mean(num_classe2);
        if(media1<media2)
            errori_classe1=[];
            for lk=1:length(classe1)
                errori_classe1=[errori_classe1, cont_tot(classe1(lk))];
                [min_err, ind_min_err]=min(errori_classe1);
                madre=classe1(ind_min_err);
            end
        else
            errori_classe2=[];
            for lk=1:length(classe2)
                errori_classe2=[errori_classe2, cont_tot(classe2(lk))];
                [min_err, ind_min_err]=min(errori_classe2);
                madre=classe2(ind_min_err);
            end
            
        end
        
    end
end
    
RRmM=RRm(madre);   


% % cont_tot
% minimo=100000;
% for u=1:4
%     if cont_tot(u)<10
%         if((numero_picchi(u)<=minimo)&&(numero_picchi(u)>45))
%             madre=u;
%             minimo=numero_picchi(u);
%         end
%     end
% end
            
% [val,madre]=min(cont_tot);
    switch(madre)
        case 1
            battiti_materni=occtmp1;
        case 2
            battiti_materni=occtmp2;
        case 3
            battiti_materni=occtmp3;
        case 4
            battiti_materni=occtmp4;
    end 

% figure
% subplot(4,1,1)
% y=FilteredECG(1,:);
% plot(y)
% hold on
% plot(battiti_materni,y(battiti_materni),'*r')
% title(['canale scelto ' num2str(madre)]) 
% subplot(4,1,2)
% y=FilteredECG(2,:);
% plot(y)
% hold on
% plot(battiti_materni,y(battiti_materni),'*r')
% subplot(4,1,3)
% y=FilteredECG(3,:);
% plot(y)
% hold on
% plot(battiti_materni,y(battiti_materni),'*r')
% subplot(4,1,4)
% y=FilteredECG(4,:);
% plot(y)
% hold on
% plot(battiti_materni,y(battiti_materni),'*r')


   
%%inserisci formula per pesare meno count_removed rispetto a count_added?

%% SOTTRAZIONE DEL CANALE MATERNO
%% FASE 1. SCELTA DEL TEMPLATE MATERNO IN OGNI CANALE


battiti_materni=4*battiti_materni;  %riferiti a 1000Hz


[r,c]=size(MQRS1);
if (r>1)
    [template1]=scelta_template(battiti_materni, MQRS1);
else
    template1=1;
end
[r,c]=size(MQRS2);
if (r>1)
    [template2]=scelta_template(battiti_materni, MQRS2);
else
    template2=1;
end
[r,c]=size(MQRS3);
if (r>1)
    [template3]=scelta_template(battiti_materni, MQRS3);
else
    template3=1;
end
[r,c]=size(MQRS4);
if (r>1)
    [template4]=scelta_template(battiti_materni, MQRS4);
else
    template4=1;
end
 


%% FASE 3.SOTTRAZIONE TEMPLATE MEDIO DI OGNI CANALE
% ---- channel selection or combination prior FQRS detection ----
% [chann] = ChannelSelectionOrCombination(occorrenze_finali1,occorrenze_finali2, occorrenze_finali3,occorrenze_finali4);   %non dovrebbe usare num_occ ma occorrenze_finali che sono quelle VERE dopo la fusione e la riduzione dei template. Attenzione che Apost è ordinato secondo num_occ

% ---- MECG cancellation ----
% for i=1:N               % run algorithm for each channel
%     FECG(:,i) = MECGcancellation(MQRS,FilteredECG(:,i)',fs,20);
% end
if(noise(1)==0)
[FECG1]=MECGcancellation(original_d(1,:),battiti_materni, Apost1, template1);
else
    FECG1=(original_d(1,:));
end
if(noise(2)==0)
[FECG2]=MECGcancellation(original_d(2,:),battiti_materni, Apost2, template2);
else
    FECG2=(original_d(2,:));
end
if(noise(3)==0)
[FECG3]=MECGcancellation(original_d(3,:),battiti_materni, Apost3, template3);
else
FECG3=(original_d(3,:));
end
if(noise(4)==0)
[FECG4]=MECGcancellation(original_d(4,:),battiti_materni, Apost4, template4);
else
    FECG4=(original_d(4,:));
end
FECG=[FECG1; FECG2; FECG3; FECG4];
for tk=1:4
   temp=resample(FECG(tk,:)',250,1000);
   FECG_down(tk,:)=temp';
end
% save FECG.mat
% save(['FECG' num2str(segnale)],'FECG')
% finestra_feto=24;

FECG1_down=resample(FECG1',250,1000);
FECG1_down=FECG1_down';

if(noise(1)==0)
[Ttotpost1,Ttot1,MQRS1,Apost1,Nt1,occorrenze_finali1]=PeakDetection...
        (FECG_down(1,:), 'bior6.8', 3, 1, 2500, 1, 'h', 1, 0, 0, 0, 1, [], 10*length(FilteredECG(1,:)), 100, 200, 0,0.8,0.8,0.75,FECG1);
end
if(noise(2)==0)
[Ttotpost2,Ttot2,MQRS2,Apost2,Nt2,occorrenze_finali2]=PeakDetection...
        (FECG_down(2,:), 'bior6.8', 3, 1, 2500, 1, 'h', 1, 0, 0, 0, 1, [], 10*length(FilteredECG(2,:)), 100, 200, 0,0.8,0.8,0.75,FECG2);
end
if(noise(3)==0)
[Ttotpost3,Ttot3,MQRS3,Apost3,Nt3,occorrenze_finali3]=PeakDetection...
        (FECG_down(3,:), 'bior6.8', 3, 1, 2500, 1, 'h', 1, 0, 0, 0, 1, [], 10*length(FilteredECG(3,:)), 100, 200, 0,0.8,0.8,0.75,FECG3);
end
if(noise(4)==0)
[Ttotpost4,Ttot4,MQRS4,Apost4,Nt41,occorrenze_finali4]=PeakDetection...
        (FECG_down(4,:), 'bior6.8', 3, 1, 2500, 1, 'h', 1, 0, 0, 0, 1, [], 10*length(FilteredECG(4,:)), 100, 200, 0,0.8,0.8,0.75,FECG4);
end
    
 
[qrs1]=balda(FECG_down(1,:),1);
[qrs2]=balda(FECG_down(2,:),1);
[qrs3]=balda(FECG_down(3,:),1);
[qrs4]=balda(FECG_down(4,:),1);





%% periodicity_cor_callerFETUS;

%% creazione sequenze di classi
[mark1] = match(qrs1*4, MQRS1, 15*4);
[mark2] = match(qrs2*4, MQRS2, 15*4);
[mark3] = match(qrs3*4, MQRS3, 15*4);
[mark4] = match(qrs4*4, MQRS4, 15*4);

%% ANALISI DELLA PERIODICITà PER SCEGLIERE IL CANALE MATERNO
th_periodicity=60;
forget=1/4;
delay=0;
% [occtmp1,occchk1,RRm1,count_removed1,count_added1,noise1]=periodicity_correction_v2(qrs1,length(y),th_periodicity,forget, delay);
% [occtmp2,occchk2,RRm2,count_removed2,count_added2,noise2]=periodicity_correction_v2(qrs2,length(y),th_periodicity,forget, delay);
% [occtmp3,occchk3,RRm3,count_removed3,count_added3,noise3]=periodicity_correction_v2(qrs3,length(y),th_periodicity,forget, delay);
% [occtmp4,occchk4,RRm4,count_removed4,count_added4,noise4]=periodicity_correction_v2(qrs4,length(y),th_periodicity,forget, delay);
[occtmpf1,occchk1,RRmf1,count_removed1,count_added1,noise1]=periodicity_cor_callerFETUS...
    (FECG_down(1,:),mark1,qrs1,length(y),th_periodicity,forget, delay);
[occtmpf2,occchk2,RRmf2,count_removed2,count_added2,noise2]=periodicity_cor_callerFETUS...
    (FECG_down(2,:),mark2,qrs2,length(y),th_periodicity,forget, delay);
[occtmpf3,occchk3,RRmf3,count_removed3,count_added3,noise3]=periodicity_cor_callerFETUS...
    (FECG_down(3,:),mark3,qrs3,length(y),th_periodicity,forget, delay);
[occtmpf4,occchk4,RRmf4,count_removed4,count_added4,noise4]=periodicity_cor_callerFETUS...
    (FECG_down(4,:),mark4,qrs4,length(y),th_periodicity,forget, delay);

RRmf=[RRmf1, RRmf2, RRmf3, RRmf4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cont_tot=zeros(4,1);
cont_tot(1)=count_removed1+1.5*count_added1;  %peso di più gli aggiunti
cont_tot(2)=count_removed2+1.5*count_added2;
cont_tot(3)=count_removed3+1.5*count_added3;
cont_tot(4)=count_removed4+1.5*count_added4;
minimo=cont_tot(1);
feto=1;
fh=1;
scelgo=cont_tot;
while(fh>0)
    [minimo,feto]=min(cont_tot);
    if(abs(RRmf(feto)-RRmM)<5)
        canale_materno=0;
        switch (feto)
            case 1
                [canale_materno]=controlla_scelta(battiti_materni,occtmpf1);
            case 2
                [canale_materno]=controlla_scelta(battiti_materni,occtmpf2);
            case 3
                [canale_materno]=controlla_scelta(battiti_materni,occtmpf3);
            case 4
                [canale_materno]=controlla_scelta(battiti_materni,occtmpf4);
        end
        %becca residuo materno periodico
        if(canale_materno==1)
            cont_tot(feto)=10000;
        else
            fh=0;
        end
    else
        fh=0;
    end
    if(minimo==10000)
        [minimo, feto]=min(scelgo);
        fh=0;
    end
end

    switch(feto)
        case 1
            battiti_fetali=occtmpf1;
            battito_medio=Apost1(1,:);
        case 2
            battiti_fetali=occtmpf2;
            battito_medio=Apost2(1,:);
        case 3
            battiti_fetali=occtmpf3;
            battito_medio=Apost3(1,:);
        case 4
            battiti_fetali=occtmpf4;
            battito_medio=Apost4(1,:);
    end 

    

%%

if(debug==0)  




%% CONTROLLO FREQUENZA MATERNA E FETALE PER EVENTUALE SCAMBIO

% total_delay=delay;
% FQRS=pic_f;
% FQRS=MQRSf(1,:);
FQRS=4*battiti_fetali;
% total_delay=0;  %ritardo balda
%% CORREZIONE FETI
% FQRS=FQRS-total_delay;
% FQRS=4*FQRS;
% upsam=resample(FilteredECG(1,:)',1000,250);  
% upsam=upsam';

if(max(battito_medio(1,:))==max(abs(battito_medio(1,:))))
    %morfologia positiva
    morf=1;
else
    morf=0;
end
[FQRS]=correggi(FQRS,FECG(feto,:),morf);

% FQRS=corretti;

% pre_FQRS=FQRS;
% pre_FQRS=pre_FQRS;
% FQRS=find(FQRS>1);

% pre_FQRS=pre_FQRS(frs(1):end);
% [FQRS]=correction(upsam,FQRS); %NEL CORREGGERE PROVARE A USARE IL SEGNALE ORIGINALE, NON QUELLO PULITO PER VEDERE SE MIGLIORA
% FQRS=FQRS;  %@1000Hz
z=find(FQRS>1);
FQRS=FQRS(z(1):end);
fetal_QRSAnn_est= round(1000*FQRS'/fs);
QT_Interval = 0;
end
end

function[canale_materno]=controlla_scelta(battiti_materni, feto)
    qrs_1=feto;
    canale_materno=0;
    soglia_val=20;
    zi=1;
    contatore_1=0;
    while (zi<length(battiti_materni))
        zj=1;
        while(zj<length(qrs_1))
            if(abs(qrs_1(zj)- battiti_materni(zi))<=soglia_val)
                contatore_1= contatore_1 + 1;
                qrs_1(zj)=[];
                break;
            end
            zj=zj+1;
        end 
        zi=zi+1;
    end
    if(contatore_1>=0.8*length(battiti_materni))
        canale_materno=1;
    end
end

function [corretti]=correggi(pic_f,YA,morf)
corretti=pic_f;
stimato=0;
indice_stimato=0;
incrementa=1;
for zm=1:length(pic_f)
stimato=0;
indice_stimato=0;
incrementa=1;
    if(morf==1)
        id=0;
        start=pic_f(zm)-25;
        if((start-8>0)&&(start+51+8<length(YA)))
            while(id<50)
                contat=0;
                for lt=1:8
                    if((YA(start+id)>=YA(start+id-lt))&&(YA(start+id)>=YA(start+id+lt)))
                        contat=contat+1;
                    end
                end
                if(contat==8)
                    stimato(incrementa)=YA(start+id);
                    indice_stimato(incrementa)=start+id;
                    incrementa=incrementa+1;

                end
                id=id+1;
            end
        else
            id=0;
            start=pic_f(zm)-25;
            if((start-8>0)&&(start+51+8<length(YA)))
                while(id<200)
                    contat=0;
                    for lt=1:8
                        if((YA(start+id)<=YA(start+id-lt))&&(YA(start+id)<=YA(start+id+lt)))
                            contat=contat+1;
                        end
                    end
                    if(contat==8)
                        stimato(incrementa)=YA(start+id);
                        indice_stimato(incrementa)=start+id;
                        incrementa=incrementa+1;

                    end
                    id=id+1;
                end
            end
        end
    end
    if(incrementa>1)
        [tr,fr]=max(abs(stimato));
        corretti(zm)=indice_stimato(fr);
    end
end
end





%% confronto con la madre per escludere se troppo simile
function[residui_materni]=escludi(qrs1,qrs2,qrs3,qrs4, mqrs)

    soglia_val=20;
    zi=1;
    l1=length(qrs1);
    l2=length(qrs2);
    l3=length(qrs3);
    l4=length(qrs4);
    contatore_1=0;
    contatore_2=0;
    contatore_3=0;
    contatore_4=0;
    residui_materni=[];
    while (zi<length(qrs1))
        zj=1;
        incremento=1;
        while(zj<length(mqrs))
            
            if(abs(qrs1(zi)- mqrs(zj))<=soglia_val)
                contatore_1= contatore_1 + 1;
                qrs1(zi)=[];
                incremento=0;
                break;
            end
            zj=zj+1;
        end
        zi=zi+incremento;
        
    end
    if(contatore_1>0.5*l1)
        residui_materni=[residui_materni 1];
    end
    
    zi=1;
    while (zi<length(qrs2))
        zj=1;
        incremento=1;
        while(zj<length(mqrs))
            if(abs(qrs2(zi)- mqrs(zj))<=soglia_val)
                contatore_2= contatore_2 + 1;
                qrs2(zi)=[];
                incremento=0;
                break;
            end
            zj=zj+1;
        end
        zi=zi+incremento;
    end
    if(contatore_2>0.5*l2)
        residui_materni=[residui_materni 2];
    end
    
    zi=1;
    while (zi<length(qrs3))
        zj=1;
        incremento=1;
        while(zj<length(mqrs))
            if(abs(qrs3(zi)- mqrs(zj))<=soglia_val)
                contatore_3= contatore_3 + 1;
                qrs3(zi)=[];
                incremento=0;
                break;
            end
            zj=zj+1;
        end
        zi=zi+incremento;
    end

    if(contatore_3>0.5*l3)
        residui_materni=[residui_materni 3];
    end
    
    zi=1;
    while (zi<length(qrs4))
        zj=1;
        incremento=1;
        while(zj<length(mqrs))
            if(abs(qrs4(zi)- mqrs(zj))<=soglia_val)
                contatore_4= contatore_4 + 1;
                qrs4(zi)=[];
                incremento=0;
                break;
            end
            zj=zj+1;
        end
        zi=zi+incremento;
    end
    
    if(contatore_4>0.5*l4)
        residui_materni=[residui_materni 4];
    end
    
    if(length(residui_materni)==4)
        [fg,ni]=min([contatore_1, contatore_2, contatore_3, contatore_4]);
        residui_materni(ni)=[];
    end
        
end

%% periodicity correction con tolleranza molto alta, scelgo il migliore


function[qrs,qrs1,qrs2,qrs3,qrs4]=balda_feto(x,delay,fat)
        [righe,colonne]=size(x);
        qrs=0;
        qrs1=0;
        qrs2=0;
        qrs3=0;
        qrs4=0;
        b0 = [1 0 -1];
        a0 = 1;
        
        b1 = [1 0 -2 0 1];
        a1 = 1;
        
        bs = 1/8*ones(1,8);
        as = 1;
        
        for jn=1:righe
            y0(jn,:)=filter(b0,a0,x(jn,:));
            y0(jn,:)=abs(y0(jn,:));
            y1(jn,:)=filter(b1,a1,x(jn,:));
            y1(jn,:)=abs(y1(jn,:));
            y0(jn,:) = [0,y0(jn,1:end-1)]; % vedere dispense a pag. 60
%             y2(jn,:) = 1.3*y0(jn,:) + 1.1*y1(jn,:);
            y2(jn,:) = 1.1*y0(jn,:) + 1.9*y1(jn,:);
            y3(jn,:)=filter(bs,as,y2(jn,:));
        end
        if(righe>1)
            yt=sum(y2(:,:));
        else
            yt=y2;
        end
        
        ys = filter(bs,as,yt);

        % Introduco un ritardo 4, che si accumula al ritardo 2 precedente
        % vedere dispense a pag. 61
        % ritardo finale = 6.
        
        % Esempio di blanking per evitare rilevazioni multiple
%         M = max(ys);
        nu=floor(length(ys)/250);
        soglia=zeros(1,length(ys));
        soglias=zeros(righe,length(ys));

        for t=1:nu-3
            soglia((t-1)*250 +1 : t*250) = 0.5* max(ys((t-1)*250 +1 : ((t+3)*250))); % la soglia è il 60% del massimo del segnale
            for jk=1:righe
                soglias(jk,(t-1)*250 +1 : t*250) = 0.5* max(y3(jk,(t-1)*250 +1 : ((t+3)*250))); % la soglia è il 60% del massimo del segnale
            end
        end
        soglia((nu-3)*250+1:end)=soglia((nu-3)*250);
        if(righe==4)
            soglias(1,(nu-3)*250+1:end)=soglias(1,(nu-3)*250);
            soglias(2,(nu-3)*250+1:end)=soglias(2,(nu-3)*250);
            soglias(3,(nu-3)*250+1:end)=soglias(3,(nu-3)*250);
            soglias(4,(nu-3)*250+1:end)=soglias(4,(nu-3)*250);
        end
%         tm=tm+delay; %delay a 1000 Hz
%         tm=round(tm./fat);


        
        j = 1;
        qrs = zeros(size(ys)); %vettore dove metterò gli indici di occorrenza di un QRS
        blanking = 0;
        
        for i=1:length(ys)-1
            if (blanking == 0)
                if (ys(i) > soglia(i)) && (ys(i+1) < ys(i)) % esempio di sogliatura.  
                    
                        qrs(j) = i-6;
                        j = j + 1;
                        blanking = 30;
                        
                end
            else
                blanking = blanking -1;
            end
        end

        qrs(j:end) = []; % Elimino la parte inutilizzata del vettore.
        
        if(righe==4)
            j = 1;
            qrs1 = zeros(size(ys)); %vettore dove metterò gli indici di occorrenza di un QRS
            ys=y3(1,:);
            soglia=soglias(1,:);
            blanking = 0;
            for i=1:length(ys)-1
                if (blanking == 0)
                    if (ys(i) > soglia(i)) && (ys(i+1) < ys(i)) % esempio di sogliatura.  

                            qrs1(j) = i-6;
                            j = j + 1;
                            blanking = 30;

                    end
                else
                    blanking = blanking -1;
                end
            end
            qrs1(j:end) = []; % Elimino la parte inutilizzata del vettore.

                    j = 1;
            qrs2 = zeros(size(ys)); %vettore dove metterò gli indici di occorrenza di un QRS
            ys=y3(2,:);
            soglia=soglias(2,:);
            blanking = 0;
            for i=1:length(ys)-1
                if (blanking == 0)
                    if (ys(i) > soglia(i)) && (ys(i+1) < ys(i)) % esempio di sogliatura.  

                            qrs2(j) = i-6;
                            j = j + 1;
                            blanking = 30;

                    end
                else
                    blanking = blanking -1;
                end
            end
            qrs2(j:end) = []; % Elimino la parte inutilizzata del vettore.

                    j = 1;
            qrs3 = zeros(size(ys)); %vettore dove metterò gli indici di occorrenza di un QRS
            ys=y3(3,:);
            soglia=soglias(3,:);
            blanking = 0;
            for i=1:length(ys)-1
                if (blanking == 0)
                    if (ys(i) > soglia(i)) && (ys(i+1) < ys(i)) % esempio di sogliatura.  

                            qrs3(j) = i-6;
                            j = j + 1;
                            blanking = 30;

                    end
                else
                    blanking = blanking -1;
                end
            end
            qrs3(j:end) = []; % Elimino la parte inutilizzata del vettore.

                    j = 1;
            qrs4 = zeros(size(ys)); %vettore dove metterò gli indici di occorrenza di un QRS
            ys=y3(4,:);
            soglia=soglias(4,:);
            blanking = 0;
            for i=1:length(ys)-1
                if (blanking == 0)
                    if (ys(i) > soglia(i)) && (ys(i+1) < ys(i)) % esempio di sogliatura.  

                            qrs4(j) = i-6;
                            j = j + 1;
                            blanking = 30;

                    end
                else
                    blanking = blanking -1;
                end
            end
            qrs4(j:end) = []; % Elimino la parte inutilizzata del vettore.
        end

    end


function [Y, proiezioni, D_ord, V_ord,W]=PCA(X)
[M, N]=size (X);
u=mean(X');
% %Calcolo della media per ogni variabile

 %CENTRATI
h=ones(1,N);
B =X - u'*h;

%matrice di covarianza
C=cov(B');
% C_prova=(B*B')/N;

%Calcolo degli autovalori e degli autovettori
[V,D]=eig(C);

%Riordiniamo gli autovalori e i rispettivi autovettori
[val, ind]=sort(diag(D));
D_ord=zeros(M,M);
V_ord=zeros(M,M);
for m=1:M
    t=ind(M+1-m);
    D_ord(m,m)=D(t,t);
    V_ord(:,m)=V(:,t);
end

%Scelta delle prime L colonne (vogliamo ridurre la dimensionalità a L=1)
L=1;
W(:,L)=V_ord(:,L);

% % %Si possono convertire i dati in Z-scores
% % s=zeros(M,1);
% % s(:,1)=sqrt(diag(C));
% % z=B ./(s*h);

%Proiezione dei dati nella nuova base (prima componene principale)
Y=W'*B;

%Proiezione di tutti i basi sulle M basi
proiezioni=V_ord'*B;
end

function [picchi]=scegli(pic_f, feti, battiti_materni)
    feti=feti+12;
    soglia_val=20;
    qrs_1=pic_f;
    qrs_2=feti;
    contatore_1=0;
    contatore_2=0;
    zi=1;
    while (zi<length(battiti_materni))
        zj=1;
        while(zj<length(qrs_1))
            if(abs(qrs_1(zj)- battiti_materni(zi))<=soglia_val)
                contatore_1= contatore_1 + 1;
                qrs_1(zj)=[];
                break;
            end
            zj=zj+1;
        end

        zj=1;
        while(zj<length(qrs_2))
            if(abs(qrs_2(zj)- battiti_materni(zi))<=soglia_val)
                contatore_2= contatore_2 + 1;
                qrs_2(zj)=[];
                break;
            end
            zj=zj+1;
        end
        zi=zi+1;  
    end
    if(contatore_1<contatore_2)
        picchi=qrs_1;
    else
        picchi=qrs_2;
    end
end

function[template]=scelta_template(battiti_materni, mqrs)
    qrs_1=mqrs(1, mqrs(1,:)>0);
    qrs_1=qrs_1+25*4-1;
    qrs_2=mqrs(2, mqrs(2,:)>0);
    qrs_2=qrs_2+25*4-1;
    soglia_val=20*4;
    zi=1;
    contatore_1=0;
    contatore_2=0;
    while (zi<length(battiti_materni))
        zj=1;
        while(zj<length(qrs_1))
            if(abs(qrs_1(zj)- battiti_materni(zi))<=soglia_val)
                contatore_1= contatore_1 + 1;
                qrs_1(zj)=[];
                break;
            end
            zj=zj+1;
        end

        zj=1;
        while(zj<length(qrs_2))
            if(abs(qrs_2(zj)- battiti_materni(zi))<=soglia_val)
                contatore_2= contatore_2 + 1;
                qrs_2(zj)=[];
                break;
            end
            zj=zj+1;
        end
        zi=zi+1;  
    end
    if(contatore_1>=contatore_2)
        template=1;
    elseif(contatore_2>contatore_1)
        template=2;
    end
end

function[qrs]=balda(x,tipo)
%     if(tipo==0)
         b0 = [1 0 -1];
        a0 = 1;
        
        b1 = [1 0 -2 0 1];
        a1 = 1;
      
        y0 = filter(b0,a0,x);
        y0 = abs(y0);

        y1 = filter(b1,a1,x);
        y1 = abs(y1);

        y0 = [0,y0(1:end-1)]; % vedere dispense a pag. 60
        y2 = 1.3*y0 + 1.1*y1;
        % accumulo un ritardo 2

        bs = 1/8*ones(1,8);
        as = 1;
        ys = filter(bs,as,y2);
        % Introduco un ritardo 4, che si accumula al ritardo 2 precedente
        % vedere dispense a pag. 61
        % ritardo finale = 6.

        % Esempio di blanking per evitare rilevazioni multiple
%         M = max(ys);
        nu=floor(length(ys)/250);
        soglia=zeros(1,length(ys));
        for t=1:nu-3
            soglia((t-1)*250 +1 : t*250) = 0.6* max(ys((t-1)*250 +1 : ((t+3)*250))); % la soglia è il 60% del massimo del segnale
        end
        soglia((nu-3)*250+1:end)=soglia((nu-3)*250);
        j = 1;
        qrs = zeros(size(ys)); %vettore dove metterò gli indici di occorrenza di un QRS
        blanking = 0;
        if(tipo==0)
            for i=1:length(ys)-1
                if (blanking == 0)
                    if (ys(i) > soglia(i)) && (ys(i+1) < ys(i)) % esempio di sogliatura.
                        sta=i;
                        avan=1;
                        tempo=0;
                        while(ys(sta)>soglia(sta))
                            tempo(avan)=ys(sta);
                            avan=avan+1;
                            sta=sta-1;
                            if(sta==0)
                                break;
                            end
                        end
                        [mj,ij]=max(tempo);
                        ij=ij-1;

                        qrs(j) = i-6-ij;
                        j = j + 1;
                        blanking = 30;

                    end
                else
                    blanking = blanking -1;
                end
            end

            qrs(j:end) = []; % Elimino la parte inutilizzata del vettore.
        else
            blanking = 0;

            for i=1:length(ys)-1
                if (blanking == 0)
                    if (ys(i) > soglia(i)) && (ys(i+1) < ys(i)) % esempio di sogliatura.

                        qrs(j) = i-6;
                        j = j + 1;
                        blanking = 30;

                    end
                else
                    blanking = blanking -1;
                end
            end

            qrs(j:end) = [];
            
        end
%         figure
%         plot(x)
%         hold on
%         plot(qrs,x(qrs),'*r')
%     else
% 
%         b0 = 1/2*[1 0 -1];  %ritardo di 1
%         a0 = 1;
%         delay0=1;
%         
%         b1 = 1/4*[1 0 -2 0 1];  %ritardo di 2
%         a1 = 1;
%         delay1=2;
% 
%         b2 = 1/44*[2 4 8 16 8 4 2];
%         delay2=3;
%         a2=1;
%         %ritardo di 3        
%       
%         y0 = filter(b0,a0,x);
%         y0 = abs(y0);
% 
%         y1 = filter(b1,a1,x);
%         y1 = abs(y1);
%         
%         y2=filter(b2,a2,y0); 
%         y2=filter(b2,a2,x);
%         
%         for f=2:length(y2)
%             y3(f)=y2(f)*y2(f-1); %ritardo di 1
%         end
%         delay3=1;
% %         y3=y2;
% %         delay3=0;
%         
%         
%         %ritardi
%         dl=delay0+delay2+delay3-delay1;
%         y1=[zeros(1,dl),y1(1:end-dl)];
%         y4=y3+0.5*y1;
%         
% 
% 
%         bs = 1/8*ones(1,8);
%         as = 1;
%         ys = filter(bs,as,y4);
%         delays=4;
%          total_delay=delay0+delay1+delay2+delay3+delays;
%         
%         figure
%         del=total_delay-delay0;
%         plot([zeros(1,del) y0(1:300-del)])
%         hold on
%         plot([zeros(1,total_delay) x(1:300-total_delay)],'m')
%         del2=total_delay-delay1;
%         plot([zeros(1,del2) y1(1:300-del2)],'r')
%         del=del-delay2;
%         plot([zeros(1,del) y2(1:300-del)],'c')
%         del=del-delay3;
%         plot([zeros(1,del) y3(1:300-del)],'g')
%         plot([zeros(1,del) y4(1:300-del)],'k')
%         plot(ys(1:300),'y')
%         title('blu der prima, rosso der sec, ciano shaped filter, verde amplificato,nero comb lin, giallo finale, magenta ingresso')
% %         
%         % Introduco un ritardo 4, che si accumula al ritardo 5 precedente
%         % ritardo finale = 9.
% 
%         % Esempio di blanking per evitare rilevazioni multiple
% %         M = max(ys);
%         nu=floor(length(ys)/250);
%         soglia=zeros(1,length(ys));
%         for t=1:nu-3
%             soglia((t-1)*250 +1 : t*250) = 0.5* max(ys((t-1)*250 +1 : ((t+3)*250))); % la soglia è il 60% del massimo del segnale
%         end
%         soglia((nu-3)*250+1:end)=soglia((nu-3)*250);
%         j = 1;
%         qrs = zeros(size(ys)); %vettore dove metterò gli indici di occorrenza di un QRS
%         blanking = 0;
% 
%         for i=1:length(ys)-1
%             if (blanking == 0)
%                 if (ys(i) > soglia(i)) && (ys(i+1) < ys(i))   
%                     
%                         qrs(j) = i;
%                         j = j + 1;
%                         blanking = 30;
%                     
%                 end
%             else
%                 blanking = blanking -1;
%             end
%         end
% 
%         qrs(j:end) = []; % Elimino la parte inutilizzata del vettore.
%         figure
%         plot(x)
%         hold on
%         plot(qrs,x(qrs),'*r')
%     end

end

function [FilteredECG, original_d] = preprocessing(ECG,F_acq)
% ---- preprocess the data ----
%controllo che non ci siano NaN
for mk=1:4
    seg=ECG(mk,:);
    l=isnan(seg);  %mette 1 logico dove ci sono NaN
    z=find(l==1);  %trovo gli 1
    z=[z,0];
    c=1;
    while c<length(z)
        if(z(c)>0)
            start=z(c);
            stop=z(c);
            if(z(c)==length(seg))
        
            else
                while(z(c+1)-z(c)==1)
                    c=c+1;
                    stop=z(c);
                end
            end
        end
        c=c+1;
        delta=stop-start;
        if(delta==0) %ho solo un campione
            delta=3;
        end
        if(((start-delta)>0)&&((stop+delta)<length(seg)+1))
            %interp sline
            a=[start-delta:start-1, stop+1:stop+delta];
        elseif((start-delta)<1)
            %spline da inizio segnale
            if(start>1)
                a=[1:start-1, stop+1:stop+delta];
            else
                a=[start,stop+1: stop + delta];
            end
        elseif((stop+delta)>length(seg))
            %spline fino alla fine del segnale
            if(stop<length(z))
                a=[start-delta:start-1, stop+1:length(seg)];
            else
                a=[start-delta:start-1, stop];
            end
        end
        b=seg(a);
        d=(start:stop);
        tt=spline(a,b,d);
        seg(start:stop)=tt;
    end
    ECG(mk,:)=seg;
    
end




F_elab=250;
ECG=ECG-mean(ECG,2)*ones(1,length(ECG));
for z=1:4
    data(z,:)=DERIVATORE(ECG(z,:),F_acq,1);
end
s=zeros(1,4);
for z=1:4
    s(z)=std(data(z,:));
end
a=10;
soglie=a*s;
for l=1:4
    ECG(l,:)=taglia_outlier2(ECG(l,:),data(l,:),soglie(l));
end
% % load coeff_FIR_46_50_1000.mat
% Hd = prova_filtro;
% Num=Hd.Numerator;
% for l=1:4
%     ECG_lp(l,:)=filter(Num,1,ECG(l,:));
% end
% delay1=floor(length(Num)/2);

ingresso=ECG;

Hd = prova_iir_lp;
Num=Hd.Numerator;
Den=Hd.Denominator;
% delay1=length(Num)-1;
for l=1:4
    ECG(l,:)=filtfilt(Num,Den,ECG(l,:));
end


 %@1000Hz
% dawn(:,1)=resample(ECG(1,:)',F_elab,F_acq);
% dawn(:,2)=resample(ECG(2,:)',F_elab,F_acq);
% dawn(:,3)=resample(ECG(3,:)',F_elab,F_acq);
% dawn(:,4)=resample(ECG(4,:)',F_elab,F_acq);
% ECG=dawn';
% % load coeff_FIR_1_2.mat
% Hd_2 = prova_filtro2;
% Num=Hd_2.Numerator;
% for l=1:4
%     ECG_fir(l,:)=filter(Num,1,ECG(l,:));
% end
% delay2=floor(length(Num)/2); %@250Hz

% Hd=prova_iir_hp;
% Hd = filtro_iir_4_8_250;
Hd = hp_4_8_1000;
Num=Hd.Numerator;
Den=Hd.Denominator;
% delay2=length(Num)-1;
for l=1:4
    ECG(l,:)=filtfilt(Num,Den,ECG(l,:));
end


% delay2=4*delay2;
 original_d= ECG;
% delay=delay1+delay2;
delay=0;

prova=original_d(1,:);
% prova=resample(prova',F_acq,F_elab);
% figure
% plot(prova')
% hold on
% plot(ingresso(1,1:end-delay),'r')
% title('filtraggio')

dawn(:,1)=resample(ECG(1,:)',F_elab,F_acq);
dawn(:,2)=resample(ECG(2,:)',F_elab,F_acq);
dawn(:,3)=resample(ECG(3,:)',F_elab,F_acq);
dawn(:,4)=resample(ECG(4,:)',F_elab,F_acq);
FilteredECG=dawn';

end

function Hd = hp_4_8_1000   %filtro iir passa alto a 1000Hz
%HP_4_8_1000 Returns a discrete-time filter object.

%
% M-File generated by MATLAB(R) 7.5 and the Signal Processing Toolbox 6.8.
%
% Generated on: 28-Aug-2013 11:15:05
%

% Butterworth Highpass filter designed using FDESIGN.HIGHPASS.

% All frequency values are in Hz.
Fs = 1000;  % Sampling Frequency

Fstop = 4;           % Stopband Frequency
Fpass = 8;           % Passband Frequency
Astop = 40;          % Stopband Attenuation (dB)
Apass = 1;           % Passband Ripple (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);
% Get the transfer function values.
[b, a] = tf(Hd);

% Convert to a singleton filter.
Hd = dfilt.df2(b, a);
end

function Hd = filtro_iir_4_8_250
%FILTRO_IIR_4_8_250 Returns a discrete-time filter object.

%
% M-File generated by MATLAB(R) 7.5 and the Signal Processing Toolbox 6.8.
%
% Generated on: 23-Aug-2013 19:19:44
%

% Butterworth Highpass filter designed using FDESIGN.HIGHPASS.

% All frequency values are in Hz.
Fs = 250;  % Sampling Frequency

Fstop = 4;           % Stopband Frequency
Fpass = 8;           % Passband Frequency
Astop = 40;          % Stopband Attenuation (dB)
Apass = 1;           % Passband Ripple (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);
% Get the transfer function values.
[b, a] = tf(Hd);

% Convert to a singleton filter.
Hd = dfilt.df2(b, a);

% [EOF]
end

function Hd = prova_iir_lp
%PROVA_IIR_LP Returns a discrete-time filter object.

%
% M-File generated by MATLAB(R) 7.5 and the Signal Processing Toolbox 6.8.
%
% Generated on: 22-Aug-2013 10:13:58
%

% Chebyshev Type II Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are in Hz.
Fs = 1000;  % Sampling Frequency

Fpass = 46;          % Passband Frequency
Fstop = 50;          % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 40;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd = design(h, 'cheby2', 'MatchExactly', match);
% Get the transfer function values.
[b, a] = tf(Hd);

% Convert to a singleton filter.
Hd = dfilt.df2(b, a);

% [EOF]
end


function Hd = prova_iir_hp
%PROVA_IIR_HP Returns a discrete-time filter object.

%
% M-File generated by MATLAB(R) 7.5 and the Signal Processing Toolbox 6.8.
%
% Generated on: 21-Aug-2013 16:09:04
%

% Butterworth Highpass filter designed using FDESIGN.HIGHPASS.

% All frequency values are in Hz.
Fs = 250;  % Sampling Frequency

Fstop = 1;           % Stopband Frequency
Fpass = 2;           % Passband Frequency
Astop = 40;          % Stopband Attenuation (dB)
Apass = 1;           % Passband Ripple (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);
% Get the transfer function values.
[b, a] = tf(Hd);

% Convert to a singleton filter.
Hd = dfilt.df2(b, a);

% [EOF]
% % % % Fs = 250;  % Sampling Frequency
% % % % 
% % % % Fstop = 1;           % Stopband Frequency
% % % % Fpass = 3;           % Passband Frequency
% % % % Astop = 40;          % Stopband Attenuation (dB)
% % % % Apass = 1;           % Passband Ripple (dB)
% % % % match = 'stopband';  % Band to match exactly
% % % % 
% % % % % Construct an FDESIGN object and call its BUTTER method.
% % % % h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
% % % % Hd = design(h, 'butter', 'MatchExactly', match);
% % % % % Get the transfer function values.
% % % % [b, a] = tf(Hd);
% % % % 
% % % % % Convert to a singleton filter.
% % % % Hd = dfilt.df2(b, a);
% % % % 
% % % % % [EOF]
end

function Hd = prova_filtro
%PROVA_FILTRO Returns a discrete-time filter object.

%
% M-File generated by MATLAB(R) 7.5 and the Signal Processing Toolbox 6.8.
%
% Generated on: 19-Aug-2013 16:46:53
%

% FIR Window Lowpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 1000;  % Sampling Frequency

Fpass = 46;              % Passband Frequency
Fstop = 50;              % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.001;           % Stopband Attenuation
flag  = 'scale';         % Sampling Flag

% Calculate the order from the parameters using KAISERORD.
[N,Wn,BETA,TYPE] = kaiserord([Fpass Fstop]/(Fs/2), [1 0], [Dstop Dpass]);


% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
Hd = dfilt.dffir(b);
end

function Hd = prova_filtro2
%PROVA_FILTRO2 Returns a discrete-time filter object.

%
% M-File generated by MATLAB(R) 7.5 and the Signal Processing Toolbox 6.8.
%
% Generated on: 19-Aug-2013 17:36:08
%

% FIR Window Highpass filter designed using the FIR1 function.

% All frequency values are in Hz.
Fs = 250;  % Sampling Frequency

Fstop = 1;               % Stopband Frequency
Fpass = 2;               % Passband Frequency
Dstop = 0.01;            % Stopband Attenuation
Dpass = 0.057501127785;  % Passband Ripple
flag  = 'scale';         % Sampling Flag

% Calculate the order from the parameters using KAISERORD.
[N,Wn,BETA,TYPE] = kaiserord([Fstop Fpass]/(Fs/2), [0 1], [Dpass Dstop]);


% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
Hd = dfilt.dffir(b);

% [EOF]
end

function [uscita]=taglia_outlier2(original,ingresso,soglia)

t=find(abs(ingresso)>soglia);
uscita=original;
if(isempty(t))
else
    j=t;
    for v=(length(t)):-1:2
        if(abs(t(v)-t(v-1))==1)
            j(v)=0;
        end
    end
    for z=1:length(j)
        indice=j(z);
        if(indice==0)
        else
            if((indice-5)<1)
                delta=original(1:(indice+10) );
                th=abs(median(delta));
                fb=indice+10;
                while(abs(original(fb)-original(fb-1))<3*th)
                    fb=fb-1; %estremo destro
                    if(fb==1)
                        fb=indice+10;
                        break;
                    end
                end
                estremo_dx=fb;
                uscita(1,1:estremo_dx)=original(estremo_dx);
            elseif((indice+10)>length(original)-1)
                delta=original((indice-5):end );
                th=abs(median(delta));
                fb=indice-5;
                while(abs(original(fb)-original(fb+1))<3*th)
                    fb=fb+1; %estremo sinistro
                    if(fb==length(original))
                        fb=indice-5;
                        break;
                    end
                end
                estremo_sx=fb;
                uscita(1,estremo_sx:end)=original(estremo_sx);
            else
                delta=original((indice-5):(indice+10) );
                th=abs(median(delta));
                fb=indice-5;
                while(abs(original(fb)-original(fb+1))<3*th)
                    fb=fb+1; %estremo sinistro
                    if(fb>=indice+10)
                        fb=indice+10;
                        break;
                    end
                end
                estremo_sx=fb;

                fb=indice+10;
                while(abs(original(fb)-original(fb-1))<3*th)
                    fb=fb-1; %estremo destro
                    if(fb<=indice-5)
                        fb=indice-5;
                        break;
                    end
                end
                estremo_dx=fb;
                c=(estremo_dx - estremo_sx)-1;
                if(c>1)
                    a=[(estremo_sx-c-1):estremo_sx, estremo_dx: (estremo_dx+c+1)];
                    b=original(a);
                    d=(estremo_sx+1 : estremo_dx - 1);
                    tt=spline(a,b,d);
                    uscita(estremo_sx+1 : estremo_dx - 1)=tt;
                elseif c==1
                    a=[(estremo_sx-1):estremo_sx, estremo_dx: (estremo_dx+1)];
                    b=original(a);
                    d=estremo_sx+1;
                    tt=spline(a,b,d);
                    uscita(estremo_sx+1 : estremo_dx - 1)=tt;
                end
            end
        end
    end
end
end

function y=DERIVATORE(x,fs,parametro)
    switch(parametro)
        case 0
            b=(1/2)*[1 -1];
            a=1;
            y=filter(b,a,x);
        case 1
            b=(1.99/2)*[1 -1];
            a=[1 -0.99];
            y=filter(b,a,x);
    end
end


function [pulito]=MECGcancellation(segnale,mqrs,Apost ,template)
    if(template==1)
        temp_sottr(1,:)=Apost(1,:);
        temp_sottr(2,:)=Apost(2,:);
        temp_sottr(3,:)=Apost(3,:);
    else
        temp_sottr(1,:)=Apost(4,:);
        temp_sottr(2,:)=Apost(5,:);
        temp_sottr(3,:)=Apost(6,:);
    end        
    dim=length(temp_sottr);
    for tm=1:3
        temp_sottr(tm,:)=temp_sottr(tm,:).*(tukeywin(dim, 0.2))';
    end
%     temp_sottr=Apost(template,:);
    pulito=segnale;
    delta_stop=0;
%     indice_debug=1;
    for lt=1:length(mqrs)
        inds=mqrs(lt);
        if((inds+dim)>length(segnale))
            ui=length(segnale);
            delta_stop=inds+dim-ui;
%             segnale=[segnale,zeros(1,delta)];
            tempn=[segnale(inds-dim:end), zeros(1,delta_stop)];
            inizio=dim-1;
%             tempn=segnale(length(segnale)-101:length(segnale));
        elseif((inds-dim)<1)
%             tempn=segnale(1:dim+1);
            tempn=[zeros(1,-inds+dim+1) segnale(1:inds+dim)];
%             inizio=dim-1+(inds-dim)-1;
            inizio=dim-1;
        else
            tempn=segnale(inds-dim:inds+dim);
            inizio=dim-1;
        end
        r1 = pearsoncoeff(tempn,temp_sottr(1,:));
        r2 = pearsoncoeff(tempn,temp_sottr(2,:));
        r3 = pearsoncoeff(tempn,temp_sottr(3,:));
        [mar(1),inr(1)]=max(r1);
        [mar(2),inr(2)]=max(r2);
        [mar(3),inr(3)]=max(r3);
        [ma,in]=max(mar);
%         start=inds-inizio-2+in;
        start=inds-inizio-2+inr(in);
        delta_start=0;
        if(start<1)
             start=1;
            delta_start=-inds+dim;
        end
        if(delta_stop>1)
            delta_stop=start+dim-length(segnale)-1;
        end
        if(delta_stop<0)
            delta_stop=0;
        end
        if(delta_stop<0)
            delta_stop=0;
        end
%         if(start>length(segnale))
%             
%         end
%         G=max(abs(segnale(start+delta_start:start+dim-1-delta_stop)))/max(abs(temp_sottr(in,delta_start+1:end-delta_stop)));
%         deb2=(temp_sottr*G);
%         deb1=segnale(start:start+dim-1);
        deb1=segnale(start:start+dim-1-delta_stop-delta_start);
%         delta_pos=abs(max((segnale(start:start+49))))-abs(max((temp_sottr)));
%         delta_neg=abs(min((segnale(start:start+49))))-abs(min((temp_sottr)));
        Gpos=max((segnale(start:start+dim-1-delta_stop-delta_start)))/max((temp_sottr(in,delta_start+1:end-delta_stop)));
        Gneg=min((segnale(start:start+dim-1-delta_stop-delta_start)))/min((temp_sottr(in,delta_start+1:end-delta_stop)));
        deb2=0;
        for kl=1:length(deb1)
            if(temp_sottr(in,kl)>0)
                deb2(kl)=temp_sottr(in,delta_start+kl)*Gpos;
            else
                deb2(kl)=temp_sottr(in,delta_start+kl)*Gneg;
            end
        end
%         if(delta_pos>delta_neg)
%                 deb2=(temp_sottr*Gpos);
%         else
%                 deb2=(temp_sottr*Gneg);
%         end
%         grt=0;
%         if(delta>0)
%            grt=start+dim-1-length(pulito); 
%         end
%         if(grt<0)
%             grt=0;
%         end
%         if(grt==0)
%             matrice_debug(indice_debug,:)=segnale(start:start+dim-1-grt);
%             
%             indice_debug=indice_debug+1;
%         end
%         deb1=segnale(start+delta_start:start+dim-1-delta_stop);
        pulito(start:start+dim-1-delta_start-delta_stop)=segnale(start:start+dim-1-delta_start-delta_stop)-(deb2(1:end));
%         figure
%         plot(segnale(start:start+dim-1-grt))
%         hold on
%         plot((deb2(1:end-grt)),'r')
%         title('in rosso template')
    end
end



function [qrs_f, DEL] =PeakDetectionFetali(FECG,Rmaterni)
    AM=8;
    DEL=floor((AM)/2)+1;
    
%     a=derivatore_elementare(FECG);
%     b=a.^2;
%     d=b;
    
    b0 = [1 0 -1];
        a0 = 1;
        
        b1 = [1 0 -2 0 1];
        a1 = 1;
      
        y0 = filter(b0,a0,FECG);
        y0 = abs(y0);

        y1 = filter(b1,a1,FECG);
        y1 = abs(y1);

        y0 = [0,y0(1:end-1)]; % vedere dispense a pag. 60
        y2 = 1.3*y0 + 1.1*y1;
        % accumulo un ritardo 2

        bs = 1/8*ones(1,8);
        as = 1;
        ys = filter(bs,as,y2);
        DEL=DEL+1;
%     for lt=1:4
%         a(lt,:)=derivatore_elementare(FECG(lt,:));
%         b(lt,:)=a(lt,:).^2;
%     end
%     d=sum(b);
%     c=media_mobile(d,AM);
    [so]=sogliatrovapiccoR_offline(ys,0.8);
%     ys=c;
    j=1;
    t=0;
    qrs_f=[];
    REFRATT=15;  %bpm max 300
    for i=1:length(ys)-1
        if(i>t)
            if(ys(i)>so(i))
                r=ys(i);
                ri=i;
                k=i;
                while((ys(k)>so(k))&&(k<(length(ys))))
                    if(ys(k)>r)
                        r=ys(k);
                        ri=k;
                    end
                    k=k+1;
                end
                qrs_f(j)=ri;
                if(j>1)
                    if(abs(qrs_f(j)-qrs_f(j-1))<=REFRATT)
%                         if(ys(qrs_f(j))>ys(qrs_f(j-1)))
                            pea=qrs_f(j)*ones(1,length(Rmaterni));
                            p=min(abs(Rmaterni - pea));
                            qea=qrs_f(j-1)*ones(1,length(Rmaterni));
                            q=min(abs(Rmaterni - qea));
                            if((p>q)&&(q<REFRATT))
                                qrs_f(j-1)=qrs_f(j);
                                qrs_f(j)=[];
                            elseif((q>p)&&(p<REFRATT))
                                qrs_f(j)=[];
                            elseif(ys(qrs_f(j))>ys(qrs_f(j-1)))  
                                qrs_f(j-1)=qrs_f(j);
                                qrs_f(j)=[];
                            elseif(ys(qrs_f(j))<ys(qrs_f(j-1)))
                                qrs_f(j)=[];
                            end
                    else
                        j=j+1;
                    end
                else
                    j=j+1;
                end
                t=k;
            end
        end
    end
end





function [Ttotpost,Ttot,occ, Apost,Nt,occorrenze_finali]=PeakDetection(signal, wavname, levels, clearapprox, buffer_len, ...
    threshold, hard_soft, deepth_th, level_independent, single_level_th_dep, level_th, scaling_type, ...
    rms_sc, N, th_base_type, LEN, par_detection,Thcross,Thc,Thfusione, original_d)
mid = LEN;
occorrenze_finali=[];
occ=0;
mid_mezzi = LEN/2;
% Thc = 0.9; % soglia di correlazione (0.9)
% Thfusione=0.9; % soglia di correlazione per fusione template
Kin = 0; % indice iniziale QRS
Kfn = 0; % indice finale QRS
Ttotdisc = 0; % conteggio QRS trovati ma non classificati
NumTempl=2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thcross=0.85;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIZIO CREAZIONE VARIABILI
NTM = 10; % numero massimo di template ammessi
Nt = 0; % numero di template presenti in quel momento
A = zeros(NTM,LEN); % matrice dei template
Templ_overwritten = 0;
Ttot = zeros(NTM,1); % vettore delle occorrenze totali
tmp = zeros(1,LEN); % buffer temporaneo per uno spike
cc = zeros(NTM,LEN);% matrice contenente le crosscorrelazioni fra template e spike trovato
Ld = cell(levels,1);
Hd = cell(levels,1);
Lr = cell(levels,1);
Hr = cell(levels,1);
% disp('creazione filtri wavelet...')
[Ld{1},Hd{1},Lr{1},Hr{1}]=wfilters(wavname); % predisposizione wavelet 'bior6.8'
for i=2:levels % creazione filtri a trous di Mallat
    Ld{i}= upsampling(Ld{i-1});
    Hd{i}= upsampling(Hd{i-1});
    Lr{i}= upsampling(Lr{i-1});
    Hr{i}= upsampling(Hr{i-1});
end
for i=2:levels % eliminazione code di zeri filtri
    t=(2^(i-1))-1;
    Ld{i}(end-t+1:end)=[];
    Hd{i}(end-t+1:end)=[];
    Lr{i}(end-t+1:end)=[];
    Hr{i}(end-t+1:end)=[];
end
% disp('definizione ritardi...')
delay=cell(levels-1,1);
delayIT = 0;
for i = levels:-1:2
    if rem(length(Ld{i}),2)==0 % filtro lunghezza pari (Ld e Hd sono della stessa lunghezza)
        delayIT=delayIT+length(Ld{i});
    else                      % filtro lunghezza dispari (Ld e Hd sono della stessa lunghezza)
        delayIT=delayIT+length(Ld{i})-1;
    end
    delay{i-1}=delayIT;
end
if rem(length(Hd{1}),2)==0 % filtro lunghezza pari (Ld e Hd sono della stessa lunghezza)
    total_delay=delayIT+length(Hd{1});
else                      % filtro lunghezza dispari (Ld e Hd sono della stessa lunghezza)
    total_delay=delayIT+length(Hd{1})-1;
end
% disp('creazione buffer dei filtri in decomposizione...')
input_buffer = zeros(buffer_len + length(Ld{1}) - 1,1);
dec_buffer_L = cell(levels,1);
dec_buffer_H = cell(levels,1);
for i = 1:levels - 1
    dec_buffer_L{i} = zeros(buffer_len + length(Ld{i+1}) - 1,1);
    dec_buffer_H{i} = zeros(buffer_len + delay{i} + length(Hr{i}) - 1,1);
end
dec_buffer_L{levels} = zeros(buffer_len + length(Ld{levels}) - 1,1);
dec_buffer_H{levels} = zeros(buffer_len + length(Hd{levels}) - 1,1);
% disp('creazione buffer dei filtri in ricomposizione...')
tmpbuffer1rec = zeros(buffer_len,1);
tmpbuffer2rec = zeros(buffer_len,1);
rec_buffer_sum = cell(levels,1);
for i = 1:levels
    rec_buffer_sum{i} = zeros(buffer_len + length(Lr{i}) - 1,1);
end
denoised = zeros(size(signal));
% disp('creazione variabili per soglia...')
buffer_for_th = zeros(deepth_th * buffer_len, levels);
th = zeros(levels, 1);
update_th = 1;
scaling = zeros(levels, 1);
switch th_base_type
    case 100
        th_base = 0.3936 + 0.1829*(log2(N));
    case 200
        th_base = sqrt(2*log(N));
    otherwise
        th_base = th_base_type;
end
if scaling_type == 5
    th_base = 1;
end
% fprintf('Processing...');
numblocks = length(signal)/buffer_len;
% s = strcat('Processing (', num2str(length(signal)/buffer_len));
% s = strcat(s, ' blocks)...');
% hwait = waitbar(0,'Processing block -0- of -20-','name',s,'visible','off');
% hwaitpos = get(hwait,'position');  % [left bottom width height]
% hwaitpos(1,1) = hwaitpos(1,1) - hwaitpos(1,3)/2 - 2;
% set(hwait, 'position', hwaitpos);
% s = strcat('Templates created (', num2str(NTM));
% s = strcat(s, ' max)...');
% hwait2 = waitbar(0,'Number of templates at now: -0-','name',s,'visible','off');
% hwaitpos = get(hwait2,'position');  % [left bottom width height]
% hwaitpos(1,1) = hwaitpos(1,1) + hwaitpos(1,3)/2 + 2;
% set(hwait2, 'position', hwaitpos);
% set(hwait2, 'visible', 'on');
% set(hwait, 'visible', 'on');
for i=1:numblocks
    input_buffer(1:end - buffer_len)=input_buffer(buffer_len + 1:end);
%     s = strcat('Processing block -', num2str(i));
%     s = strcat(s,'- of -20-');
%     waitbar(i/numblocks, hwait, s);
    input_buffer(length(Ld{1}):end)=signal((i-1)*buffer_len+1:(i)*buffer_len);
    tmpl = filter(Ld{1},1,input_buffer);
    tmph = filter(Hd{1},1,input_buffer);
    for lev = 1:levels - 1
        dec_buffer_L{lev}(1:end - buffer_len) = dec_buffer_L{lev}(buffer_len + 1:end);
        dec_buffer_H{lev}(1:end - buffer_len) = dec_buffer_H{lev}(buffer_len + 1:end);
        dec_buffer_L{lev}(end - buffer_len + 1:end) = tmpl(length(Ld{lev}):end);
        dec_buffer_H{lev}(end - buffer_len + 1:end) = tmph(length(Hd{lev}):end);
        tmpl = filter(Ld{lev+1},1,dec_buffer_L{lev});
        tmph = filter(Hd{lev+1},1,dec_buffer_L{lev});
    end
    dec_buffer_L{levels}(1:end - buffer_len) = dec_buffer_L{levels}(buffer_len + 1:end);
    dec_buffer_H{levels}(1:end - buffer_len) = dec_buffer_H{levels}(buffer_len + 1:end);
    dec_buffer_L{levels}(end - buffer_len + 1:end) = tmpl(length(Lr{levels}):end);
    dec_buffer_H{levels}(end - buffer_len + 1:end) = tmph(length(Hr{levels}):end);
    if i == 1
        for lev = 1:levels
            for j = 1:deepth_th
                buffer_for_th((j-1)*buffer_len+1:j*buffer_len, lev) = dec_buffer_H{lev}(end-buffer_len+1:end);
            end
        end
        if threshold ~= 0 % se è previsto di sogliare
            if level_independent ~= 0 % se la soglia non dipende dal segnale alle varie scale
                for lev = 1:levels
                    scaling(lev) = 1; % nessuno scalamento alla soglia th di base
                end
            else
                if single_level_th_dep ~= 0
                    scaling(1) = scal_th(scaling_type, level_th, buffer_for_th, rms_sc);
                    for lev = 2:levels
                        scaling(lev) = scaling(lev-1); % propago lo stesso scalamento su tutte le scale
                    end
                else
                    for lev = 1:levels
                        scaling(lev) = scal_th(scaling_type, lev, buffer_for_th, rms_sc);
                    end
                end
            end
            for lev = 1:levels
                th(lev) = th_base * scaling(lev);
            end
        end
    end
    if clearapprox ~= 0
        dec_buffer_L{levels} = zeros(buffer_len + length(Ld{levels}) - 1,1);
    end
    % applica soglia
    if (hard_soft == 's') || (hard_soft == 'S' ) % se specificato soft (s oppure S) fa soft
        for lev=1:levels
            dec_th = (abs(dec_buffer_H{lev}(end-buffer_len+1:end)) - th(lev));
            dec_th = (dec_th+abs(dec_th))/2;
            dec_buffer_H{lev}(end-buffer_len+1:end) = ...
                sign(dec_buffer_H{lev}(end-buffer_len+1:end)).*dec_th;
        end
    else % altrimenti hard
        for lev=1:levels
            dec_buffer_H{lev}(end-buffer_len+1:end) = ...
                dec_buffer_H{lev}(end-buffer_len+1:end).*(abs(dec_buffer_H{lev}(end-buffer_len+1:end)) > th(lev));
        end
    end
    % ricomposizione
    rec_buffer_sum{levels} = dec_buffer_L{levels};
    for lev = levels-1:-1:1
        % filtra Lr su tmpl e lo copia su tmp1 (256)
        tmpl = filter(Lr{lev+1},1,rec_buffer_sum{lev+1});
        tmpbuffer1rec = tmpl(end-buffer_len+1:end);
        % filtra Hr su tmpr e lo copia su tmp2 (256)
        tmph = filter(Hr{lev+1},1,dec_buffer_H{lev+1}(1:length(Hr{lev+1})-1+buffer_len));
        tmpbuffer2rec = tmph(end-buffer_len+1:end);
        % aggiorna il rec_buffer_sum del livello sottostante
        rec_buffer_sum{lev}(1:end - buffer_len) = rec_buffer_sum{lev}(buffer_len + 1:end);
        % somma tmp1 e tmp2 in rec_buffer_sum (posizione opportuna)
        rec_buffer_sum{lev}(end - buffer_len + 1:end)=0.5*(tmpbuffer1rec + tmpbuffer2rec);
    end
    % filtra Lr su tmpl e lo copia su tmp1 (256)
    tmpl = filter(Lr{1},1,rec_buffer_sum{1});
    tmpbuffer1rec = tmpl(end-buffer_len+1:end);
    % filtra Hr su tmpr e lo copia su tmp2 (256)
    tmph = filter(Hr{1},1,dec_buffer_H{1}(1:length(Hr{1})-1+buffer_len));
    tmpbuffer2rec = tmph(end-buffer_len+1:end);
    % somma tmp1 e tmp2 in rec_buffer_sum (posizione opportuna)
    denoised((i-1)*buffer_len+1:(i)*buffer_len) = 0.5*(tmpbuffer1rec + tmpbuffer2rec);
end
%  QRS DETECTOR
if(par_detection~=1)
    [ys, so]=qrs_detector_wd(denoised);
    delay_stadi_QRSdet = 7;  %2 della der seconda + 5 media mobile
else
    a=derivatore_elementare(denoised);
    b=a.^2;
    ys=media_mobile(b,8);
    [so]=sogliatrovapiccoR_offline(ys,0.9);  %AUMENTATA SOGLIA A 0.9
    delay_stadi_QRSdet = 5;  %1 della der prima + 4 media mobile
end
j=1;
t=0;
QRSs=[];
for i=1:length(ys)-1
    if(i>t)
        if(ys(i)>so(i))
            r=ys(i);
            ri=i;
            k=i;
            while((ys(k)>so(k))&&(k<(length(ys))))
                if(ys(k)>r)
                    r=ys(k);
                    ri=k;
                end
                k=k+1;
            end
            QRSs(j)=ri;
            j=j+1;
            t=k;
        end
    end
end
% figure
% plot(ys)
% hold on
% plot(so,'g')
% if tmatch_on_den~=0
%     %% Template matching sul segnale denoised
%     BUFF = denoised;
%     QRSs=QRSs-delay_stadi_QRSdet;
% else
%     %% Template matching sul segnale originale
%     BUFF = signal;
%     QRSs=QRSs-delay_stadi_QRSdet-total_delay; %+ eventuale centratura
% end
% leggere posizioni picchi QRS detector

BUFF =original_d;
QRSs=QRSs-delay_stadi_QRSdet-total_delay;
QRSs=QRSs*4;
indice_debug=1;
cb=1;
if(length(QRSs)>5)
while(QRSs(cb)<1000)
    cb=cb+1;
end
for idx_QRSs=cb:length(QRSs)
    if QRSs(idx_QRSs)-LEN/2<1
        Kin=1;
    else
        Kin=QRSs(idx_QRSs)-LEN/2;
    end
    if QRSs(idx_QRSs)+LEN/2-1>length(BUFF)
        Kfn=length(BUFF);
    else
        Kfn=QRSs(idx_QRSs)+LEN/2-1;
    end
    if(Kin==1)
%         Kfn=QRSs(idx_QRSs)+LEN-1;
        Kfn=LEN;
    end
    if(Kfn==length(BUFF));
%         Kin=QRSs(idx_QRSs)-LEN+1;
        Kin=length(BUFF)-LEN+1;
    end
    delta=Kfn-Kin;
    prova1=0;
    prova2=1000000;
    for tb=Kin+round(mid/4):Kfn-round(mid/4)
        if(BUFF(tb)>prova1)
            prova1=BUFF(tb);
            m_pos=tb;
        end
        if(BUFF(tb)<prova2)
            prova2=BUFF(tb);
            m_neg=tb;
        end
    end
   prova2=abs(prova2);
%    if(prova2>1.2*prova1)
   if(prova2>1.4*prova1)
       m=m_neg;
   else
       m=m_pos;
   end
%         
%     prova1=max([zeros(1,0.1*LEN) BUFF(Kin+0.1*LEN:Kfn-0.1*LEN) zeros(1,0.1*LEN)]);
%     prova2=min([zeros(1,0.1*LEN) BUFF(Kin+0.1*LEN:Kfn-0.1*LEN) zeros(1,0.1*LEN)]);
%     prova2=abs(prova2);
%     if(prova2>1.2*prova1)
%         [c,m] = (min([zeros(1,0.1*LEN) BUFF(Kin+0.1*LEN:Kfn-0.1*LEN) zeros(1,0.1*LEN)]));
%     else
%         [c,m] = (max([zeros(1,0.1*LEN) BUFF(Kin+0.1*LEN:Kfn-0.1*LEN) zeros(1,0.1*LEN)]));
%     end
%     [c,m] = max(abs([zeros(1,5) BUFF(Kin+5:Kfn-5) zeros(1,5)]));
    tmp = zeros(size(tmp));
    ampt=floor((5/8)*mid);
    sts=floor(mid/8);
    if((m-ampt+2>0)&&(m+ampt+1<=length(BUFF)))        
        tmp(1:2*ampt)=BUFF(m-ampt+2:m+ampt+1);
    end
    invld=0;
    %%analizzo se tmp è centrato male
    if(1.2*max(abs(tmp(ampt-0.1*mid:ampt+0.1*mid)))<(max(abs(tmp))))
        %%template non valido
        invld=1;
    end
    if(invld==0)
    
%     tmp(mid-m:mid-m+delta)=BUFF(Kin:Kfn);
        maxcorr=zeros(NTM,1);
        indexmaxcorr=zeros(NTM,1);
        if Nt == 0; % nessuna classe prima
            A(1,:)=tmp(sts:sts+mid-1);
            Ttot(1,1)=1;
            Nt = 1;
%         s = strcat('Number of templates at now: -', num2str(Nt));
%         s = strcat(s,'-');
%         waitbar(Nt/NTM,hwait2,s);
%         occurrences(i)=Nt;
        else
        % ho almeno un template quindi faccio la cross-correlazione fra tmp e tutti
            for k=1:Nt
%             tmp2=[zeros(1,50) tmp(51:end-50) zeros(1,50)];
                cctemp = pearsoncoeff(tmp,A(k,:));
                cc(k,1:length(cctemp)) = cctemp; %% serve solo per debug
                [maxcorr(k),indexmaxcorr(k)]=max(cctemp);
            end
            [mm,im]= max(maxcorr);
            if mm < Thc
                if Nt < NTM % numero di classi minore di quello massimo: aggiungo una classe
                
                    A(Nt+1,:)=tmp(sts:sts+LEN-1);
                    Ttot(Nt+1,1)=1;
                    Nt=Nt+1;
%                 s = strcat('Number of templates at now: -', num2str(Nt));
%                 s = strcat(s,'-');
%                 waitbar(Nt/NTM,hwait2,s);
%                 occurrences(i)=Nt;
                else % numero di classi pari o maggiore di quello massimo: devo togliere il template con occorrenza minore
                    [mmin, mminidx] = min(Ttot);
%                 s=strcat('Removing template with number of occurrences = ', num2str(mmin));
%                 disp(s);
                    Templ_overwritten = Templ_overwritten + 1;
                % e aggiungere il nuovo
                    A(mminidx,:)=tmp(sts:sts+LEN-1);
                    Ttot(mminidx,1)=1;
%                 occurrences(i)=mminidx+100;
                end
            else
            % tiro via la parte centrale del segnale, allineata sulla correlazione 'estrattoy'
                estrattoy = tmp(indexmaxcorr(im):indexmaxcorr(im)+LEN-1);
%             estrattoy = tmp(mid_mezzi:mid_mezzi+LEN-1);

            % eseguo la media sincronizzata
            
                A(im,:) = synchronizedAVG(A(im,:),estrattoy,Ttot(im,1));
                debug_t(indice_debug,:)=estrattoy;
                vettore_debug(indice_debug)=im; %numero classe
                indice_debug=indice_debug+1;
            
            % aggiorno occorrenze di quella classe
                Ttot(im,1) = Ttot(im,1)+1;
%             occurrences(i)=im;
            end
        end
    end
end
% close(hwait);
% fprintf('-----------------------------------------------------------\n');
% fprintf('AVVIO PROCESSO DI FUSIONE TEMPLATE SIMILI CON SOGLIA %d\n',Thfusione);
indexA = 1;

while indexA < Nt
    tmp(sts:sts+LEN-1) = A(indexA,:);
    maxcorr=zeros(NTM,1);
    indexmaxcorr=zeros(NTM,1);
    for k=indexA+1:Nt
        %plot(A(indexA,:)); hold on; plot(A(k,:),'r');
        cctemp = pearsoncoeff(tmp,A(k,:));
        cc(k,1:length(cctemp)) = cctemp; %% serve solo per debug
        %%permetto solo piccoli spostamenti per allineare
        
        [maxcorr(k),indexmaxcorr(k)]=max(cctemp);
        %close
    end
    [mm,im]= max(maxcorr);
    if mm > Thfusione
%         fprintf('Fuso il template %d (occ. %d) con quello %d (occ. %d)\n',indexA,Ttot(indexA,1),im,Ttot(im,1));
%         fprintf('  - Risultato fusione in pos. %d (occ. %d)\n',indexA,Ttot(indexA,1)+Ttot(im,1));
        estrattoy = tmp(indexmaxcorr(im):indexmaxcorr(im)+LEN-1);
        %figure; plot(estrattoy); hold on; plot(A(im,:),'r')
        A(indexA,:) = synchronizedWAVG(estrattoy,A(im,:),Ttot(indexA,1),Ttot(im,1));
        %plot(A(indexA,:),'g');
        Ttot(indexA,1) = Ttot(indexA,1)+Ttot(im,1);
        A(im,:)=A(Nt,:);
        Ttot(im,1)=Ttot(Nt,1);
        A(Nt,:)=0;
        Ttot(Nt,1)=0;
        Nt=Nt-1;
        indexA=1;
    else
        indexA=indexA+1;
    end
end
% fprintf('Numero di template finali: %d\n', Nt);
% fprintf('Occorrenze per template finali:\n');
% for i=1:Nt
%     fprintf('Template %d: %d\n',i,Ttot(i,1));
% end
Apost = zeros(NumTempl,LEN);

Atemp=A;
Ttottemp = Ttot;
Ttotpost = zeros(NumTempl,1);
% applica riduzione numero di classi
for i = 1 : NumTempl
    [m,idx] = max(Ttottemp);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    Ttotpost(i) = m;
    Apost(i,:)=Atemp(idx,:);
    Ttottemp(idx) = [];
    Atemp(idx,:) = [];
end
clear Atemp
clear Ttottemp
% fprintf('-----------------------------------------------------------\n');
% fprintf('Riduzione numero di template.\n');
% fprintf('Numero finale di template: %d\n', length(Ttotpost));
% fprintf('Numero finale di occorrenze totali: %d\n', sum(Ttotpost));
% fprintf('Occorrenze per template:\n');
for i=1:length(Ttotpost)
%     fprintf('Template %d: %d\n',i,Ttotpost(i,1));
end
% fprintf('-----------------------------------------------------------\n\n\n\n\n');
% close(hwait2);
% Template matching continuo
% if tmatch_on_den~=0
%     %% Template matching sul segnale denoised
%     BUFF = denoised;
% else
%     %% Template matching sul segnale originale senza contare il WD
%     BUFF = signal;
% end
BUFF=original_d;
for k=1:NumTempl  %numtempl è 2
    cctemp=[];
    [mA,iA]=max(Apost(k,:));
    cctemp = pearsoncoeff(BUFF,Apost(k,:));
    cctemp=cctemp.*(abs(cctemp) > Thcross); %hard thresholding della crosscorrelazione
    j=1;
    h=1;
    while j<length(cctemp)
        if cctemp(j)>0
            start_m=j;
            while cctemp(j)>0 && j<length(cctemp)
                end_m=j;
                j = j+1;
            end
            [mp,ip]=max(cctemp(start_m:end_m));
            occ(k,h)=ip+start_m-1;
            occorrenze_finali(k,1)=h;
            h=h+1;
        end
        j=j+1;
    end
end
else
    Ttotpost=0;
    Ttot=0;
    occ=0;
    Apost=0;
    Nt=0;
    occorrenze_finali=0;
end

end

function [r] = pearsoncoeff(x,y)
sx = size(x);
sy = size(y);
if min(size(x)) > 1
    error('Impossibile usare la funzione su dati matriciali!')
    return
end
if min(size(y)) > 1
    error('Impossibile usare la funzione su dati matriciali!')
    return
end
if sx(1) > sx(2)
    x=x';
end
if sy(1) > sy(2)
    y=y';
end
if length(x) < length(y) % il vettore x deve essere il più lungo
    tmp=y;
    y=x;
    x=tmp;
end
ly =length(y);
lr = length(x)-ly+1;
r = zeros(1,lr);
stdy = std(y);
meany = mean(y);
yc = y - meany;
for i=1:lr
   stdx = std(x(i:i+length(y)-1));
   meanx = mean(x(i:i+length(y)-1));
   xc = x(i:i+length(y)-1) - meanx;
   num = sum((1/(ly-1))*xc.*yc);
   r(i)=num/(stdy*stdx);
end
end


function scaling = scal_th(scaling_type, level, buffer_for_threshold, rms_sc)
switch scaling_type
    case 1 
        scaling = median(abs(buffer_for_threshold(:,level))) / 0.6745;
    case 2
        scaling = 1.4826 * quantile(buffer_for_threshold(:,level),.75);
    case 3
        scaling = std(buffer_for_threshold(:,level));
    case 4
        scaling = 1.4826 * median(abs(buffer_for_threshold(:,level)-median(buffer_for_threshold(:,level)))); 
    case 5
        scaling = rms_sc(level) * norm(buffer_for_threshold(:,level))/sqrt(length(buffer_for_threshold(:,level)));
    otherwise scaling = 1;
end
end

function res = standardize(x)
res = (x-mean(x))/std(x);
end

function u = synchronizedAVG(oldu,y,n)
loldu = length(oldu);
ly = length(y);
if loldu ~= ly
    error('Arrays must have the same length');
end
if (min(size(oldu)) ~= 1) || (min(size(y)) ~= 1)
    error('Only arrays, not matrices!');
end
% u=((oldu*n)+y)/(n+1);
if(n>60)
    u=y;
else
    u=((oldu*n)+y)/(n+1);
end
end

function u = synchronizedWAVG(firstsig,secondsig,wfirst,wsecond)
N=wfirst+wsecond;
u=(firstsig*wfirst/N)+(secondsig*wsecond/N);
end

function [y]=upsampling(x)
l=length(x);
n=2*l;
y=zeros(1,n);
for u=1:l
    y(2*u-1)=x(u);
end
end

function [soglia]=sogliatrovapiccoR_offline(segnale,G)
soglia=[];
i=1;
while(i+1024)<length(segnale)
    if((i-1024)>0)
%         B=sort(segnale(i-1023:i+1024));
%         r=G*rms(B(10:end-10));
        r=G*mean(segnale(i-1023:i+1024));
    else
%         B=sort(segnale(i:i+2*1024-1));
%         r=G*rms(B(10:end-10));
    r=G*mean(segnale(i:i+2*1024-1));
    end
    if r<0.02;
        r=0.02;
    end
    soglia(i:i+1024)=r;
    i=i+256;
end
soglia(i:length(segnale))=r;
end

function [yf]=media_mobile(x,N)
b=(1/N)* ones(1,N);
a=1;
yf=filter(b,a,x);
end

function [chann] = ChannelSelectionOrCombination(num_occ1,num_occ2, num_occ3,num_occ4)
    A=[num_occ1(2);num_occ2(2); num_occ3(2); num_occ4(2)];
%  A=[num_occ1(1);num_occ2(1); num_occ3(1); num_occ4(1)];
    [k,chann]=max(A);
end

function r=rms(x)

r=norm(x)/sqrt(length(x));

end

function [FQRS]=correction(SelectedResidual,FQRS)
%@1000Hz
R=10*2; %numero pari
for r=1:length(FQRS)
    in=FQRS(r);
    max_loc=[];
    ind_max=[];
    min_loc=[];
    ind_min=[];
    mx=1;
    mi=1;
    if(((in+R)<length(SelectedResidual))&&((in-R)>0))
        for hn=0:(2*R-4)  %2*R+1 elementi
            %ricerca massimi locali
            if((SelectedResidual(in-R+hn)<=SelectedResidual(in-R+hn+1))&&(SelectedResidual(in-R+hn+1)<=SelectedResidual(in-R+hn+2))...
                    &&(SelectedResidual(in-R+hn+3)<=SelectedResidual(in-R+hn+2))&&(SelectedResidual(in-R+hn+4)<=SelectedResidual(in-R+hn+3)))
                max_loc(mx)=SelectedResidual(in-R+hn+2);
                ind_max(mx)=in-R+hn+2;
                mx=mx+1;
            end
            if((SelectedResidual(in-R+hn)>=SelectedResidual(in-R+hn+1))&&(SelectedResidual(in-R+hn+1)>=SelectedResidual(in-R+hn+2))...
                    &&(SelectedResidual(in-R+hn+3)>=SelectedResidual(in-R+hn+2))&&(SelectedResidual(in-R+hn+4)>=SelectedResidual(in-R+hn+3)))
                min_loc(mx)=SelectedResidual(in-R+hn+2);
                ind_min(mx)=in-R+hn+2;
                mi=mi+1;
            end
        end
    end
    %condizione vettore vuoto (analizzare morfologia picco  qrs fetale
    %rispetto a picco di rumore: dovrebbe essere più ampio! quindi analisi
    %su più campioni invece che solo 2 a dx e 2 a sx?
    MM=0;
    MI=0;
    if(isempty(max_loc))
    else
    [MM,MI]=max(abs(max_loc));
    end
    mm=0;
    mi=0;
    if(isempty(min_loc))
    else
    [mm,mi]=max(abs(min_loc));
    end
    if((MM>mm)&&(MM>0))
        FQRS(r)=ind_max(MI);
    elseif((mm>MM)&&(mm>0))
        FQRS(r)=ind_min(mi);
    end
end

end

function [c, vett_soglia]=qrs_detector_wd(x)
    [c]=estrazioneFS(x);
    [vett_soglia]=sogliatrovapiccoR_offline(c,0.8);
end

function [c]=estrazioneFS(x)
    y0=derivatore_elementare(x);
    y0=abs(y0);
    b1 = [1 0 -2 0 1];
    a1 = 1;
    y1 = filter(b1,a1,x);
    y1 = abs(y1);
    y0 = [0,y0(1:end-1)]; 
    y2 = 1.3*y0 + 1.1*y1; 
    c=media_mobile(y2,10);
end

function [yf]=derivatore_elementare(x)
yf=filter([.5 -.5],1,x);
end

% 
% function[occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_v2(occ,stop_cond,thRR,gamma,starting_delay )
%     occtmp=0;
%     occchk=0;
%     RRm=0;
%     count_removed=0;
%     count_added=0;
%     if nargin<4
%         thRR=20;
%         gamma=1/4;
%     end
%     noise=0;
%     if size(occ,1)>size(occ,2)
%         occ=occ';
%     end
%     occtmp=occ(1,occ>0); %taglio la coda di eventuali zeri
%     RRist=occtmp(2:end)-occtmp(1:end-1); %trovo gli RR istantanei
%     RRm=round(median(RRist)); %mediana RR starting value per autoregressivo
%     %RRm=100;%%%%%%%%%%%%%%%%%%%%%%%%%Solo per test.
%     delta=abs(RRist-RRm); %identifico lo scostamento di ogni distanza RR dal valore mediano
%     occchk=zeros(size(delta)); %creo vettore di check su tale scostamento(0=fuori tolleranza)
%     occchk(delta<=thRR)=1; %come sopra (1=ok)
%     start=find(occchk, 1, 'first'); %prendo il primo elemento ok e parto dalì
%     if start+2 < length(occchk)
%         while occchk(start)+occchk(start+1)+occchk(start+2) ~= 3
%             start = start+1;
%             if start+2 >= length(occchk)
%                 noise = 1;
%                 break;
%             end
%         end
%     end
%     start = start + 1;
% 
% 
%     %prima ragioniamo in avanti (da start in poi), poi, se necessario,andremo
%     %ai punti che precedono start.
%     if noise == 0
%         count_removed = 0;
%         count_added = 0;
%         i=start;
%         RRm = occtmp(i+1)-occtmp(i);
%         while i<length(occtmp)-2 %fino alla fine del vettore
%           
% %             stem(occtmp(1:end-1),occchk,'r')
%             if occchk(i) ~= 0 %se la distanza fra l'i+1 e l'i è ok per qualche motivo...
%                 %non può essercene uno migliore intorno solo all'inizio, ma a
%                 %regime sì. Quindi...
%                 %provo a vedere con l'autoregressivo dove dovrebbe trovarsi la prossima occorrenza
%                 nextR = occtmp(i) + RRm;
% 
%                 %provo a mettermi nell'intorno di tolleranza sulla posizione
%                 %stimata del prossimo picco per "fare pulizia" se è il caso
%                 j=i+1;
%                 best_occ = [];
%                 best_idx = [];
%                 while abs(occtmp(j)-nextR) <=  thRR
%                     if isempty(best_occ) %inizializzazione
%                         best_occ = abs(occtmp(j)-nextR);
%                         best_idx = j;
%                     else %se c'è più di un picco nella fascia di tolleranza
%                         if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
%                             best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
%                             best_idx = j; %conservo l'indice
%                         end
%                     end
%                     j=j+1;
%                     if j>length(occtmp)
%                         break
%                     end
%                 end
%                 %cancello tutti i picchi nell'intorno del migliore (uno almeno
%                 %lo trovo, sennò non ci sarebbe stato un occchk=1
%                 for aux=i+1:j-1
%                     if aux ~= best_idx
%                         occtmp(aux)=60000; %marco per poi rimuovere
%                         occchk(aux)=60000; %marco per poi rimuovere
%                         count_removed = count_removed+1; %conto i rimossi
%                     end
%                 end
%                 occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                 occchk(occchk==60000)=[]; %rimuovo quelli marcati
%                 %sistemo i vettori occchk e occtmp
%                 if(i<length(occtmp)-2)
%                     break;
%                 end
%                 if j-1>=i+2 %se ne ho trovato solo uno lo lascio indicato com'è
%                     if abs(occtmp(i+2)-occtmp(i+1)-RRm) < thRR
%                         occchk(i+1)=1;
%                     else
%                         occchk(i+1)=0;
%                     end
%                 end
% 
%                 %aggiorno l'autoregressivo
%                 RRm=round(RRm*(1-gamma)+(occtmp(i+1)-occtmp(i))*gamma);
%             else %se la distanza fra l'i+1 e l'i è fuori tolleranza...
%                 %provo a vedere con l'autoregressivo dove dovrebbe trovarsi la prossima occorrenza
%                 nextR = occtmp(i) + RRm;
% 
%                 %quindi tolgo tutto quello fra il primo picco "buono" e il margine
%                 %inferiore di tolleranza rispetto a dove dovrebbe essere il
%                 %prossimo picco
%                 j=i+1; %mi posiziono sulla prossima occorrenza
%                 while occtmp(j) < nextR - thRR %fino a che la distanza fra due occorrenze non supera la zona "buona" intorno al prossimo picco presunto
%                     occtmp(j)=60000; %marco per poi rimuovere
%                     occchk(j)=60000; %marco per poi rimuovere
%                     count_removed = count_removed+1; %conto i rimossi
%                     j=j+1; %provo ad andre avanti di uno
%                     if j>length(occtmp)
%                         break
%                     end
%                 end
%                 occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                 occchk(occchk==60000)=[]; %rimuovo quelli marcati
%                 if abs(occtmp(i+1)-occtmp(i)-RRm) < thRR
% %                 if abs(occtmp(i+1)-occtmp(i)-nextR) < thRR
%                     occchk(i)=2;
%                 else
%                     occchk(i)=0;
%                 end
% 
%                 nextR = occtmp(i) +RRm;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                 %provo a mettermi nell'intorno di tolleranza sulla posizione
%                 %stimata del prossimo picco per "fare pulizia" se è il caso
%                 j=i+1;
%                 best_occ = [];
%                 best_idx = [];
%                 while abs(occtmp(j)-nextR) <=  thRR
%                     if isempty(best_occ) %inizializzazione
%                         best_occ = abs(occtmp(j)-nextR);
%                         best_idx = j;
%                     else %se c'è più di un picco nella fascia di tolleranza
%                         if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
%                             best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
%                             best_idx = j; %conservo l'indice
%                         end
%                     end
%                     j=j+1;
%                     if j>length(occtmp)
%                         break
%                     end
%                 end
%                 if isempty(best_occ) %vuol dire che non c'erano elementi nella zona intorno alla posizione presunta
%                     %devo aggiungere un picco nella posizione presunta
%                     occtmp=[occtmp(1:i),nextR,occtmp(i+1:end)];
%                     occchk=[occchk(1:i),3,occchk(i+1:end)]; %col codice 3 identifico un picco che è stato aggiunto
%                     count_added = count_added+1;
%                     if abs(occtmp(i+2)-occtmp(i+1)-RRm) < thRR
%                         occchk(i+1)=3;
%                     else
%                         occchk(i+1)=0;
%                     end
%                     occchk(i)=3;
%                 else %ho trovato un picco molto prossimo alla posizione presunta
%                     %cancello tutti gli altri trovati
%                     for aux=i+1:j-1
%                         if aux ~= best_idx
%                             occtmp(aux)=60000; %marco per poi rimuovere
%                             occchk(aux)=60000; %marco per poi rimuovere
%                             count_removed = count_removed+1; %conto i rimossi
%                         end
%                     end
%                     occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                     occchk(occchk==60000)=[]; %rimuovo quelli marcati
%                     if(i<length(occtmp)-2)
%                         break;
%                     end
%                     %sistemo i vettori occchk e occtmp
%                     if j-1>=i+2 %se ne ho trovato solo uno lo lascio indicato com'è
%                         if abs(occtmp(i+2)-occtmp(i+1)-RRm) < thRR
%                             occchk(i+1)=1;
%                         else
%                             occchk(i+1)=0;
%                         end
%                     end
%                 end
%             end
%             i=i+1;
%         end
%         %coda se non sono state trovate occorrenze fino alla fine del segnale
%         while stop_cond - occtmp(length(occtmp)) > RRm
%             nextR = occtmp(length(occtmp)) + RRm;
%             occtmp=[occtmp,nextR];
%             occchk=[occchk,3]; %col codice 3 identifico un picco che è stato aggiunto
%             count_added = count_added+1;
%         end
% 
% 
%         % adesso torniamo indietro
%         i=start;
%         RRm = occtmp(i+1)-occtmp(i);
%         %RRm=round(median(RRist)); %mediana RR starting value per autoregressivo
%         %RRm=100;
%         while i>3 %fino alla testa del vettore (corretto! se metto 1 da errore, indici negativi)
%             
% %             stem(occtmp(1:end-1),occchk,'r')
%             if occchk(i-1) ~= 0 %se la distanza fra l'i-1 e l'i è ok per qualche motivo...
%                 %non può essercene uno migliore intorno solo all'inizio, ma a
%                 %regime sì. Quindi...
%                 %provo a vedere con l'autoregressivo dove dovrebbe trovarsi la precedente occorrenza
%                 nextR = occtmp(i) - RRm;
%                 if(nextR<starting_delay)
%                     break
%                 end
% 
%                 %provo a mettermi nell'intorno di tolleranza sulla posizione
%                 %stimata del prossimo picco per "fare pulizia" se è il caso
%                 j=i-1;
%                 best_occ = [];
%                 best_idx = [];
%                 while abs(occtmp(j)-nextR) <=  thRR
%                     if isempty(best_occ) %inizializzazione
%                         best_occ = abs(occtmp(j)-nextR);
%                         best_idx = j;
%                     else %se c'è più di un picco nella fascia di tolleranza
%                         if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
%                             best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
%                             best_idx = j; %conservo l'indice
%                         end
%                     end
%                     j=j-1;
%                     if j<1
%                         break
%                     end
%                 end
%                 %cancello tutti i picchi nell'intorno del migliore (uno almeno
%                 %lo trovo, sennò non ci sarebbe stato un occchk=1
%                 for aux=i-1:-1:j+1
%                     if aux ~= best_idx
%                         occtmp(aux)=60000; %marco per poi rimuovere
%                         occchk(aux)=60000; %marco per poi rimuovere
%                         count_removed = count_removed+1; %conto i rimossi
%                     end
%                 end
%                 occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                 occchk(occchk==60000)=[]; %rimuovo quelli marcati
%                 %sistemo i vettori occchk e occtmp
%                 if j+1<i-1 %se ne ho trovato solo uno lo lascio indicato com'è
%                     if abs(occtmp(i-1)-occtmp(i-2)-RRm) < thRR
%                         occchk(i-2)=1;
%                     else
%                         occchk(i-2)=0;
%                     end
%                 end
% 
%                 %aggiorno l'autoregressivo
%                 RRm=round(RRm*(1-gamma)+(occtmp(i)-occtmp(i-1))*gamma);
%             else %se la distanza fra l'i e l'i-1 è fuori tolleranza...
%                 %provo a vedere con l'autoregressivo dove dovrebbe trovarsi
%                 %la prossima occorrenza
%                 nextR = occtmp(i) - RRm;
% 
%                 %quindi tolgo tutto quello fra il primo picco "buono" e il margine
%                 %inferiore di tolleranza rispetto a dove dovrebbe essere il
%                 %prossimo picco
%                 j=i-1; %mi posiziono sulla prossima occorrenza
%                 while occtmp(j) >= nextR + thRR %fino a che la distanza fra due occorrenze non supera la zona "buona" intorno al prossimo picco presunto
%                     occtmp(j)=60000; %marco per poi rimuovere
%                     occchk(j)=60000; %marco per poi rimuovere
%                     count_removed = count_removed+1; %conto i rimossi
%                     j=j-1; %provo ad andre avanti di uno
%                     if j<1
%                         break
%                     end
%                 end
%                 if abs(occtmp(i)-occtmp(j)-RRm) <thRR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     occchk(j)=2;
%                 else
%                     occchk(j)=0;
%                 end
%                 occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                 occchk(occchk==60000)=[]; %rimuovo quelli marcati
% 
%                 nextR = occtmp(i) - RRm;
% 
%                 %provo a mettermi nell'intorno di tolleranza sulla posizione
%                 %stimata del prossimo picco per "fare pulizia" se è il caso
%                 j=i-1;
%                 best_occ = [];
%                 best_idx = [];
%                 while abs(occtmp(j)-nextR) <=  thRR
%                     if isempty(best_occ) %inizializzazione
%                         best_occ = abs(occtmp(j)-nextR);
%                         best_idx = j;
%                     else %se c'è più di un picco nella fascia di tolleranza
%                         if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
%                             best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
%                             best_idx = j; %conservo l'indice
%                         end
%                     end
%                     j=j-1;
%                     if j<1
%                         break
%                     end
%                 end
%                 if isempty(best_occ) %vuol dire che non c'erano elementi nella zona intorno alla posizione presunta
%                     %devo aggiungere un picco nella posizione presunta
%                     occtmp=[occtmp(1:i-1),nextR,occtmp(i:end)];
%                     occchk=[occchk(1:i-1),3,occchk(i:end)]; %col codice 3 identifico un picco che è stato aggiunto%%%%%%%%%%%%%%
%                     count_added = count_added+1;
%                     if abs(occtmp(i)-occtmp(i-1)-RRm) < thRR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         occchk(i-1)=3;
%                     else
%                         occchk(i-1)=0;
%                     end
%                     i=i+1;
%                 else %ho trovato un picco molto prossimo alla posizione presunta
%                     %cancello tutti gli altri trovati
%                     for aux=i-1:-1:j+1
%                         if aux ~= best_idx
%                             occtmp(aux)=60000; %marco per poi rimuovere
%                             occchk(aux)=60000; %marco per poi rimuovere
%                             count_removed = count_removed+1; %conto i rimossi
%                         end
%                     end
%                     occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                     occchk(occchk==60000)=[]; %rimuovo quelli marcati
%                     %sistemo i vettori occchk e occtmp
%                     if j+1<i-1
%                         i=j+2;
%                         if abs(occtmp(i)-occtmp(i-1)-RRm) < thRR
%                             occchk(i-1)=1;
%                         else
%                             occchk(i-1)=0;
%                         end
%                     end
%                 end
%             end
%             i=i-1;
%         end
%         %coda se non inizia subito il vettore delle occorrenze
%         while occtmp(1) > RRm
%             nextR = occtmp(1) - RRm;
%             if(nextR<starting_delay)
%                 break
%             end
%             occtmp=[nextR,occtmp];
%             occchk=[3,occchk]; %col codice 3 identifico un picco che è statoaggiunto
%             count_added = count_added+1;
%         end
% 
%     end
% end


% 
% function[occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_v2(occ,stop_cond,thRR,gamma,starting_delay)
%     occtmp=0;
%     occchk=0;
%     RRm=0;
%     count_removed=0;
%     count_added=0;
%     noise=0;
%     if size(occ,1)>size(occ,2)
%         occ=occ';
%     end
%     occtmp=occ(1,occ>0); %taglio la coda di eventuali zeri
%     RRist=occtmp(2:end)-occtmp(1:end-1); %trovo gli RR istantanei
%     RRm=round(median(RRist)); %mediana RR starting value per autoregressivo
%     %RRm=100;%%%%%%%%%%%%%%%%%%%%%%%%%Solo per test.
%     delta=abs(RRist-RRm); %identifico lo scostamento di ogni distanza RR dal valore mediano
%     occchk=zeros(size(delta)); %creo vettore di check su tale     scostamento(0=fuori tolleranza)
%     occchk(delta<=thRR)=1; %come sopra (1=ok)
%     start=find(occchk, 1, 'first'); %prendo il primo elemento ok e parto dalì
%     if start+2 < length(occchk)
%         while occchk(start)+occchk(start+1)+occchk(start+2) ~= 3
%             start = start+1;
%             if start+2 >= length(occchk)
%                 noise = 1;
%                 break;
%             end
%         end
%     end
% 
%     start = start+1;
% 
%     %prima ragioniamo in avanti (da start in poi), poi, se necessario,andremo
%     %ai punti che precedono start.
%     if noise == 0
%         count_removed = 0;
%         count_added = 0;
%         i=start;
%         RRm = occtmp(i+1)-occtmp(i);
%         while i<length(occtmp)-2 %fino alla fine del vettore
% 
%             if occchk(i) ~= 0 %se la distanza fra l'i+1 e l'i è ok per qualche motivo...
%                 %non può essercene uno migliore intorno solo all'inizio, ma a
%                 %regime sì. Quindi...
%                 %provo a vedere con l'autoregressivo dove dovrebbe trovarsi la prossima occorrenza
%                 nextR = occtmp(i) + RRm;
% 
%                 %provo a mettermi nell'intorno di tolleranza sulla posizione
%                 %stimata del prossimo picco per "fare pulizia" se è il caso
%                 j=i+1;
%                 best_occ = [];
%                 best_idx = [];
%                 while abs(occtmp(j)-nextR) <= thRR
%                     if isempty(best_occ) %inizializzazione
%                         best_occ = abs(occtmp(j)-nextR);
%                         best_idx = j;
%                     else %se c'è più di un picco nella fascia di tolleranza
%                         if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
%                             best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
%                             best_idx = j; %conservo l'indice
%                         end
%                     end
%                     j=j+1;
%                     if j>length(occtmp)
%                         break
%                     end
%                 end
%             %cancello tutti i picchi nell'intorno del migliore (uno almeno
% %lo trovo, sennò non ci sarebbe stato un occchk=1
%                 for aux=i+1:j-1
%                     if aux ~= best_idx
%                         occtmp(aux)=60000; %marco per poi rimuovere
%                         occchk(aux)=60000; %marco per poi rimuovere
%                         count_removed = count_removed+1; %conto i rimossi
%                     end
%                 end
%                 occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                 occchk(occchk==60000)=[]; %rimuovo quelli marcati
%                 %sistemo i vettori occchk e occtmp
%                 if(i>length(occtmp)-2)
%                     break;
%                 end
%                 if j-1>=i+2 %se ne ho trovato solo uno lo lascio indicato com'è
%                     if abs(occtmp(i+2)-occtmp(i+1)-RRm) <= thRR
%                         occchk(i+1)=1;
%                     else
%                         occchk(i+1)=0;
%                     end
%                 end
% 
%                 %aggiorno l'autoregressivo
%                 RRm=round(RRm*(1-gamma)+(occtmp(i+1)-occtmp(i))*gamma);
%             else %se la distanza fra l'i+1 e l'i è fuori tolleranza...
%             %provo a vedere con l'autoregressivo dove dovrebbe trovarsi  la prossima occorrenza
%                 nextR = occtmp(i) + RRm;
% 
%                 %quindi tolgo tutto quello fra il primo picco "buono" e il margine
%                 %       inferiore di tolleranza rispetto a dove dovrebbe essere il
%                 %prossimo picco
%                 j=i+1; %mi posiziono sulla prossima occorrenza
%                 while occtmp(j) < nextR - thRR %fino a che la distanza fra due occorrenze non supera la zona "buona" intorno al prossimo picco presunto
%                     occtmp(j)=60000; %marco per poi rimuovere
%                     occchk(j)=60000; %marco per poi rimuovere
%                     count_removed = count_removed+1; %conto i rimossi
%                     j=j+1; %provo ad andre avanti di uno
%                     if j>length(occtmp)
%                         break
%                     end
%                 end
%                 occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                 occchk(occchk==60000)=[]; %rimuovo quelli marcati
%                 if abs(occtmp(i+1)-occtmp(i)-RRm) <= thRR
% % if abs(occtmp(i+1)-occtmp(i)-nextR) < thRR
%                     occchk(i)=2;
%                 else
%                     occchk(i)=0;
%                 end
% 
%                 nextR = occtmp(i) +RRm;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %provo a mettermi nell'intorno di tolleranza sulla posizione
% %stimata del prossimo picco per "fare pulizia" se è il caso
%                 j=i+1;
%                 best_occ = [];
%                 best_idx = [];
%                 while abs(occtmp(j)-nextR) <= thRR
%                     if isempty(best_occ) %inizializzazione
%                         best_occ = abs(occtmp(j)-nextR);
%                         best_idx = j;
%                     else %se c'è più di un picco nella fascia di tolleranza
%                         if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
%                             best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
%                             best_idx = j; %conservo l'indice
%                         end
%                     end
%                     j=j+1;
%                     if j>length(occtmp)
%                     break
%                     end
%                 end
%                 if isempty(best_occ) %vuol dire che non c'erano elementi nella zona intorno alla posizione presunta
% %devo aggiungere un picco nella posizione presunta
%                     occtmp=[occtmp(1:i),nextR,occtmp(i+1:end)];
%                     occchk=[occchk(1:i),3,occchk(i+1:end)]; %col codice 3 identifico un picco che è stato aggiunto
%                     count_added = count_added+1;
%                     if abs(occtmp(i+2)-occtmp(i+1)-RRm) <= thRR
%                         occchk(i+1)=3;
%                     else
%                         occchk(i+1)=0;
%                     end
%                     occchk(i)=3;
%                     else %ho trovato un picco molto prossimo alla posizione presunta
% %cancello tutti gli altri trovati
%                     for aux=i+1:j-1
%                         if aux ~= best_idx
%                             occtmp(aux)=60000; %marco per poi rimuovere
%                             occchk(aux)=60000; %marco per poi rimuovere
%                             count_removed = count_removed+1; %conto i rimossi
%                         end
%                     end
%                     occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                     occchk(occchk==60000)=[]; %rimuovo quelli marcati
%                     if(i>length(occtmp)-2)
%                         break;
%                     end
%                     %sistemo i vettori occchk e occtmp
%                     if j-1>=i+2 %se ne ho trovato solo uno lo lascio indicato com'è
%                         if abs(occtmp(i+2)-occtmp(i+1)-RRm) <= thRR
%                             occchk(i+1)=1;
%                         else
%                             occchk(i+1)=0;
%                         end
%                     end
%                 end
%             end
%             i=i+1;
%         end
%         %coda se non sono state trovate occorrenze fino alla fine del  segnale
%         while stop_cond - occtmp(length(occtmp)) > RRm
%             nextR = occtmp(length(occtmp)) + RRm;
%             occtmp=[occtmp,nextR];
%             occchk=[occchk,3]; %col codice 3 identifico un picco che è stato  aggiunto
%             count_added = count_added+1;
%         end
% 
% 
%         % adesso torniamo indietro
%         i=start;
%         RRm = occtmp(i+1)-occtmp(i);
%         %RRm=round(median(RRist)); %mediana RR starting value per autoregressivo
%         %RRm=100;
%         while i>3 %fino alla testa del vettore (corretto! se metto 1 da errore, indici negativi)
% 
%             if occchk(i-1) ~= 0 %se la distanza fra l'i-1 e l'i è ok per qualche motivo...
%                 %non può essercene uno migliore intorno solo all'inizio, ma a
%                 %regime sì. Quindi...
%                 %provo a vedere con l'autoregressivo dove dovrebbe trovarsi la precedente occorrenza
%                 nextR = occtmp(i) - RRm;
%                 if(nextR<starting_delay)
%                     break
%                 end
% 
%                 %provo a mettermi nell'intorno di tolleranza sulla posizione
%                 %stimata del prossimo picco per "fare pulizia" se è il caso
%                 j=i-1;
%                 best_occ = [];
%                 best_idx = [];
%                 while abs(occtmp(j)-nextR) <= thRR
%                     if isempty(best_occ) %inizializzazione
%                         best_occ = abs(occtmp(j)-nextR);
%                         best_idx = j;
%                     else %se c'è più di un picco nella fascia di tolleranza
%                         if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
%                             best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
%                             best_idx = j; %conservo l'indice
%                         end
%                     end
%                     j=j-1;
%                     if j<1
%                         break
%                     end
%                 end
%                 %cancello tutti i picchi nell'intorno del migliore (uno almeno
%                 %lo trovo, sennò non ci sarebbe stato un occchk=1
%                 for aux=i-1:-1:j+1
%                     if aux ~= best_idx
%                         occtmp(aux)=60000; %marco per poi rimuovere
%                         occchk(aux)=60000; %marco per poi rimuovere
%                         count_removed = count_removed+1; %conto i rimossi
%                     end
%                 end
%                 occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                 occchk(occchk==60000)=[]; %rimuovo quelli marcati
%                 %sistemo i vettori occchk e occtmp
%                 if j+1<i-1 %se ne ho trovato solo uno lo lascio indicato com'è
%                     if abs(occtmp(i-1)-occtmp(i-2)-RRm) <= thRR
%                         occchk(i-2)=1;
%                     else
%                         occchk(i-2)=0;
%                     end
%                 end
% 
%                 %aggiorno l'autoregressivo
%                 RRm=round(RRm*(1-gamma)+(occtmp(i)-occtmp(i-1))*gamma);
%             else %se la distanza fra l'i e l'i-1 è fuori tolleranza...
%                 %provo a vedere con l'autoregressivo dove dovrebbe trovarsi
%                 %la prossima occorrenza
%                 nextR = occtmp(i) - RRm;
% 
%                 %quindi tolgo tutto quello fra il primo picco "buono" e il margine
%                 %inferiore di tolleranza rispetto a dove dovrebbe essere il
%                 %prossimo picco
%                 j=i-1; %mi posiziono sulla prossima occorrenza
%                 while occtmp(j) >= nextR + thRR %fino a che la distanza fra due occorrenze non supera la zona "buona" intorno al prossimo picco presunto
%                     occtmp(j)=60000; %marco per poi rimuovere
%                     occchk(j)=60000; %marco per poi rimuovere
%                     count_removed = count_removed+1; %conto i rimossi
%                     j=j-1; %provo ad andre avanti di uno
%                     if j<1
%                         break
%                     end
%                 end
%                 if abs(occtmp(i)-occtmp(j)-RRm) <= thRR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     occchk(j)=2;
%                 else
%                     occchk(j)=0;
%                 end
%                 occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                 occchk(occchk==60000)=[]; %rimuovo quelli marcati
% 
%                 nextR = occtmp(i) - RRm;
% 
%                 %provo a mettermi nell'intorno di tolleranza sulla posizione
%                 %stimata del prossimo picco per "fare pulizia" se è il caso
%                 j=i-1;
%                 best_occ = [];
%                 best_idx = [];
%                 while abs(occtmp(j)-nextR) <= thRR
%                     if isempty(best_occ) %inizializzazione
%                         best_occ = abs(occtmp(j)-nextR);
%                         best_idx = j;
%                     else %se c'è più di un picco nella fascia di tolleranza
%                         if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
%                             best_occ = abs(occtmp(j)-nextR); %questo diventa                            il migliore
%                             best_idx = j; %conservo l'indice
%                         end
%                     end
%                     j=j-1;
%                     if j<1
%                         break
%                     end
%                 end
%                 if isempty(best_occ) %vuol dire che non c'erano elementi nella zona intorno alla posizione presunta
%                     %devo aggiungere un picco nella posizione presunta
%                     occtmp=[occtmp(1:i-1),nextR,occtmp(i:end)];
%                     occchk=[occchk(1:i-1),3,occchk(i:end)]; %col codice 3identifico un picco che è stato aggiunto%%%%%%%%%%%%%%
%                     count_added = count_added+1;
%                     if abs(occtmp(i)-occtmp(i-1)-RRm) <= thRR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         occchk(i-1)=3;
%                     else
%                         occchk(i-1)=0;
%                     end
%                     i=i+1;
%                 else %ho trovato un picco molto prossimo alla posizione presunta
%                     %cancello tutti gli altri trovati
%                     for aux=i-1:-1:j+1
%                         if aux ~= best_idx
%                             occtmp(aux)=60000; %marco per poi rimuovere
%                             occchk(aux)=60000; %marco per poi rimuovere
%                             count_removed = count_removed+1; %conto i rimossi
%                         end
%                     end
%                     occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
%                     occchk(occchk==60000)=[]; %rimuovo quelli marcati
%                     %sistemo i vettori occchk e occtmp
%                     if j+1<i-1
%                         i=j+2;
%                         if abs(occtmp(i)-occtmp(i-1)-RRm) <= thRR
%                             occchk(i-1)=1;
%                         else
%                             occchk(i-1)=0;
%                         end
%                     end
%                 end
%             end
%             i=i-1;
%         end
%         %coda se non inizia subito il vettore delle occorrenze
%         while occtmp(1) > RRm
%             nextR = occtmp(1) - RRm;
%             if(nextR<starting_delay)
%                 break
%             end
%             occtmp=[nextR,occtmp];
%             occchk=[3,occchk]; %col codice 3 identifico un picco che è statoaggiunto
%             count_added = count_added+1;
%         end
% 
%     end
% end
% 
% 

function[occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_cor_caller(ecg_ref,mark,occ,stop_cond,thRR,gamma,starting_delay)

    [preval, idx_preval, len_longest, type] = prevalence(mark);
    if preval == 1
        other = 2;
    else
        other = 1;
    end

    switch type
        case 1
            % c'è solo una classe
            [occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_2classes(ecg_ref,...
                mark,occ,stop_cond,thRR,gamma,starting_delay,preval,idx_preval,len_longest);
        case 20
            % situazione generica a due classi
            [occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_2classes(ecg_ref,...
                mark,occ,stop_cond,thRR,gamma,starting_delay,preval,idx_preval,len_longest);
        case 21
            % c'è una classe molto più rappresentata dell'altra
            i = find(mark==other);
            mark_preval = mark;
            occ_preval = occ;
            mark_preval(i) = [];
            occ_preval(i) = [];
            [idx_preval, length_longest_preval] =longest_seq(mark_preval,preval);
            [occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_2classes(ecg_ref,...
                mark_preval,occ_preval,stop_cond,thRR,gamma,starting_delay,preval,idx_preval,length_longest_preval);
        case 22
            % siamo in condizioni di alternanza
            % pulisce delle occorrenze dell'altro template
            i = find(mark==other);
            mark_preval = mark;
            occ_preval = occ;
            mark_preval(i) = [];
            occ_preval(i) = [];
            [idx_preval, length_longest_preval] =longest_seq(mark_preval,preval);
            [occtmp1,occchk1,RRm1,count_removed1,count_added1,noise1]=periodicity_correction_2classes(ecg_ref,...
                mark_preval,occ_preval,stop_cond,thRR,gamma,starting_delay,preval,idx_preval,length_longest_preval);


            % pulisce delle occorrenze dell'altro template
            i = find(mark==preval);
            mark_other = mark;
            occ_other = occ;
            mark_other(i) = [];
            occ_other(i) = [];
            [idx_other, length_longest_other] = longest_seq(mark_other,other);
            [occtmp2,occchk2,RRm2,count_removed2,count_added2,noise2]=periodicity_correction_2classes(ecg_ref,...
                mark_other,occ_other,stop_cond,thRR,gamma,starting_delay,preval,idx_other,length_longest_other);

            % confronta le sequenze per capire quale è materna
            if RRm1>RRm2
                occtmp = occtmp1;
                occchk = occchk1;
                RRm = RRm1;
                count_removed = count_removed1;
                count_added = count_added1;
                noise = noise1;
            else
                occtmp = occtmp2;
                occchk = occchk2;
                RRm = RRm2;
                count_removed = count_removed2;
                count_added = count_added2;
                noise = noise2;
            end
        otherwise % type = 0 (non ci sono proprio template riconosciuti)
            [occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_v2(ecg_ref,occ,stop_cond,thRR,gamma,starting_delay,1,1);
    end
end


function[occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_cor_callerFETUS(ecg_ref,mark,occ,stop_cond,thRR,gamma,starting_delay)

    [preval, idx_preval, len_longest, type] = prevalence(mark);
    if preval == 1
        other = 2;
    else
        other = 1;
    end

    switch type
        case 1
            % c'è solo una classe
            [occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_2classes(ecg_ref,...
                mark,occ,stop_cond,thRR,gamma,starting_delay,preval,idx_preval,len_longest);
        case 20
            % situazione generica a due classi
            [occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_2classes(ecg_ref,...
                mark,occ,stop_cond,thRR,gamma,starting_delay,preval,idx_preval,len_longest);
        case 21
            % c'è una classe molto più rappresentata dell'altra
            i = find(mark==other);
            mark_preval = mark;
            occ_preval = occ;
            mark_preval(i) = 0;
            %occ_preval(i) = [];
            [idx_preval, length_longest_preval] =longest_seq(mark_preval,preval);
            [occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_2classes(ecg_ref,...
                mark_preval,occ_preval,stop_cond,thRR,gamma,starting_delay,preval,idx_preval,length_longest_preval);
        case 22
            % siamo in condizioni di alternanza
            % pulisce delle occorrenze dell'altro template
            i = find(mark==other);
            mark_preval = mark;
            occ_preval = occ;
            mark_preval(i) = 0;
            %occ_preval(i) = [];
            [idx_preval, length_longest_preval] =longest_seq(mark_preval,preval);
            [occtmp1,occchk1,RRm1,count_removed1,count_added1,noise1]=periodicity_correction_2classes(ecg_ref,...
                mark_preval,occ_preval,stop_cond,thRR,gamma,starting_delay,preval,idx_preval,length_longest_preval);


            % pulisce delle occorrenze dell'altro template
            i = find(mark==preval);
            mark_other = mark;
            occ_other = occ;
            mark_other(i) = 0;
            %occ_other(i) = [];
            [idx_other, length_longest_other] = longest_seq(mark_other,other);
            [occtmp2,occchk2,RRm2,count_removed2,count_added2,noise2]=periodicity_correction_2classes(ecg_ref,...
                mark_other,occ_other,stop_cond,thRR,gamma,starting_delay,preval,idx_other,length_longest_other);

            % confronta le sequenze per capire quale è fetale
            if RRm1>RRm2
                occtmp = occtmp2;
                occchk = occchk2;
                RRm = RRm2;
                count_removed = count_removed2;
                count_added = count_added2;
                noise = noise2;
            else
                occtmp = occtmp1;
                occchk = occchk1;
                RRm = RRm1;
                count_removed = count_removed1;
                count_added = count_added1;
                noise = noise1;
            end
        otherwise % type = 0 (non ci sono proprio template riconosciuti)
            [occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_v2(ecg_ref,occ,stop_cond,thRR,gamma,starting_delay,1,1);
    end
end





%%
function[occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_2classes(ecg_ref,...
    mark,occ,stop_cond,thRR,gamma,starting_delay,preval,idx_preval,len_longest)
    occtmp=0;
    occchk=0;
    RRm=0;
    count_removed=0;
    count_added=0;
    noise=0;
    if size(occ,1)>size(occ,2)
        occ=occ';
    end
    occtmp=occ(1,occ>0); %taglio la coda di eventuali zeri
    RRist = occtmp(2:end)-occtmp(1:end-1); %trovo gli RR istantanei solo per la mediana

    if preval == -1 % non c'è una classe prevalente, ma ci sono 2 classi
        current = mark(idx_preval);
        [start, len] =reg_analysis(occtmp);
        if len_longest >= 15
            RRm = round(median(RRist(idx_preval:idx_preval+len_longest-2)));
        else
            if len > 1
                RRm=round(median(RRist(start:start+len-2))); %mediana RR startingvalue per autoregressivo
            else 
                RRm=RRist(idx_preval);
            end
        end
    else
        current = -1;
        [y_n, idx_preval_t, len_longest_t] = preval_zeros(mark);
        if y_n == 1
            [start, len] =reg_analysis(occtmp);
            idx_preval = start;
            len_longest = len;
            RRm=round(median(RRist(start:start+len-2))); %mediana RR starting value per autoregressivo
        else
            if len_longest > 1
                RRm = round(median(RRist(idx_preval:idx_preval+len_longest-2)));
            else
                RRm = RRist(idx_preval);
            end
        end
    end
    RRmorig=RRm;
    RRist=occtmp(2:end)-occtmp(1:end-1); %trovo gli RR istantanei
    delta=abs(RRist-RRm); %identifico lo scostamento di ogni distanza RR dal valore mediano
    occchk=zeros(size(delta)); %creo vettore di check su tale scostamento(0=fuori tolleranza)
    occchk(delta<=thRR)=1; %come sopra (1=ok)
    %     start=find(occchk, 1, 'first'); %prendo il primo elemento ok e partodalì
    %     if start+2 < length(occchk)
    %         while occchk(start)+occchk(start+1)+occchk(start+2) ~= 3
    %             start = start+1;
    %             if start+2 >= length(occchk)
    %                 noise = 1;
    %                 break;
    %             end
    %         end
    %     end

    start = idx_preval;

    %prima ragioniamo in avanti (da start in poi), poi, se necessario,andremo
    %ai punti che precedono start.
    if noise == 0
        count_removed = 0;
        count_added = 0;
        i=start;
        while i<length(occtmp)-2 %fino alla fine del vettore
%             close all
%             plot(ecg_ref)
%             hold on
%             plot(occ(1:end-1),ecg_ref(occ(1:end-1)),'og');
%             plot(occtmp(1:end-1),occchk*100,'*r')
            if occchk(i) ~= 0 %se la distanza fra l'i+1 e l'i è ok per qualchemotivo...
                %non può essercene uno migliore intorno solo all'inizio, ma a
                %regime sì. Quindi...
                %provo a vedere con l'autoregressivo dove dovrebbe trovarsi laprossima occorrenza
                nextR = occtmp(i) + RRm;

                %provo a mettermi nell'intorno di tolleranza sulla posizione
                %stimata del prossimo picco per "fare pulizia" se è il caso
                j=i+1;
                best_occ = [];
                best_idx = [];
                while abs(occtmp(j)-nextR) <=  thRR
                    if isempty(best_occ) %inizializzazione
                        best_occ = abs(occtmp(j)-nextR);
                        best_idx = j;
                    else %se c'è più di un picco nella fascia di tolleranza
                        if abs(occtmp(j)-nextR) < best_occ %se il picco è piùvicino alla posizione presunta del migliore precedente
                            best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
                            best_idx = j; %conservo l'indice
                        end
                    end
                    j=j+1;
                    if j>length(occtmp)
                        break
                    end
                end
                %cancello tutti i picchi nell'intorno del migliore (uno almeno
                %lo trovo, sennò non ci sarebbe stato un occchk=1) ma cerco
                %di vedere se non sto togliendo uno un po' spostato
                %coerente con i precedenti in favore di uno non coerente
                if (~isempty(best_idx))
                    if current < 0
                        if (mark(best_idx) ~= preval) &&(any(mark(i+1:j-1)==preval))
                            % dovrei scegliere il preval più vicino
                            for aux=i+1:j-1
                                if mark(aux) ~= preval
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000; %marco per poi rimuovere
                                end
                            end
                            [best_occ, pos] = min(abs(occtmp(i+1:j-1)-nextR));
                            pos = i + pos;
                            for aux=i+1:j-1
                                if aux ~= pos
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000;
                                    count_removed = count_removed+1; %conto irimossi
                                end
                            end
                            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                            occchk(occchk==60000)=[]; %rimuovo quelli marcati
                            mark(mark==60000)=[];
                        else
                            %cancello tutti i picchi nell'intorno del migliore(uno almeno
                            %lo trovo, sennò non ci sarebbe stato un occchk=1
                            for aux=i+1:j-1
                                if aux ~= best_idx
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000; %marco per poi rimuovere
                                    count_removed = count_removed+1; %conto irimossi
                                end
                            end
                            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                            occchk(occchk==60000)=[]; %rimuovo quelli marcati
                            mark(mark==60000)=[];
                        end
                    else
                        if (mark(best_idx) ~= current) &&(any(mark(i+1:j-1)==current))
                            % dovrei scegliere il preval più vicino
                            for aux=i+1:j-1
                                if mark(aux) ~= current
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000; %marco per poi rimuovere
                                end
                            end
                            [best_occ, pos] = min(abs(occtmp(i+1:j-1)-nextR));
                            pos = i + pos;
                            for aux=i+1:j-1
                                if aux ~= pos
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000;
                                    count_removed = count_removed+1; %conto i rimossi
                                end
                            end
                            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                            occchk(occchk==60000)=[]; %rimuovo quelli marcati
                            mark(mark==60000)=[];
                        else
                            %cancello tutti i picchi nell'intorno del migliore(uno almeno
                            %lo trovo, sennò non ci sarebbe stato un occchk=1
                            current = mark(best_idx);
                            for aux=i+1:j-1
                                if aux ~= best_idx
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000; %marco per poi rimuovere
                                    count_removed = count_removed+1; %conto irimossi
                                end
                            end
                            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                            occchk(occchk==60000)=[]; %rimuovo quelli marcati
                            mark(mark==60000)=[];
                        end
                    end
                end
                %sistemo i vettori occchk e occtmp
                if(i>length(occtmp)-2)
                    break;
                end
                if j-1>i+1 %se ne ho trovato solo uno lo lascio indicato com'è
                    if abs(occtmp(i+2)-occtmp(i+1)-RRm) <= thRR
                        occchk(i+1)=2;
                    else
                        occchk(i+1)=0;
                    end
                end

                %aggiorno l'autoregressivo
                RRm=round(RRm*(1-gamma)+(occtmp(i+1)-occtmp(i))*gamma);


            else %se la distanza fra l'i+1 e l'i è fuori tolleranza...
                %provo a vedere con l'autoregressivo dove dovrebbe trovarsi la prossima occorrenza
                nextR = occtmp(i) + RRm;

                %quindi tolgo tutto quello fra il primo picco "buono" e ilmargine
                %inferiore di tolleranza rispetto a dove dovrebbe essere il
                %prossimo picco
                j=i+1; %mi posiziono sulla prossima occorrenza%================================================================))))))))))))
                while occtmp(j) < nextR - thRR %fino a che la distanza fra dueoccorrenze non supera la zona "buona" intorno al prossimo picco presunto
                    occtmp(j)=60000; %marco per poi rimuovere
                    occchk(j)=60000; %marco per poi rimuovere
                    mark(j)=60000; %marco per poi rimuovere
                    count_removed = count_removed+1; %conto i rimossi
                    j=j+1; %provo ad andre avanti di uno
                    if j>length(occchk)
                        break
                    end
                end
                occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                occchk(occchk==60000)=[]; %rimuovo quelli marcati
                mark(mark==60000)=[]; %rimuovo quelli marcati
                if i > length(occtmp)-2
                    break
                end
                %                 if abs(occtmp(i+1)-occtmp(i)-RRm) <= thRR
                % %                 if abs(occtmp(i+1)-occtmp(i)-nextR) < thRR
                %                     occchk(i)=2;
                %                 else
                %                     occchk(i)=0;
                %                 end
                if abs(occtmp(i+1)-occtmp(i)-RRm) <= thRR
                    %                 if abs(occtmp(i+1)-occtmp(i)-nextR) < thRR
                    occchk(i)=2;
                else
                    occchk(i)=0;
                end

                nextR = occtmp(i)+RRm;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %provo a mettermi nell'intorno di tolleranza sulla posizione
                %stimata del prossimo picco per "fare pulizia" se è il caso
                j=i+1;
                best_occ = [];
                best_idx = [];
                while abs(occtmp(j)-nextR) <=  thRR
                    if isempty(best_occ) %inizializzazione
                        best_occ = abs(occtmp(j)-nextR);
                        best_idx = j;
                    else %se c'è più di un picco nella fascia di tolleranza
                        if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
                            best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
                            best_idx = j; %conservo l'indice
                        end
                    end
                    j=j+1;
                    if j>length(occtmp)
                        break
                    end
                end
                if isempty(best_occ) %vuol dire che non c'erano elementi nella zona intorno alla posizione presunta
                    %devo aggiungere un picco nella posizione presunta
                    occtmp=[occtmp(1:i),nextR,occtmp(i+1:end)];
                    occchk=[occchk(1:i),3,occchk(i+1:end)]; %col codice 3 identifico un picco che è stato aggiunto
                    mark=[mark(1:i),-1,mark(i+1:end)]; %col codice -1 identificoun picco che è stato aggiunto quindi non ha passato template matching
                    count_added = count_added+1;
                    if abs(occtmp(i+2)-occtmp(i+1)-RRm) <= thRR
                        occchk(i+1)=3;
                    else
                        occchk(i+1)=0;
                    end
                    occchk(i)=3;
                else %ho trovato un picco molto prossimo alla posizione presunta
                    %cancello tutti gli altri trovati
                    %cancello tutti i picchi nell'intorno del migliore (uno almeno
                    %lo trovo, sennò non ci sarebbe stato un occchk=1) ma cerco
                    %di vedere se non sto togliendo uno un po' spostato
                    %coerente con i precedenti in favore di uno non coerente
                    if ~isempty(best_idx)
                        if current < 0
                            if (mark(best_idx) ~= preval) &&(any(mark(i+1:j-1)==preval))
                                % dovrei scegliere il preval più vicino
                                for aux=i+1:j-1
                                    if mark(aux) ~= preval
                                        occtmp(aux)=60000; %marco per poirimuovere
                                        occchk(aux)=60000; %marco per poirimuovere
                                        mark(aux)=60000; %marco per poirimuovere
                                    end
                                end
                                [best_occ, pos] =min(abs(occtmp(i+1:j-1)-nextR));
                                pos = i + pos;
                                for aux=i+1:j-1
                                    if aux ~= pos
                                        occtmp(aux)=60000; %marco per poirimuovere
                                        occchk(aux)=60000; %marco per poirimuovere
                                        mark(aux)=60000;
                                        count_removed = count_removed+1; %contoi rimossi
                                    end
                                end
                                occtmp(occtmp==60000)=[]; %rimuovo quellimarcati
                                occchk(occchk==60000)=[]; %rimuovo quellimarcati
                                mark(mark==60000)=[];
                            else
                                %cancello tutti i picchi nell'intorno del migliore (uno almeno
                                %lo trovo, sennò non ci sarebbe stato un occchk=1
                                for aux=i+1:j-1
                                    if aux ~= best_idx
                                        occtmp(aux)=60000; %marco per poi rimuovere
                                        occchk(aux)=60000; %marco per poi rimuovere
                                        mark(aux)=60000; %marco per poi rimuovere
                                        count_removed = count_removed+1; %conto i rimossi
                                    end
                                end
                                occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                                occchk(occchk==60000)=[]; %rimuovo quelli marcati
                                mark(mark==60000)=[];
                            end
                        else
                            if (mark(best_idx) ~= current) && (any(mark(i+1:j-1)==current))
                                % dovrei scegliere il preval più vicino
                                for aux=i+1:j-1
                                    if mark(aux) ~= current
                                        occtmp(aux)=60000; %marco per poi rimuovere
                                        occchk(aux)=60000; %marco per poi rimuovere
                                        mark(aux)=60000; %marco per poi rimuovere
                                    end
                                end
                                [best_occ, pos] =min(abs(occtmp(i+1:j-1)-nextR));
                                pos = i + pos;
                                for aux=i+1:j-1
                                    if aux ~= pos
                                        occtmp(aux)=60000; %marco per poi rimuovere
                                        occchk(aux)=60000; %marco per poi rimuovere
                                        mark(aux)=60000;
                                        count_removed = count_removed+1; %conto i rimossi
                                    end
                                end
                                occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                                occchk(occchk==60000)=[]; %rimuovo quelli marcati
                                mark(mark==60000)=[];
                            else
                                current = mark(best_idx);
                                %cancello tutti i picchi nell'intorno del migliore (uno almeno
                                %lo trovo, sennò non ci sarebbe stato un occchk=1
                                for aux=i+1:j-1
                                    if aux ~= best_idx
                                        occtmp(aux)=60000; %marco per poi rimuovere
                                        occchk(aux)=60000; %marco per poi rimuovere
                                        mark(aux)=60000; %marco per poi rimuovere
                                        count_removed = count_removed+1; %conto i rimossi
                                    end
                                end
                                occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                                occchk(occchk==60000)=[]; %rimuovo quelli marcati
                                mark(mark==60000)=[];
                            end
                        end
                    end
                    %sistemo i vettori occchk e occtmp
                    if(i>length(occtmp)-2)
                        break;
                    end
                    if j-1>i+1
                        if abs(occtmp(i+1)-occtmp(i)-RRm) <= thRR
                            occchk(i)=2;
                        else
                            occchk(i)=0;
                        end
                    end

                end
            end
            i=i+1;
        end

        occtmp(i+1:end)=[];
        occchk(i:end)=[];
        mark(i:end)=[];

        %coda se non sono state trovate occorrenze fino alla fine del segnale
        while stop_cond - occtmp(length(occtmp)) > RRm
            nextR = occtmp(length(occtmp)) + RRm;
            occtmp=[occtmp,nextR];
            occchk=[occchk,3]; %col codice 3 identifico un picco che è stato aggiunto
            mark=[mark,-1];%col codice -1 identifico un picco che è stato aggiunto
            count_added = count_added+1;
        end



        % adesso torniamo indietro
        i=start;
        RRm = RRmorig;
        while i>3 %fino alla testa del vettore (corretto! se metto 1 da errore,indici negativi)
%             close all
%             plot(ecg_ref)
%             hold on
%             plot(occ(1:end-1),ecg_ref(occ(1:end-1)),'og');
%             plot(occtmp(1:end-1),occchk*100,'*r')
            if occchk(i-1) ~= 0 %se la distanza fra l'i-1 e l'i è ok per qualche motivo...
                %non può essercene uno migliore intorno solo all'inizio, ma a
                %regime sì. Quindi...
                %provo a vedere con l'autoregressivo dove dovrebbe trovarsi la precedente occorrenza
                nextR = occtmp(i) - RRm;
                if(nextR<starting_delay)
                    break
                end

                %provo a mettermi nell'intorno di tolleranza sulla posizione
                %stimata del prossimo picco per "fare pulizia" se è il caso
                j=i-1;
                best_occ = [];
                best_idx = [];
                while abs(occtmp(j)-nextR) <=  thRR
                    if isempty(best_occ) %inizializzazione
                        best_occ = abs(occtmp(j)-nextR);
                        best_idx = j;
                    else %se c'è più di un picco nella fascia di tolleranza
                        if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
                            best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
                            best_idx = j; %conservo l'indice
                        end
                    end
                    j=j-1;
                    if j<1
                        break
                    end
                end
                %cancello tutti i picchi nell'intorno del migliore (uno almeno
                %lo trovo, sennò non ci sarebbe stato un occchk=1) ma cerco
                %di vedere se non sto togliendo uno un po' spostato
                %coerente con i precedenti in favore di uno non coerente
                if ~isempty(best_idx)
                    if preval < 0
                        if (mark(best_idx) ~= preval) &&(any(mark(j+1:i-1)==preval))
                            % dovrei scegliere il preval più vicino
                            for aux=i-1:-1:j+1
                                if mark(aux) ~= preval
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000; %marco per poi rimuovere
                                end
                            end
                            [best_occ, pos] = min(abs(occtmp(j+1:i-1)-nextR));
                            pos = i + pos;
                            for aux = i-1:-1:j+1
                                if aux ~= pos
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000;
                                    count_removed = count_removed+1; %conto irimossi
                                end
                            end

                            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                            occchk(occchk==60000)=[]; %rimuovo quelli marcati
                            mark(mark==60000)=[];

                        else
                            %cancello tutti i picchi nell'intorno del migliore (uno almeno
                            %lo trovo, sennò non ci sarebbe stato un occchk=1
                            for aux=i-1:-1:j+1
                                if aux ~= best_idx
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000; %marco per poi rimuovere
                                    count_removed = count_removed+1; %conto irimossi
                                end
                            end
                            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                            occchk(occchk==60000)=[]; %rimuovo quelli marcati
                            mark(mark==60000)=[];
                        end
                    else
                        if (mark(best_idx) ~= current) &&(any(mark(j+1:i-1)==current))
                            % dovrei scegliere il preval più vicino
                            for aux=i-1:-1:j+1
                                if mark(aux) ~= current
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000; %marco per poi rimuovere
                                end
                            end
                            [best_occ, pos] = min(abs(occtmp(j+1:i-1)-nextR));
                            pos = i + pos;
                            for aux = i-1:-1:j+1
                                if aux ~= pos
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000;
                                    count_removed = count_removed+1; %conto i rimossi
                                end
                            end

                            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                            occchk(occchk==60000)=[]; %rimuovo quelli marcati
                            mark(mark==60000)=[];

                        else
                            %cancello tutti i picchi nell'intorno del migliore (uno almeno
                            %lo trovo, sennò non ci sarebbe stato un occchk=1
                            current = mark(best_idx);
                            for aux=i-1:-1:j+1
                                if aux ~= best_idx
                                    occtmp(aux)=60000; %marco per poi rimuovere
                                    occchk(aux)=60000; %marco per poi rimuovere
                                    mark(aux)=60000; %marco per poi rimuovere
                                    count_removed = count_removed+1; %conto i rimossi
                                end
                            end
                            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                            occchk(occchk==60000)=[]; %rimuovo quelli marcati
                            mark(mark==60000)=[];
                        end
                    end
                end
                if i < 3
                    break
                end
                %sistemo i vettori occchk e occtmp
                if j+1<i-1
                    if abs(occtmp(i-1)-occtmp(i-2)-RRm) <= thRR
                        occchk(i-2)=2;
                    else
                        occchk(i-2)=0;
                    end
                end

                %aggiorno l'autoregressivo
                RRm=round(RRm*(1-gamma)+(occtmp(i)-occtmp(i-1))*gamma);
            else %se la distanza fra l'i e l'i-1 è fuori tolleranza...
                %provo a vedere con l'autoregressivo dove dovrebbe trovarsi
                %la prossima occorrenza
                nextR = occtmp(i) - RRm;

                %quindi tolgo tutto quello fra il primo picco "buono" e il margine
                %inferiore di tolleranza rispetto a dove dovrebbe essere il
                %prossimo picco
                j=i-1; %mi posiziono sulla prossima occorrenza
                while occtmp(j) >= nextR + thRR %fino a che la distanza fra due occorrenze non supera la zona "buona" intorno al prossimo picco presunto
                    occtmp(j)=60000; %marco per poi rimuovere
                    occchk(j)=60000; %marco per poi rimuovere
                    mark(j)=60000; %marco per poi rimuovere
                    count_removed = count_removed+1; %conto i rimossi
                    j=j-1; %provo ad andre avanti di uno
                    if j<1
                        break
                    end
                end
                if j>1
                    if abs(occtmp(i)-occtmp(j)-RRm) <=thRR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        occchk(j)=2;
                    else
                        occchk(j)=0;
                    end
                end
                occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                occchk(occchk==60000)=[]; %rimuovo quelli marcati
                mark(mark==60000)=[]; %rimuovo quelli marcati

                nextR = occtmp(i) - RRm;

                %provo a mettermi nell'intorno di tolleranza sulla posizione
                %stimata del prossimo picco per "fare pulizia" se è il caso
                j=i-1;
                best_occ = [];
                best_idx = [];
                while abs(occtmp(j)-nextR) <=  thRR
                    if isempty(best_occ) %inizializzazione
                        best_occ = abs(occtmp(j)-nextR);
                        best_idx = j;
                    else %se c'è più di un picco nella fascia di tolleranza
                        if abs(occtmp(j)-nextR) < best_occ %se il picco è piùvicino alla posizione presunta del migliore precedente
                            best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
                            best_idx = j; %conservo l'indice
                        end
                    end
                    j=j-1;
                    if j<1
                        break
                    end
                end
                if isempty(best_occ) %vuol dire che non c'erano elementi nella zona intorno alla posizione presunta
                    %devo aggiungere un picco nella posizione presunta
                    occtmp=[occtmp(1:i-1),nextR,occtmp(i:end)];
                    occchk=[occchk(1:i-1),3,occchk(i:end)]; %col codice 3 identifico un picco che è stato aggiunto%%%%%%%%%%%%%%
                    mark=[mark(1:i-1),-1,mark(1:end)]; %col codice -1 identifico un picco che è stato aggiunto quindi non ha passato template matching
                    count_added = count_added+1;
                    if abs(occtmp(i)-occtmp(i-1)-RRm) <= thRR
                        occchk(i-1)=3;
                    else
                        occchk(i-1)=0;
                    end
                    i=i+1;
                else %ho trovato un picco molto prossimo alla posizione presunta
                    %cancello tutti gli altri trovati
                    %cancello tutti i picchi nell'intorno del migliore (uno almeno
                    %lo trovo, sennò non ci sarebbe stato un occchk=1) ma cerco
                    %di vedere se non sto togliendo uno un po' spostato
                    %coerente con i precedenti in favore di uno non
                    %coerente
                    if ~isempty(best_idx)
                        if preval < 0
                            if (mark(best_idx) ~= preval) &&(any(mark(j+1:i-1)==preval))
                                % dovrei scegliere il preval più vicino
                                for aux=i-1:-1:j+1
                                    if mark(aux) ~= preval
                                        occtmp(aux)=60000; %marco per poi rimuovere
                                        occchk(aux)=60000; %marco per poi rimuovere
                                        mark(aux)=60000; %marco per poi rimuovere
                                    end
                                end
                                [best_occ, pos] =min(abs(occtmp(j+1:i-1)-nextR));
                                pos = i + pos;

                                for aux=i-1:-1:j+1
                                    if aux ~= best_idx
                                        occtmp(aux)=60000; %marco per poi rimuovere
                                        occchk(aux)=60000; %marco per poi rimuovere
                                        mark(aux)=60000;
                                        count_removed = count_removed+1; %conto i rimossi
                                    end
                                end
                                occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                                occchk(occchk==60000)=[]; %rimuovo quelli marcati
                                mark(mark==60000)=[];
                            else
                                %cancello tutti i picchi nell'intorno del migliore (uno almeno
                                %lo trovo, sennò non ci sarebbe stato un occchk=1
                                for aux=i-1:-1:j+1
                                    if aux ~= best_idx
                                        occtmp(aux)=60000; %marco per poi rimuovere
                                        occchk(aux)=60000; %marco per poi rimuovere
                                        mark(aux)=60000; %marco per poi rimuovere
                                        count_removed = count_removed+1; %conto i rimossi
                                    end
                                end
                                occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                                occchk(occchk==60000)=[]; %rimuovo quelli marcati
                                mark(mark==60000)=[];
                            end
                        else
                            if (mark(best_idx) ~= current) &&(any(mark(j+1:i-1)==current))
                                % dovrei scegliere il preval più vicino
                                for aux=i-1:-1:j+1
                                    if mark(aux) ~= current
                                        occtmp(aux)=60000; %marco per poi rimuovere
                                        occchk(aux)=60000; %marco per poi rimuovere
                                        mark(aux)=60000; %marco per poi rimuovere
                                    end 
                                end
                                [best_occ, pos] =min(abs(occtmp(j+1:i-1)-nextR));
                                pos = i + pos;

                                for aux=i-1:-1:j+1
                                    if aux ~= best_idx
                                        occtmp(aux)=60000; %marco per poi rimuovere
                                        occchk(aux)=60000; %marco per poi rimuovere
                                        mark(aux)=60000;
                                        count_removed = count_removed+1; %conto i rimossi
                                    end
                                end
                                occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                                occchk(occchk==60000)=[]; %rimuovo quelli marcati
                                mark(mark==60000)=[];
                            else
                                current = mark(best_idx);
                                %cancello tutti i picchi nell'intorno del migliore (uno almeno
                                %lo trovo, sennò non ci sarebbe stato un occchk=1
                                for aux=i-1:-1:j+1
                                    if aux ~= best_idx
                                        occtmp(aux)=60000; %marco per poi rimuovere
                                        occchk(aux)=60000; %marco per poi rimuovere
                                        mark(aux)=60000; %marco per poi rimuovere
                                        count_removed = count_removed+1; %conto i rimossi
                                    end
                                end
                                occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                                occchk(occchk==60000)=[]; %rimuovo quelli marcati
                                mark(mark==60000)=[];
                            end
                        end
                    end
                    %sistemo i vettori occchk e occtmp
                    if j+1<i-1
                        i=j+2;
                        if abs(occtmp(i)-occtmp(i-1)-RRm) <= thRR
                            occchk(i-1)=2;
                        else
                            occchk(i-1)=0;
                        end
                    end
                end
            end
            i=i-1;
        end
        %coda se non inizia subito il vettore delle occorrenze
        while occtmp(1) > RRm
            nextR = occtmp(1) - RRm;
            if(nextR<starting_delay)
                break
            end
            occtmp=[nextR,occtmp];
            occchk=[3,occchk]; %col codice 3 identifico un picco che è statoaggiunto
            mark=[mark,-1];%col codice -1 identifico un picco che è stato aggiunto
            count_added = count_added+1;
        end

    end
%             close all
%             plot(ecg_ref)
%             hold on
%             plot(occ(1:end-1),ecg_ref(occ(1:end-1)),'og');
%             plot(occtmp(1:end-1),occchk*100,'*r')
    
end



% [preval] = prevalence(a)
%
% a partire dallo spezzone più lungo di "mark" che non contiene zeri,
% restituisce l'etichetta del template più rappresentato

%%
function [preval, idx_preval, len_longest, type] = prevalence(a)

    % cerca la prima occorrenza non nulla, che etichetta come seed
    ix = find(a, 1, 'first');
    if isempty(ix)
        preval=0;
        idx_preval=1;
        len_longest = 1;
        type = 0; % vuol dire che non ci sono template che sono stati trovati
        return;
    end

    % in caso contrario, la prima occorrenza ci dà il seme
    seed = a(ix);
    idx_preval_seed = ix;

    if ~any((a~=seed)&(a~=0)) % se oltre a zeri c'è rappresentata solo unaclasse
        preval = seed;
        [idx_preval, len_longest] = longest_seq(a, seed);
        type = 1; % vuol dire che c'è solo una classe
    else % se invece ci sono due classi le etichetta come seed e other
        if seed == 1
            other = 2;
        else
            other = 1;
        end
        % a questo punto controllo se esiste una classe prevalente oppure se
        % siamo in una condizione nella quale c'è sostanziale alternanza
        seed_vec_idx = find(a==seed);
        other_vec_idx = find(a==other);

        if min(length(seed_vec_idx),length(other_vec_idx)) <0.3*max(length(seed_vec_idx),length(other_vec_idx))
            % vuol dire che c'è una classe molto più rappresentata dell'altra
            type = 21; % lo indico così
            if length(seed_vec_idx) > length(other_vec_idx)
                [idx_preval, len_longest] = longest_seq(a, seed);
                preval = seed;
            else
                [idx_preval, len_longest] = longest_seq(a, other);
                preval = other;
            end
        else % cioè se invece non esiste questa relazione, cerco se esiste alternanza
            temp = a(a~=0); % tolgo gli zeri
            count_alternances = 0;
            for i = 1:length(temp)-1
                if temp(i+1)~=temp(i)
                    count_alternances = count_alternances + 1;
                end
            end
            if abs(count_alternances-length(temp)) < 0.7*length(temp)
                % siamo in condizioni di alternanza quindi è meglio analizzare
                % le due sequenze separatamente
                type = 22; % indico l'alternanza
                if length(seed_vec_idx) > length(other_vec_idx)
                    preval = seed;
                    idx_preval = 0;
                    len_longest = 0;
                else
                    preval = other;
                    idx_preval = 0;
                    len_longest = 0;
                end
            else
                % non siamo in condizioni di alternanza. Non possiamo sfruttare
                % ulteriori informazioni in analisi di periodicità
                type = 20;
                no_preval = a;
                no_preval(no_preval~=0)=1;
                [idx_preval, len_longest] = longest_seq(no_preval, 1);
                preval=-1;
            end
        end
    end
end

%%
function [y_n, idx_preval, len_longest] = preval_zeros(a)

    % cerca la prima occorrenza non nulla, che etichetta come seed
    if ~any(a~=0)
        y_n=1;
        idx_preval=1;
        len_longest = length(a);
        return;
    end

    % in caso contrario, la prima occorrenza ci dà il seme
    a(a~=0)=1;
    zeros_vec_idx = find(a==0);
    other_vec_idx = find(a==1);

    if length(other_vec_idx) < 0.4*length(zeros_vec_idx)
        y_n=1;
        [idx_preval, len_longest] = longest_seq(a, 0);
    else
        y_n=0;
        [idx_preval, len_longest] = longest_seq(a, 1);
    end
end

%%
% function [] = match(QRSs, MQRS, thmm)
%
% thmm è la tolleranza in matching

function [mark] = match(QRSs, MQRS, thmm)  %thmm 15---> ora è 15*4

    if nargin <3
        thmm=15;
    end
    mark=zeros(size(QRSs));
%     tt=MQRS+25;
    tt=MQRS+25*4;
%     tt(tt==25)=1000000;
    tt(tt==25*4)=1000000;

    i=1;
    while i<=length(mark)
        m = min(min(tt));
        if abs(QRSs(i)-m) < thmm
            [r,c]=find(tt==m);
            mark(i) = r(1,1);
            tt(r,c)=1000000;
            i=i+1;
        elseif m <= QRSs(i)
            [r,c]=find(tt==m);
            tt(r,c)=1000000;
        else
            i=i+1;
        end
    end
end

%%
%%
function [start_longest, length_longest] = longest_seq(a, key)

    a_max = 0;
    i_start_max = 0;
    a_count = 0;
    i_start = 0;

    for i=1:length(a)
        if a(i) == key
            a_count = a_count+1;
            if a_count == 1
                i_start = i;
            end
        else
            if a_max < a_count
                a_max = a_count;
                i_start_max = i_start;
            end
            a_count = 0;
            i_start = 0;
        end
    end
    
if a_count > a_max
    i_start_max = i_start;
    a_max = a_count;
end

    start_longest = i_start_max;
    length_longest = a_max;
end

%%
function [start, len] =reg_analysis(a)

    d = a(2:end)-a(1:end-1);
    c = d(2:end)-d(1:end-1);

    c(c < 10) = 1;
    c(c >= 10) = 0;

    [start, len] = longest_seq(c, 1);
    len = len+2;
end

%%
function M = dist_matrix(d)
    l = length(d);
    M=zeros(l,l);
    for i = 1:l
        for j = 1:l
            M(i,j)=abs(d(i)-d(j));
        end
    end
end

function[occtmp,occchk,RRm,count_removed,count_added,noise]=periodicity_correction_v2(ecg_ref,...
occ,stop_cond,thRR,gamma,starting_delay,idx_preval,len_longest)
    occtmp=0;
    occchk=0;
    RRm=0;
    count_removed=0;
    count_added=0;
    noise=0;
    if size(occ,1)>size(occ,2)
    occ=occ';
    end
    occtmp=occ(1,occ>0); %taglio la coda di eventuali zeri
    RRist=occtmp(2:end)-occtmp(1:end-1); %trovo gli RR istantanei
    RRm=round(median(RRist)); %mediana RR starting value per autoregressivo
    %RRm=100;%%%%%%%%%%%%%%%%%%%%%%%%%Solo per test.
    delta=abs(RRist-RRm); %identifico lo scostamento di ogni distanza RR dalvalore mediano
    occchk=zeros(size(delta)); %creo vettore di check su tale    scostamento(0=fuori tolleranza)
    occchk(delta<=thRR)=1; %come sopra (1=ok)
    %     start=find(occchk, 1, 'first'); %prendo il primo elemento ok e parto da lì
    %     if start+2 < length(occchk)
    %         while occchk(start)+occchk(start+1)+occchk(start+2) ~= 3
    %             start = start+1;
    %             if start+2 >= length(occchk)
    %                 noise = 1;
    %                 break;
    %             end
    %         end
    %     end

    start = idx_preval;

    %prima ragioniamo in avanti (da start in poi), poi, se necessario,andremo
    %ai punti che precedono start.
    if noise == 0
    count_removed = 0;
    count_added = 0;
    i=start;
    if len_longest >= 5
        RRm = occtmp(i+1)-occtmp(i);
    end
    while i<length(occtmp)-2 %fino alla fine del vettore
% 
%         close all
%         plot(ecg_ref)
%         hold on
%         plot(occ(1:end-1),ecg_ref(occ(1:end-1)),'og');
%         plot(occtmp(1:end-1),occchk*100,'*r')
        if occchk(i) ~= 0 %se la distanza fra l'i+1 e l'i è ok per qualchemotivo...
            %non può essercene uno migliore intorno solo all'inizio, ma a
            %regime sì. Quindi...
            %provo a vedere con l'autoregressivo dove dovrebbe trovarsi la prossima occorrenza
            nextR = occtmp(i) + RRm;

            %provo a mettermi nell'intorno di tolleranza sulla posizione
            %stimata del prossimo picco per "fare pulizia" se è il caso
            j=i+1;
            best_occ = [];
            best_idx = [];
            while abs(occtmp(j)-nextR) <=  thRR
                if isempty(best_occ) %inizializzazione
                    best_occ = abs(occtmp(j)-nextR);
                    best_idx = j;
                else %se c'è più di un picco nella fascia di tolleranza
                    if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
                        best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
                        best_idx = j; %conservo l'indice
                    end
                end
                j=j+1;
                if j>length(occtmp)
                    break
                end
            end
            %cancello tutti i picchi nell'intorno del migliore (uno almeno
            %lo trovo, sennò non ci sarebbe stato un occchk=1
            for aux=i+1:j-1
                if aux ~= best_idx
                    occtmp(aux)=60000; %marco per poi rimuovere
                    occchk(aux)=60000; %marco per poi rimuovere
                    count_removed = count_removed+1; %conto i rimossi
                end
            end
            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
            occchk(occchk==60000)=[]; %rimuovo quelli marcati
            %sistemo i vettori occchk e occtmp
            if(i>length(occtmp)-2)
                break;
            end
            if j-1>=i+2 %se ne ho trovato solo uno lo lascio indicato com'è
                if abs(occtmp(i+2)-occtmp(i+1)-RRm) <= thRR
                    occchk(i+1)=1;
                else
                    occchk(i+1)=0;
                end
            end

            %aggiorno l'autoregressivo
            RRm=round(RRm*(1-gamma)+(occtmp(i+1)-occtmp(i))*gamma);
        else %se la distanza fra l'i+1 e l'i è fuori tolleranza...
            %provo a vedere con l'autoregressivo dove dovrebbe trovarsi la prossima occorrenza
            nextR = occtmp(i) + RRm;

            %quindi tolgo tutto quello fra il primo picco "buono" e il margine
            %inferiore di tolleranza rispetto a dove dovrebbe essere il
            %prossimo picco
            j=i+1; %mi posiziono sulla prossima occorrenza
            while occtmp(j) < nextR - thRR %fino a che la distanza fra due occorrenze non supera la zona "buona" intorno al prossimo picco presunto
                occtmp(j)=60000; %marco per poi rimuovere
                occchk(j)=60000; %marco per poi rimuovere
                count_removed = count_removed+1; %conto i rimossi
                j=j+1; %provo ad andre avanti di uno
                if j>length(occchk)
                    break
                end
            end
            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
            occchk(occchk==60000)=[]; %rimuovo quelli marcati
            if abs(occtmp(i+1)-occtmp(i)-RRm) <= thRR
                %                 if abs(occtmp(i+1)-occtmp(i)-nextR) < thRR
                occchk(i)=2;
            else
                occchk(i)=0;
            end

            nextR = occtmp(i) +RRm;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %provo a mettermi nell'intorno di tolleranza sulla posizione
            %stimata del prossimo picco per "fare pulizia" se è il caso
            j=i+1;
            best_occ = [];
            best_idx = [];
            while abs(occtmp(j)-nextR) <=  thRR
                if isempty(best_occ) %inizializzazione
                    best_occ = abs(occtmp(j)-nextR);
                    best_idx = j;
                else %se c'è più di un picco nella fascia di tolleranza
                    if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
                        best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
                        best_idx = j; %conservo l'indice
                    end
                end
                j=j+1;
                if j>length(occtmp)
                    break
                end
            end
            if isempty(best_occ) %vuol dire che non c'erano elementi nella zona intorno alla posizione presunta
                %devo aggiungere un picco nella posizione presunta
                occtmp=[occtmp(1:i),nextR,occtmp(i+1:end)];
                occchk=[occchk(1:i),3,occchk(i+1:end)]; %col codice 3 identifico un picco che è stato aggiunto
                count_added = count_added+1;
                if abs(occtmp(i+2)-occtmp(i+1)-RRm) <= thRR
                    occchk(i+1)=3;
                else
                    occchk(i+1)=0;
                end
                occchk(i)=3;
            else %ho trovato un picco molto prossimo alla posizione presunta
                %cancello tutti gli altri trovati
                for aux=i+1:j-1
                    if aux ~= best_idx
                        occtmp(aux)=60000; %marco per poi rimuovere
                        occchk(aux)=60000; %marco per poi rimuovere
                        count_removed = count_removed+1; %conto i rimossi
                    end
                end
                occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                occchk(occchk==60000)=[]; %rimuovo quelli marcati
                if(i>length(occtmp)-2)
                    break;
                end
                %sistemo i vettori occchk e occtmp
                if j-1>=i+2 %se ne ho trovato solo uno lo lascio indicato com'è
                    if abs(occtmp(i+2)-occtmp(i+1)-RRm) <= thRR
                        occchk(i+1)=1;
                    else
                        occchk(i+1)=0;
                    end
                end
            end
        end
        i=i+1;
    end
    %coda se non sono state trovate occorrenze fino alla fine del segnale
    while stop_cond - occtmp(length(occtmp)) > RRm
        nextR = occtmp(length(occtmp)) + RRm;
        occtmp=[occtmp,nextR];
        occchk=[occchk,3]; %col codice 3 identifico un picco che è stato aggiunto
        count_added = count_added+1;
    end


    % adesso torniamo indietro
    i=start;
    if len_longest >= 5
        RRm = occtmp(i+1)-occtmp(i);
    end
    %RRm=round(median(RRist)); %mediana RR starting value per autoregressivo
    %RRm=100;
    while i>3 %fino alla testa del vettore (corretto! se metto 1 da errore, indici negativi)
% 
%         close all
%         plot(ecg_ref)
%         hold on
%         plot(occ(1:end-1),ecg_ref(occ(1:end-1)),'og');
%         plot(occtmp(1:end-1),occchk*100,'*r')
        if occchk(i-1) ~= 0 %se la distanza fra l'i-1 e l'i è ok per qualche motivo...
            %non può essercene uno migliore intorno solo all'inizio, ma a
            %regime sì. Quindi...
            %provo a vedere con l'autoregressivo dove dovrebbe trovarsi la precedente occorrenza
            nextR = occtmp(i) - RRm;
            if(nextR<starting_delay)
                break
            end

            %provo a mettermi nell'intorno di tolleranza sulla posizione
            %stimata del prossimo picco per "fare pulizia" se è il caso
            j=i-1;
            best_occ = [];
            best_idx = [];
            while abs(occtmp(j)-nextR) <=  thRR
                if isempty(best_occ) %inizializzazione
                    best_occ = abs(occtmp(j)-nextR);
                    best_idx = j;
                else %se c'è più di un picco nella fascia di tolleranza
                    if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
                        best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
                        best_idx = j; %conservo l'indice
                    end
                end
                j=j-1;
                if j<1
                    break
                end
            end
            %cancello tutti i picchi nell'intorno del migliore (uno almeno
            %lo trovo, sennò non ci sarebbe stato un occchk=1
            for aux=i-1:-1:j+1
                if aux ~= best_idx
                    occtmp(aux)=60000; %marco per poi rimuovere
                    occchk(aux)=60000; %marco per poi rimuovere
                    count_removed = count_removed+1; %conto i rimossi
                end
            end
            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
            occchk(occchk==60000)=[]; %rimuovo quelli marcati
            %sistemo i vettori occchk e occtmp
            if j+1<i-1 %se ne ho trovato solo uno lo lascio indicato com'è
                if abs(occtmp(i-1)-occtmp(i-2)-RRm) <= thRR
                    occchk(i-2)=1;
                else
                    occchk(i-2)=0;
                end
            end

            %aggiorno l'autoregressivo
            RRm=round(RRm*(1-gamma)+(occtmp(i)-occtmp(i-1))*gamma);
        else %se la distanza fra l'i e l'i-1 è fuori tolleranza...
            %provo a vedere con l'autoregressivo dove dovrebbe trovarsi
            %la prossima occorrenza
            nextR = occtmp(i) - RRm;

            %quindi tolgo tutto quello fra il primo picco "buono" e il margine
            %inferiore di tolleranza rispetto a dove dovrebbe essere il
            %prossimo picco
            j=i-1; %mi posiziono sulla prossima occorrenza
            while occtmp(j) >= nextR + thRR %fino a che la distanza fra due occorrenze non supera la zona "buona" intorno al prossimo picco presunto
                occtmp(j)=60000; %marco per poi rimuovere
                occchk(j)=60000; %marco per poi rimuovere
                count_removed = count_removed+1; %conto i rimossi
                j=j-1; %provo ad andre avanti di uno
                if j<1
                    break
                end
            end
            if abs(occtmp(i)-occtmp(j)-RRm) <= thRR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                occchk(j)=2;
            else
                occchk(j)=0;
            end
            occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
            occchk(occchk==60000)=[]; %rimuovo quelli marcati

            nextR = occtmp(i) - RRm;

            %provo a mettermi nell'intorno di tolleranza sulla posizione
            %stimata del prossimo picco per "fare pulizia" se è il caso
            j=i-1;
            best_occ = [];
            best_idx = [];
            while abs(occtmp(j)-nextR) <=  thRR
                if isempty(best_occ) %inizializzazione
                    best_occ = abs(occtmp(j)-nextR);
                    best_idx = j;
                else %se c'è più di un picco nella fascia di tolleranza
                    if abs(occtmp(j)-nextR) < best_occ %se il picco è più vicino alla posizione presunta del migliore precedente
                        best_occ = abs(occtmp(j)-nextR); %questo diventa il migliore
                        best_idx = j; %conservo l'indice
                    end
                end
                j=j-1;
                if j<1
                    break
                end
            end
            if isempty(best_occ) %vuol dire che non c'erano elementi nella zona intorno alla posizione presunta
                %devo aggiungere un picco nella posizione presunta
                occtmp=[occtmp(1:i-1),nextR,occtmp(i:end)];
                occchk=[occchk(1:i-1),3,occchk(i:end)]; %col codice 3 identifico un picco che è stato aggiunto%%%%%%%%%%%%%%
                count_added = count_added+1;
                if abs(occtmp(i)-occtmp(i-1)-RRm) <= thRR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    occchk(i-1)=3;
                else
                    occchk(i-1)=0;
                end
                i=i+1;
            else %ho trovato un picco molto prossimo alla posizione presunta
                %cancello tutti gli altri trovati
                for aux=i-1:-1:j+1
                    if aux ~= best_idx
                        occtmp(aux)=60000; %marco per poi rimuovere
                        occchk(aux)=60000; %marco per poi rimuovere
                        count_removed = count_removed+1; %conto i rimossi
                    end
                end
                occtmp(occtmp==60000)=[]; %rimuovo quelli marcati
                occchk(occchk==60000)=[]; %rimuovo quelli marcati
                %sistemo i vettori occchk e occtmp
                if j+1<i-1
                    i=j+2;
                    if abs(occtmp(i)-occtmp(i-1)-RRm) <= thRR
                        occchk(i-1)=1;
                    else
                        occchk(i-1)=0;
                    end
                end
            end
        end
        i=i-1;
    end
    %coda se non inizia subito il vettore delle occorrenze
    while occtmp(1) > RRm
        nextR = occtmp(1) - RRm;
        if(nextR<starting_delay)
            break
        end
        occtmp=[nextR,occtmp];
        occchk=[3,occchk]; %col codice 3 identifico un picco che è statoaggiunto
        count_added = count_added+1;
    end

    end
end

