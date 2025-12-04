function [Tq std_rr iqr_rr med_rr ampl] = MyBuildSeq_3(S, Fs, bs)
%
% function: attempts to build sequence of annotatios
%

v_mu    = 350:10:480;
v_sigma = 10:5:50;
rr_iqr  = zeros(numel(v_mu)*numel(v_sigma), 1); 
rr_med  = zeros(numel(v_mu)*numel(v_sigma), 1); 
err     = zeros(numel(v_mu)*numel(v_sigma), 1); 
ampl    = zeros(numel(v_mu)*numel(v_sigma), 1); 
v       = zeros(numel(v_mu)*numel(v_sigma), 2); 
j = 1;
for mj=1:numel(v_mu)
    for sj=1:numel(v_sigma)
        Tq = BuildChain(S, Fs, v_mu(mj), v_sigma(sj));
        if( isempty(Tq)==0 )
            rr        = diff(Tq);
            rr_iqr(j) = iqr(rr);
            rr_med(j) = median(rr);
            
            err(j)    = std(rr); 
            ampl(j)   = std(S(Tq));
            v(j,1)    = mj;
            v(j,2)    = sj;
            j = j+1; 
        end
    end
end

rr_iqr(j:end) = [];
rr_med(j:end) = [];
err(j:end)    = [];
ampl(j:end)   = [];
v(j:end, :)   = [];

[~, I]  = min(err);
sigma       = v_sigma( v(I,2) );
mu          = v_mu( v(I,1) );
Tq          = BuildChain(S, Fs, mu, sigma);

std_rr = err(I);
iqr_rr = rr_iqr(I);
med_rr = rr_med(I);
ampl   = ampl(I);