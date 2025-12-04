function out = correctFQRS(FECG, FQRS)
% FQRS correction function
%   (c) Jakub Kuzilek
% 	http://bio.felk.cvut.cz/~kuziljak/     E-mail: jakub.kuzilek@gmail.com
%   Version: 1.0                    Last update:   24/08/2013.
% 					                (Version: 1.0, 23/08/2013)
%
%======================================================
%
% PURPOSE:    This function detects FQRS in ECG recording. Developed for
%             purpose of Physionet Challenge 2013.
%
% MANDATORY INPUT ARGUMENTS
%   FECG ..... fetal ECG Mx1, M - length of data (one lead)
%   FQRS ..... fetal HRV Mx1, M - length of data (one lead)
% OPTIONAL INPUT ARGUMENTS
%   none
% OUTPUT ARGUMENTS
%   out .... positions of FQRS
%=========================================================
indx = find(findMeanHigh(diff(FQRS),10));

FHR = diff(FQRS);
FHR(indx) = [];

meanFHR = mean(FHR);

passed = [];

it = 1;
% DECISION OF SELECTED is problem or PREV or NEXT is problem?
while length(indx)> length(FQRS)*0.1 && it < 1000
    
    pri = indx(1);
    if pri > 1 && pri < length(FQRS)
        prev = FQRS(pri) - FQRS(pri-1);
        next = FQRS(pri+1) - FQRS(pri);
        %         if (prev < 0.8*meanFHR || prev > 1.2*meanFHR) && (next < 0.8*meanFHR || next > 1.2*meanFHR)
        if prev < 0.8*meanFHR || prev > 1.2*meanFHR 
            pos = FQRS(pri)-meanFHR;
            if pos > 1 && pos <= length(FECG)
                [~,lok] = max(FECG(max([pos-100 1]):min([pos+100 length(FECG)])));
                pos = pos-100+lok;
                FQRS(pri-1) = round(pos);
                indx = find(findMeanHigh(diff(FQRS),10));
                passed = [passed FQRS(pri)];
                indx(findPassed(FQRS(indx),passed)) = [];
            elseif pos<=1
                FQRS(1:pri-1) = [];
            end
        elseif next < 0.8*meanFHR || next > 1.2*meanFHR
            pos = FQRS(pri)+meanFHR;
            if pos - 100 > 1 && pos+100 <= length(FECG)
                [~,lok] = max(FECG(pos-100:pos+100));
                pos = pos-100+lok;
                FQRS(pri+1) = round(pos);
                indx = find(findMeanHigh(diff(FQRS),10));
                passed = [passed FQRS(pri)];
                indx(findPassed(FQRS(indx),passed)) = [];
            end
        elseif (prev < 0.8*meanFHR || prev > 1.2*meanFHR) && (next < 0.8*meanFHR || next > 1.2*meanFHR)
            pos = (FQRS(pri-1)+FQRS(pri+1))/2;
            if pos - 100 > 1 && pos+100 <= length(FECG)
                [~,lok] = max(FECG(pos-100:pos+100));
                pos = pos-100+lok;
                FQRS(pri) = round(pos);
                indx = find(findMeanHigh(diff(FQRS),10));
                passed = [passed FQRS(pri)];
                indx(findPassed(FQRS(indx),passed)) = [];
            end
        end
    end
%     figure
%     plot(FECG)
%     hold on
%     stem(FQRS, FECG(FQRS),'r')
%     stem(FQRS(indx), FECG(FQRS(indx)),'g')
    
    it = it+1;
    
end

out = FQRS;


function out = findPassed(in,passed)
out = [];
for m = 1:length(in)
    if any(passed == in(m))
        out(end+1) = m;
    end
end
% for m = 1:length(indx)
%     ind = indx(m);
%     if ind > 1 && ind < length(FQRS);
%         next = FQRS(ind+1);
%         if next-meanFHR-100>1 && next-meanFHR+100 <= length(FECG)
%             [~,lok] = max(FECG(next-meanFHR-100:next-meanFHR+100));
%             FQRS(ind) = next-meanFHR-100+lok;
%         end
%     end
% end
%
% out = FQRS;