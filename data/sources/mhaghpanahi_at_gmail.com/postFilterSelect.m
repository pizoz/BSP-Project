function [FECG,order] = postFilterSelect(Fout,blocksize) 

F1=LPFilter(Fout,45/1000);
F2=BPFilter(Fout,15/1000,45/1000);

F_order=zeros(1,4);
avg_kurt=zeros(1,4);
for i=1:4
    choice=[Fout(i,:)' F1(i,:)' F2(i,:)'];
    [F_order(i),avg_kurt(i)] = FSelect(choice,blocksize);
end

compare= hist(F_order,1:3);
[dummy, order]=sort(compare,'descend');

if (sum(avg_kurt(F_order==order(1)))<sum(avg_kurt(F_order==order(2))))
    dummy=order(1);
    order(1)=order(2);
    order(2)=dummy;
end


switch order(1)
    case 1
        FECG = Fout;
    case 2
        FECG = F1;
    case 3
        FECG = F2;
end

end