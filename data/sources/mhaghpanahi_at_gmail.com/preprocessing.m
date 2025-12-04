function Y=preprocessing(I,fs,wr)
% simple preprocessing of the raw data using OSET tools (slightly modified)
    [m,n]=size(I);
    
    if m>n % make sure I is a row matrix
        I=I';
        [m,n]=size(I);
    end
    
    Y=zeros(m,n);
    [num1,den1] = iirnotch(50/(fs/2),1/(fs/2));
    [num2,den2] = iirnotch(100/(fs/2),1/(fs/2));
    [num3,den3] = iirnotch(60/(fs/2),1/(fs/2));
    [num4,den4] = iirnotch(120/(fs/2),1/(fs/2));
    for j = 1:m,
        I(j,:) = filter(num1,den1,I(j,:));
        I(j,:) = filter(num2,den2,I(j,:));
        I(j,:) = filter(num3,den3,I(j,:));
        I(j,:) = filter(num4,den4,I(j,:));
        b = movMedian(I(j,:),wr*fs+1);
        b = movMedian(b,wr*fs+1);
        b = LPFilter(b,15/fs);

        Y(j,:) = I(j,:) - b;

    end
end