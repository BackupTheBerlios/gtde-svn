function corr = myXCorrEUSIPCO(signal1,signal2,maxLag)

    %%% Check input

    corr = zeros(1,2*maxLag+1);
    
    chunkSize = length(signal1)-maxLag;
    % Positive values of tau
    for xcPos = 2:maxLag+1,
        corr(maxLag+xcPos) = sum(signal1(xcPos:chunkSize-1+xcPos).*signal2(1:chunkSize));
    end
    % Negative values of tau
    for xcPos = 2:maxLag+1,
        corr(maxLag+2-xcPos) = sum(signal2(xcPos:chunkSize-1+xcPos).*signal1(1:chunkSize));
    end
    % At 0
    corr(maxLag+1) = sum(signal1(1:chunkSize).*signal2(1:chunkSize));

end