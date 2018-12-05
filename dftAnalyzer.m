[wave,fs]=audioread(filename,'native');
SOUND = wave(:,1);

%use fft and parseval's thm to get energy output
T = period;
numChunks = ceil(length(SOUND)/(fs*T));
leftover = mod(length(SOUND),(fs*T));
dftoutput = zeros([1 numChunks]);
for i = 1:1:numChunks %processing individual chunks
    if (i == numChunks) && (leftover ~= 0)
        tempifs = zeros([1 (fs*T)]);
        for k = 1:leftover
            tempifs(k) = SOUND((numChunks-1)*(fs*T)+k);
        end
    else
        tempifs = SOUND((i-1)*(fs*T)+1:i*(fs*T));
    end
    tempifs = abs(fft(tempifs)).^2;
    %apply weighting to first half of fft output vector
    for k = 1:1:(fs*T)
        tempifs(k) = tempifs(k)*10^(aweighting((k-1)*(fs/length(tempifs)))/20);
    end
    %sum first half of fft vector and then double it
    %also divide by vector size per parseval's thm
    dftoutput(i) = 2*sum(tempifs(1:(fs*T)))/(fs*T);
end

dftoutput = 10*log(dftoutput/1);
t = 0:T:(numChunks-1)*T;
plot(t,dftoutput);

function weightdB = aweighting(ft)
    Ra = (12200^2*ft^4)/((ft^2+20.6^2)*(ft^2+12200^2)*sqrt((ft^2+107.7^2)*(ft^2+737.9^2)));
    weightdB = 2 + 20*log(Ra);
end