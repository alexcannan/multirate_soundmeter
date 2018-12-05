fs = 25600;
%Generate band-pass filters
for i = 1:1:3
    fc(i) = (fs/2)*2^((i-1)/3);
    swSemi = 3; %semitone half-width of stopband
    pwSemi = 0.9; %semitone half-width of passband
    cenF = fc(i);
    Fst1 = cenF*2^(-swSemi/12);
    Fst2 = cenF*2^(swSemi/12);
    Fp1 = cenF*2^(-pwSemi/12);
    Fp2 = cenF*2^(pwSemi/12);
    ds = 40; %40 dB stop attenuation
    dp = 0.01; %0.01 passband ripple
    d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Fst1/fs,Fp1/fs,Fp2/fs,Fst2/fs,ds,dp,ds);
    hd(i) = design(d,'equiripple');
    %freqz(hd);
end

%output individual and summed response of b-p filters
% figure(1); clf;
% hold on;
order = 512;
fL1 = hd(1).Numerator; fL2 = hd(2).Numerator; fL3 = hd(3).Numerator;
bpf1 = abs(fft(hd(1).Numerator, order));
f = 0:fs/(order-1):fs;
bpf1db = 20*log(bpf1);
% plot(f,bpf1db);
bpf2 = abs(fft(hd(2).Numerator, order));
bpf2db = 20*log(bpf2);
% plot(f,bpf2db);
bpf3 = abs(fft(hd(3).Numerator, order));
bpf3db = 20*log(bpf3);
% plot(f,bpf3db);
title('Individual filter responses')
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([0 12800]);
ylim([-100 10]);
% figure(2); clf;
% plot(f,20*log(bpf1+bpf2+bpf3));
title('Summed filter response')
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([0 12800]);
ylim([-100 10]);

%Generate half-band LPF
LPFnum = firhalfband('minorder',0.48,0.0001);
%fvtool(LPFnum,'Fs',fs);

[wave,fsfile]=audioread(filename,'native');
SOUND = wave(:,1);

seconds = 0;
outputdata = analyzer(SOUND,hd,LPFnum,period);

%generate weight matrix
weights = zeros([7 3]);
for o = 7:-1:1
    for b = 1:3
        ft = 100*2^(o-1)*2^((b-1)/3);
        Ra = (12200^2*ft^4)/((ft^2+20.6^2)*(ft^2+12200^2)*sqrt((ft^2+107.7^2)*(ft^2+737.9^2)));
        weightdB(o,b) = 2 + 20*log(Ra);
    end
end

%apply weights to outputdata
for i = 1:size(outputdata,1)
for o = 7:-1:1
    for b = 1:3
        outputdata(i,o,b) = outputdata(i,o,b)*10^(weightdB(o,b)/20);
    end
end
end


%sum weighted output over time
multirateoutput = sum(sum(outputdata,3),2);
%dB scale relative to unity
multirateoutput = 10*log(multirateoutput/1);
% plot(multirateoutput);

function meter = analyzer(ifs,bpf,lpf,period)
    T = period;
    fs = 25600;
    numChunks = ceil(length(ifs)/(fs*T));
    leftover = mod(length(ifs),(fs*T));
    disp(strcat("Beginning analysis of ",num2str(numChunks)," chunk input"));
    if ~(isvector(ifs) && (length(ifs) > fs))
        error('Input must be vector with length > 25600')
    end
    meter = zeros([numChunks 7 3]);
    for i = 1:1:numChunks %processing individual chunks
        disp(strcat("Analyzing chunk ",num2str(i),"/",num2str(numChunks)));
        %fill temp chunk. Pad with zeros if not enough to fill.
        if (i == numChunks) && (leftover ~= 0)
            tempifs = zeros([1 (fs*T)]);
            for k = 1:leftover
                tempifs(k) = ifs((numChunks-1)*(fs*T)+k);
            end
        else
            tempifs = ifs((i-1)*(fs*T)+1:i*(fs*T));
        end
        for o = 7:-1:1  %processing octaves
            for b = 1:3 %processing bands
                meter(i,o,b) = (2^(7-o))*sum(conv(tempifs,bpf(b).Numerator,'same').^2);
            end
            tempifs = conv(tempifs,lpf,'same');
            tempifs = decimator(tempifs);
        end
    end
end



%function decimates input vector and returns vector of half length
function decim = decimator(ifs)
    lengthI = length(ifs);
    %disp(strcat("Discovered input of length ",num2str(lengthI)));
    decim = zeros([1 lengthI/2]);
    for d = 0:1:lengthI/2-1
        decim(d+1) = ifs(2*d+1);
    end
end


