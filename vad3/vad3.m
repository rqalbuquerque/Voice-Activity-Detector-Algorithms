
% VAD algorithm
%
% Based on:
% 	Authors: J. Ramirez, J. C. Segura, C. BenÃ­tez, Ã?. de la Torre e A. Rubio
% 	Title: Efficient voice activity detection algorithms using long-term speech information
%

function vadResult = vad3(x, fs, len)

    debug = 0;

	%% vars
	frameTime = 0.02;		% 20ms
	frameLen = 0;			% samples
    nFrames = 0;
	percOverlap = 0.5;		% 
	nFrames = 0;			% total frames without
    order = 6;
	nSilenceFrames = 5;     %
	nSilenceTime = 0.02;	% 100ms of initial silence interval
    
	nfft = 256;
    K = 3;
    alpha = 0.95;
    E0 = 1.5941e-07;        % cleanest
    E1 = 0.9519;            % noisest
    gamma0 = 6;
    gamma1 = 2.5;
    nonSpeechCount = 0;
    hangoverCount = 0;
    maxHangCount = 2;
    maxLtsd = 0;
    biasLtsd = 0;
	
	window = [];
    vadResult = [];

	%% init vars
    if(frameTime > 0)
		frameLen = ceil(frameTime*fs);
	else 
		frameLen = 1024;
    end
    
    %% initial adjust
    nSilenceSamples = ceil(nSilenceTime*fs);
	avgValue = sum(x(1:nSilenceSamples))/(nSilenceSamples);
    x = x - avgValue + eps;
    
    nIgnoredSamples = find(abs(x)>eps,1,'first') - 1;
    nIgnoredFrames = floor(nIgnoredSamples/(frameLen/2)+1);
    
	%% window
	window = hamming(frameLen);
	
	%% get ovelaped frames
	frames = buffer(x,frameLen,round(frameLen*percOverlap));
	nFrames = size(frames,2);
    
	%% initial avg noise
	N = avg_spectrum(frames,window,1,nSilenceFrames,nfft,nIgnoredFrames);
	
    %% Optimal threshold calculation
    [threshold] = calibration(frames,window,nSilenceFrames,E0,E1,gamma0,gamma1,nIgnoredFrames);
    
    %% vad decision
    vadFrames = zeros(1,nFrames);
    ltsdFrames = zeros(1,nFrames);
    threshFrames = ones(1,nFrames)*threshold;
    
    for i=order+1:nFrames-order
        amp = ltse(frames,window,i,order,nfft);
        ltsdFrames(i) = 10*log10( (1.0/(nfft+1)) * sum( (amp.^2)./(N.^2)) );
    end
    
    [pks,locs] = findpeaks(ltsdFrames);
    biasLtsd = inf;
    for i=2:length(locs)
        peak1 = locs(i-1);
        peak2 = locs(i);
        minBetweenPeaks = min(ltsdFrames(peak1:peak2));
        if peak1>20 && peak2<length(ltsdFrames)-20 && minBetweenPeaks < biasLtsd
            biasLtsd = minBetweenPeaks;
        end
    end
    
    ltsdFrames = ltsdFrames-biasLtsd;
    
    if debug == 1
        fig = figure(1),clf(1),subplot(2,1,1),plot(x);
        subplot(2,1,2),plot(ltsdFrames), axis([1 length(ltsdFrames) 0 max(ltsdFrames)]);
        hold on, subplot(2,1,2),plot(ones(1,length(ltsdFrames))*threshold);
    end
    
    for i=order+1:nFrames-order
        
        ltsd = ltsdFrames(i);
        hangtest = 0;
        
       % hangover rule
        if(ltsd > maxLtsd)
           maxLtsd = ltsd; 
        end
        
        if(ltsd > threshold)
            vadFrames(i) = 1;
            nonSpeechCount = 0;
            hangtest=1;
        else               
            nonSpeechCount = nonSpeechCount+1;
            
            if(maxLtsd < 25 && hangtest==1)
                
                if(hangoverCount == maxHangCount)
                   hangoverCount = 0; 
                   hangtest = 0;
                elseif(vadFrames(i-1) == 1 && hangoverCount==0)

                    hangoverCount = hangoverCount+1;

                    if(ltsd+2.5 > threshold)
                       vadFrames(i) = 1; 
                    end
                elseif(hangoverCount > 0)
                    hangoverCount = hangoverCount+1; 
                    if(ltsd+2.5 > threshold)
                       vadFrames(i) = 1; 
                    end
                end
                
            end
            
            % update avg noise spectrum
            if(nonSpeechCount>K)
                Nk = avg_spectrum(frames,window,i-K+1,i,nfft,0);
                N = alpha*N + (1-alpha)*Nk;
            end
            
        end
    end
    
    % debug
    if debug == 2
        fig = figure(1);
        subplot(3,1,1),plot(x),axis([1 length(x) min(x) max(x)]),title('signal');
        subplot(3,1,2),plot(vadFrames),axis([1 length(vadFrames) 0 1.3]),title('VAD');
        subplot(3,1,3),plot(ltsdFrames),axis([1 length(ltsdFrames) min(ltsdFrames) max(max(ltsdFrames),1)]),title('LTSD');
        hold on,plot(threshFrames);
        pause
        clf(fig);
    end
    
    %% results
	vadResult = zeros(1,len);
    overlapLen = frameLen*percOverlap;
	for i=1:nFrames
		vadResult((overlapLen*(i-1))+1:min(overlapLen*i,len)) = vadFrames(i);
	end
end

% função de calibracao do threshold que utiliza os valores de 
%  energia do sinal
function [gamma,avgE] = calibration(frames,window,n,E0,E1,gamma0,gamma1,nIgnoredFrames)
    
    avgE = 0;
    
    if n+nIgnoredFrames > size(frames,2)
        error('not input range valid!');
    end
    for i=nIgnoredFrames+1:n+nIgnoredFrames
        f = frames(:,i);
        E = energy(f.*window);
        avgE = avgE + E;        
    end
    
    avgE = avgE/n;
    
    if(avgE <= E0)
       gamma = gamma1; 
    elseif(avgE >= E1)
        gamma = gamma0;
    else
        gamma = ((gamma0-gamma1)/(E0-E1))*avgE + gamma0 - (gamma0-gamma1)/(1-(E1/E0));
    end
    
end

function maxAmp = ltse(frames,window,ind,order,nfft)

	len = size(frames,2);
	maxAmp = zeros(nfft+1,1);

    for i=max(1,ind-order):min(len,ind+order)

        f = frames(:,i);
        F = fft(f.*window,2*nfft);
        amp = abs(F(1:nfft+1));
        maxAmp = max(maxAmp,amp);

    end

end

function avg = avg_spectrum(frames,window,ini,fin,nfft,nIgnoredFrames)

	avg = zeros(nfft+1,1);
    len = size(frames,2);
    count = 0;
     
    ini = ini+nIgnoredFrames;
    fin = fin+nIgnoredFrames;

    if(ini>0 && fin<=len)
        
        for i=ini:fin
            
            f = frames(:,i);
            F = fft(f.*window,2*nfft);
            amp = abs(F(1:nfft+1));
            avg = avg + amp;
            
        end
        
        avg = avg/(fin-ini+1);
    else
       error('not valid input range!'); 
    end
    
end




