
% VAD algorithm
%
% Based on:
% 	Authors: Y. Ma e A. Nishihara
% 	Title: Efficient voice activity detection algorithm using long-term	spectral flatness measure
%

function vadResult = vad4(x, fs, len)
    
    debug = 0;

	%% vars
	frameTime = 0.02;		% ms
	frameLen = 0;			% samples
    nFrames = 0;
    threshold = 0;
	percOverlap = 0.5;		% overlapping
	nFrames = 0;			% total frames without
    nfft = 512;             % fft order
	nSilenceTime = 0.1;		% 100ms of initial silence interval
    LN_initialFrames = 100; 
    lambda = 0.1;

    M = 10;
    R = 30;
    
	window = [];
    L = [];
    S = [];
    vadResult = [];
    bufferSpeech = [];
    bufferNonSpeech = [];
    
    %% initial adjust
    nSilenceSamples = ceil(nSilenceTime*fs);
    avgValue = sum(x(1:nSilenceSamples))/nSilenceSamples;
    x = x - avgValue + eps;
    
    %% repeat noise interval to better results
    %x = [repmat(x(1:nSilenceSamples),10,1); x; repmat(x(end-nSilenceSamples+1:end),10,1)];
    
	%% init vars
	if(frameTime > 0)
		frameLen = ceil(frameTime*fs);
	else 
		frameLen = 512;
    end
	
	%% window
	window = hanning(frameLen);
	
	%% get ovelaped frames
	frames = buffer(x,frameLen,round(frameLen*percOverlap));
	nFrames = size(frames,2);
    
    %% get power spectrum estimation
    S = welch_bartlett_spectrum_all_frames(frames,window,M,nfft,fs);
    
    %% get long-term spectral flatness measure
    L = lsfm_all_frames(frames,S,M,R,nfft,fs);
    
    %% initial decision
    LN = L(R+M:R+M+LN_initialFrames)';
    threshold = min(LN);
    vadFrames = zeros(nFrames,1);
    bufferNonSpeech = LN;
    bufferSpeech = LN;
    threshFrames = zeros(nFrames,1);
    for i=R+M+LN_initialFrames+1:nFrames
        
        threshold = lambda*min(bufferNonSpeech) + (1-lambda)*max(bufferSpeech);
        threshFrames(i) = threshold;
        lsfm = L(i);
        
        if(lsfm < threshold)
           vadFrames(i) = 1; 
           bufferSpeech = [bufferSpeech(2:LN_initialFrames) lsfm];
        else
           bufferNonSpeech = [bufferNonSpeech(2:LN_initialFrames) lsfm]; 
        end
        
    end
    
    % debug
    if debug
        subplot(4,1,1),plot(x),axis([1 length(x) min(x) max(x)]),title('signal');
        subplot(3,1,2),plot(L),axis([1 length(L) min(L) max(L)]);
        hold on,plot(threshFrames),axis([1 length(L) min(L) max(L)]),title('LSFM - Threshold');
        hold off,subplot(3,1,3),plot(vadFrames),axis([1 nFrames 0 1.2]),title('Initial VAD Decisions');
    end
    
    %% results
	vadResult = zeros(1,len);
    overlapLen = frameLen*percOverlap;
	for i=1:nFrames
		vadResult((overlapLen*(i-1))+1:min(overlapLen*i,len)) = vadFrames(i);
    end

end

function S = welch_bartlett_spectrum_all_frames(frames,window,M,nfft,fs)

    nFrames = size(frames,2);
    ks = nfft*(500/fs);
    ke = nfft*(4000/fs);
    S = zeros(ke-ks+1,nFrames);
    
    if(M<=nFrames)
        
        for i=M:nFrames
           
            s = zeros(ke-ks+1,1);
            
            for j=i-M+1:i

                f = frames(:,j);
                amp = abs(fft(f.*window,nfft));
                s = s + abs(amp(ks:ke)).^2;
                
            end
            
            S(:,i) = s/M;
        end
  
    end

end
  
function L = lsfm_all_frames(frames,S,M,R,nfft,fs)

    nFrames = size(frames,2);
    ks = nfft*(500/fs);
    ke = nfft*(4000/fs);
    L = zeros(nFrames,1);
        
    for i=R+M:nFrames

        am = zeros(ke-ks+1,1);
        gm = ones(ke-ks+1,1);

        for j=i-R+1:i
            am = am + S(:,j);
            gm = gm.*S(:,j);
        end

        am = am/R;
        gm = gm.^(1.0/R);   

        L(i) = sum(log10(gm./am));
    end
    
end

  























