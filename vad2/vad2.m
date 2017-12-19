
% VAD algorithm
%
% Based on:
% 	Authors: S. Kirill, V. Ekaterina and S. and Boris
% 	Title: Dynamical Energy-Based Speech/Silence Detector for Speech Enhancement Applications
%

function vadResult = vad2(x, fs, len)
    
    debug = 0;

	%% vars
	frameTime = 0.01;		% 10ms
	frameLen = 0;			% samples
	nFrames = 0;			% total frames without
	nSilenceTime = 0.05;		% 100ms of initial silence interval
	nSilenceSamples = 0;
    percOverlap = 0;

	lambda = 0.975;
	threshold = 0;
	factorEmin = 0.25;
	EminInit = 0;
	Emin = 0;
	Emax = 0;
	deltaIni = 1.001;
	delta = deltaIni;
	hangthreshold = 8;
    minUttLen = 5;
    
    % initial adjust
    nIgnoredSamples = find(abs(x)>0,1,'first') - 1;
    nSilenceSamples = ceil(nSilenceTime*fs);
	avgValue = sum(x(1+nIgnoredSamples:nSilenceSamples+nIgnoredSamples))/(nSilenceSamples);
    x = x - avgValue;
    
    
	%% segmentation
	if(frameTime > 0)
		frameLen = ceil(frameTime*fs);
	else 
		frameLen = 256;
    end
        
    frames = buffer(x,frameLen,round(frameLen*percOverlap));
    nFrames = size(frames,2);
	    
    %threshold = energy(x(1:nSilenceSamples))/nSilenceSamples;
    %threshold = rms_energy(x(nIgnoredSamples+1:nIgnoredSamples+nSilenceSamples));
    nSilenceFrames = floor((nSilenceSamples)/frameLen);
    rmseSilenceFrames = zeros(1,nSilenceFrames);
    for i=1:nSilenceFrames
        rmseSilenceFrames(i) = rms_energy(frames(i,:));
    end
	
    %% init vars
	Emax = max(rmseSilenceFrames);
	EminInit = factorEmin*min(rmseSilenceFrames);
	Emin = EminInit;
	lambda = (Emax-Emin)/Emax;
    threshold = sum(rmseSilenceFrames)/nSilenceFrames;
    
    %% debug
    thresholdArray = zeros(1,nFrames);
    EmaxArray = zeros(1,nFrames);
    EminArray = zeros(1,nFrames);
    energyArray = zeros(1,nFrames);
    
	%% whole processing
    lastVadFrames = zeros(1,minUttLen);
	vadFrames = zeros(1,nFrames);
	inactiveFrameCount = 0;
	for i=2:nFrames
		
		% frame
		frame = frames(:,i);
		rmse = rms_energy(frame);
		delta = delta*1.001;
		
		% updating vars
		if(rmse > Emax)
			Emax = rmse;
		end
		
		if(rmse < Emin)
			
			if(abs(rmse)<=eps)
				Emin = EminInit;
			else 
				Emin = rmse;
			end
			
			delta = deltaIni;
			
		else
			delta = deltaIni;
		end
		
		% update threshold
		threshold = (1-lambda)*Emax + lambda*Emin;
		
		% vad decisions
		if(rmse > threshold)
			vadFrames(i) = 1;
			inactiveFrameCount = 0;
        else
            
            % hangover
            countLastResults = sum(lastVadFrames);
            if(countLastResults>=minUttLen && inactiveFrameCount==0 && vadFrames(i-1) == 1)
                vadFrames(i) = 1;
                inactiveFrameCount = 1;
            elseif(inactiveFrameCount == hangthreshold)
                inactiveFrameCount = 0;
            elseif(inactiveFrameCount > 0)
				vadFrames(i) = 1;
				inactiveFrameCount = inactiveFrameCount+1;
			end
		
        end
		
        lastVadFrames = [lastVadFrames(1:minUttLen-1) vadFrames(i)];
        
		% increase minimum energy
		Emin = delta*Emin;
		
        energyArray(i) = rmse;
        EmaxArray(i) = Emax;
        EminArray(i) = Emin;
        thresholdArray(i) = threshold;
        
    end
	
    %% debug
    if debug==1
        figure(1),subplot(2,2,1),plot(energyArray),axis([0 nFrames 0 0.2]),title('Energy');
        subplot(2,2,2),plot(thresholdArray),axis([0 nFrames 0 0.2]),title('threshold');
        subplot(2,2,3),plot(EmaxArray),axis([0 nFrames 0 0.2]),title('Emax');
        subplot(2,2,4),plot(EminArray),axis([0 nFrames 0 0.2]),title('Emin');
    end
    
    %% results
	vadResult = zeros(1,len);
    vadFrames = join_short_intervals(vadFrames,3,nFrames);
	for i=1:nFrames
		vadResult((frameLen*(i-1))+1:min(frameLen*i,len)) = vadFrames(i);
	end
	
end