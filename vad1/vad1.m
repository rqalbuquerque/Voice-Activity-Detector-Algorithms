
% VAD algorithm
%
% Based on:
% 	Authors: R. G. Bachu, S. Kopparthi, B. Adapa and B. D. Barkana
% 	Title: Voiced/Unvoiced Decision for Speech Signals Based on Zero-Crossing Rate and Energy
%

function vadResult = vad1(x, fs, len)

    debug = 0;

	%% vars
	frameTime = 0.05;		% 50ms
	frameLen = 0;			% samples
	nFrames = 0;			% total frames without 
	nSilenceFrames = 1;     % init silence interval
    percAvg = 0.3;
	
	avgSilenceEnergy = 0;	% average energy of silence interval
	silencezeroCrossRate = 0;
	
	threshold1 = 0;
	thresholdzeroCrossRate = 0;
	thresholdEnergy = 0;
	factorThresh1 = 0.0006;
    factorThreshEnergy = 0.0012;
    factorThreshzeroCrossRate = 1.5;
	
	minIntervalBetweenUtts = 20;
	nUtterances = 0;
	searchWinVicinity = 1;
	maxLevelSubFrames = 6;
	
	window = [];			% window
	e = [];					% energy
	silenceEnergy = [];		% silence interval energy
	endPoints = [];
    
	if(frameTime > 0)
		frameLen = ceil(frameTime*fs);
	else 
		frameLen = 256;
    end 
    
    %% remove gain value of signal
    nSilenceSamples = ceil(nSilenceFrames*frameLen);
	avgValue = sum(x(1:nSilenceSamples))/(nSilenceSamples);
    x = x - avgValue + eps;
        
    %% segmentation
    window = hamming(frameLen);   
    frames = buffer(x,frameLen,0);
    nFrames = size(frames,2);
	
	%% End-point Detection
	e = zeros(1,nFrames);
    zcr = zeros(1,nFrames);
    for i=1:nFrames
		frame = frames(:,i);
		e(i) = energy(frame.*window);
        zcr(i) = zeroCrossRate(frame);
    end
    
    nIgnoredSamples = find(abs(x)>eps,1,'first') - 1;
    if nIgnoredSamples > 100
        percAvg = 0.1*percAvg;
    end
    
    avgSilenceEnergy = energy(x(1+nIgnoredSamples:nIgnoredSamples+nSilenceSamples).*window);
    % to avoid get the initials samples of speech with high value
    while avgSilenceEnergy > 1.2 && nIgnoredSamples>100
        avgSilenceEnergy = percAvg*avgSilenceEnergy;
    end
	silencezeroCrossRate = zeroCrossRate(x(1+nIgnoredSamples:nIgnoredSamples+nSilenceSamples));
    peakEnergy = max(e);
	threshold1 = avgSilenceEnergy + factorThresh1*(peakEnergy - avgSilenceEnergy);
	
    cross_inds = (e - threshold1)>0;
    %endPoints = [find(cross_inds,1) find(cross_inds,1,'last')];
	ini = 1;
    while ini < nFrames
        
        while(ini < nFrames && cross_inds(ini)<1)
            ini = ini+1;
        end
        
        fin = min(ini+1,nFrames);
        while(fin < nFrames && cross_inds(fin)>0)
           fin = fin+1; 
        end
        fin = fin-1;
        
        if(ini < fin)
            endPoints = [endPoints; ini fin];
        end
        
        ini = fin+1;
    end

	% refina a localizacao das elocucoes utilizando zeroCrossRate da vizinhanca
	nPoints = size(endPoints,1);
    for i=1:nPoints
        
		iniUtt = endPoints(i,1);
		finUtt = endPoints(i,2);
        
        searchWindowIni = max(iniUtt-searchWinVicinity,1);
        if i-1>0 
            searchWindowIni = max(endPoints(i-1,2)+1,searchWindowIni);
        end
		searchWindowFin = min(iniUtt+searchWinVicinity,nFrames);
        if (i+1)<=nPoints
            searchWindowFin = min(endPoints(i+1,1)-1,searchWindowFin);
        end
		zeroCrossRateVicinity = zeros(1,searchWindowFin-searchWindowIni+1);
		for j=searchWindowIni:searchWindowFin
			zeroCrossRateVicinity(j-searchWindowIni+1) = zeroCrossRate(frames(:,j));
        end
		strongValues = abs(silencezeroCrossRate-zeroCrossRateVicinity);
		if(max(strongValues) > 5*silencezeroCrossRate)
			indVicinity = find(max(strongValues)==strongValues);
            indexesFrames = searchWindowIni:searchWindowFin;
			endPoints(i,1) = indexesFrames(max(indVicinity));
		end

        searchWindowIni = max(finUtt-searchWinVicinity,1);
        if i-1>0 
            searchWindowIni = max(endPoints(i-1,2)+1,searchWindowIni);
        end
		searchWindowFin = min(finUtt+searchWinVicinity,nFrames);
        if i+1<=nPoints
            searchWindowFin = min(endPoints(i+1,1)-1,searchWindowFin);
        end
		zeroCrossRateVicinity = zeros(1,searchWindowFin-searchWindowIni+1);
		for j=searchWindowIni:searchWindowFin
			zeroCrossRateVicinity(j-searchWindowIni+1) = zeroCrossRate(frames(:,j));
        end
        indexesFrames = searchWindowIni:searchWindowFin;
		strongValues = abs(silencezeroCrossRate-zeroCrossRateVicinity);
		if(max(strongValues) > 5*silencezeroCrossRate)
			indVicinity = find(max(strongValues)==strongValues);
            indexesFrames = searchWindowIni:searchWindowFin;
			endPoints(i,2) = indexesFrames(min(indVicinity));
		end
		
    end
	
	%% vad decisions
	vadFrames = [];
	vadInds = [];
	thresholdzeroCrossRate = factorThreshzeroCrossRate*silencezeroCrossRate;
	thresholdEnergy = avgSilenceEnergy + factorThreshEnergy*(peakEnergy - avgSilenceEnergy);
    nUtterances = size(endPoints,1);
     
    % debug
    if(debug)
        figure(1), subplot(3,1,1),plot(x), title('signal');
        subplot(3,1,2),plot(e), title('energy'),axis([1 length(e) 0 max(e)]);
        subplot(3,1,2),hold on,plot((ones(1,length(e)))*thresholdEnergy);
        subplot(3,1,3),plot(zcr), title('zcr'),axis([1 length(e) 0 max(e)]);
        pause
    end
    
	for i=1:nUtterances
        
		for j=endPoints(i,1):endPoints(i,2)
            
			stack = [frameLen*(j-1)+1 min(frameLen*j,len)];
			
            while(~isempty(stack))
				ini = stack(1,1);
				fin = stack(1,2);
				stack(1:2) = [];
				frame = x(ini:fin);
                window = hamming(fin-ini+1);
                
				frameEnergy = energy(frame.*window(1:fin-ini+1));
				framezeroCrossRate = zeroCrossRate(frame);
				
				if(framezeroCrossRate < thresholdzeroCrossRate && frameEnergy > thresholdEnergy...
                        || frameEnergy>5*thresholdEnergy)
					vadFrames = [1; vadFrames];
					vadInds = [ini fin; vadInds];
				elseif(framezeroCrossRate > thresholdzeroCrossRate && frameEnergy < thresholdEnergy)
					vadFrames = [0; vadFrames];
					vadInds = [ini fin; vadInds];
				else
					if(length(frame) > 2^maxLevelSubFrames)
						half = ceil((fin-ini)/2);
						stack = [ini half+ini-1; half+ini fin];
					else
						if(frameEnergy > thresholdEnergy)
							vadFrames = [1; vadFrames];
							vadInds = [ini fin; vadInds];
						else
							vadFrames = [0; vadFrames];
							vadInds = [ini fin; vadInds];
						end
					end
				end
			end
			
		end
	
	end
	
	subFramesLen = length(vadInds);
	vadResult = zeros(1,len);
	for i=1:subFramesLen
		ini = vadInds(i,1);
		fin = vadInds(i,2);
		vadResult(ini:fin) = vadFrames(i);
	end
	
end