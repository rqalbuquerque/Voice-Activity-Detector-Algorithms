% Apply a VAD based on PESQ VAD zeroing inactive speech samples

function vad_result = vadPESQ(signal,fs,len)
	
	global MINJOINTINTERVAL Downsample MINSPEECHLGTH JOINSPEECHLGTH MINUTTLENGTH

	fs_vis=8000;
	signal=resample(signal, fs_vis,fs);

	% const, some values from PESQ
	Downsample = 32; 
	MINSPEECHLGTH= 4;
	MINUTTLENGTH= 5;
	MINJOINTINTERVAL = 4;
	JOINSPEECHLGTH= 50;
	
	% VAD of PESQ
	lenSignal = length(signal);
	[VAD,logVAD] = apply_VAD(signal, lenSignal);
	% figure,plot(logVAD),title('logVAD');
	
	% Convert array logVAD to binLogVAD
	lenVAD = length(VAD);
	binLogVAD = logVad_to_binLogVAD(logVAD,lenVAD);
	% figure,plot(binLogVAD),title('binLogVAD');

	% Join utterances with short separation interval
	binLogVAD = join_utterances(binLogVAD,lenVAD);
	% figure,plot(binLogVAD),title('binLogVAD');
	
	% Find start, end and number of utterances
	% Based on the PESQ algorithm to find utterances
	[Nutterances, UttSearch_Start, UttSearch_End] = find_utterances_VAD( binLogVAD , lenVAD);
	
	% The innactive regions are updated to zero value
	vad_result = get_binary_results(signal,lenSignal,Nutterances,UttSearch_Start,UttSearch_End);

	% clean used variables
	clearvars -global MINJOINTINTERVAL Downsample MINSPEECHLGTH JOINSPEECHLGTH MINUTTLENGTH
end


function modified_binLogVAD = join_utterances(binLogVAD,nSamples)
	
	global MINJOINTINTERVAL
	
	modified_binLogVAD = binLogVAD(:);
    speech_flag = 0;
	last_end = 1;
	this_start = 1;
    for count= 1:nSamples
        value= binLogVAD(count);
		
        if( (value> 0) && (speech_flag== 0) ) 
            speech_flag= 1;
            this_start= count;
			
			if(this_start - last_end <= MINJOINTINTERVAL)
				modified_binLogVAD(last_end:this_start) = 1;
			end
        end
	
        if( ((value== 0) || (count == (nSamples))) && (speech_flag == 1) ) 
            speech_flag = 0;
			last_end = count;
        end
    end
		
end


function result = get_binary_results(signal,lenSignal,Nutterances,UttSearch_Start,UttSearch_End)

	global Downsample

	result = zeros(lenSignal,1);
	for i=1:Nutterances
		initUtt = (UttSearch_Start(i)-1)*Downsample + 1;
		endUtt = UttSearch_End(i)*Downsample;
		% next utterance
		result(initUtt:endUtt) = 1;
	end
end


function [Nutterances, UttSearch_Start, UttSearch_End] = find_utterances_VAD( VAD, VAD_length )    
    Utt_num = 1;
    speech_flag = 0;
    
    for count= 1:VAD_length
        VAD_value= VAD(count);
		
        if( (VAD_value> 0) && (speech_flag== 0) ) 
            speech_flag= 1;
            this_start= count;
            UttSearch_Start(Utt_num)= count;
        end
    
        if( ((VAD_value== 0) || (count == (VAD_length))) && (speech_flag == 1) ) 
            speech_flag = 0;
            UttSearch_End(Utt_num) = count;
			
            Utt_num= Utt_num + 1;            
        end
    end
    Nutterances = Utt_num-1;

end

function array_updated = logVad_to_binLogVAD(logVAD, nSamples)
	
	global MINUTTLENGTH
	
	array_updated = zeros(1,nSamples);
    speech_flag = 0;
    for count= 1:nSamples
        value= logVAD(count);
		
        if( (value> 0) && (speech_flag== 0) ) 
            speech_flag= 1;
            this_start= count;
        end
    
        if( ((value== 0) || (count == (nSamples))) && (speech_flag == 1) ) 
            speech_flag = 0;
			
			if(count - this_start < MINUTTLENGTH)
				array_updated(this_start:count) = 0;
			else
				array_updated(this_start:count) = 1;
			end
        end
    end
end

function [VAD, logVAD]= apply_VAD( signal, Nsamples)
    
    global Downsample MINSPEECHLGTH JOINSPEECHLGTH
	
    Nwindows= floor( Nsamples/ Downsample);
    %number of 4ms window
    
    VAD= zeros( 1, Nwindows);
    for count= 1: Nwindows
        VAD( count)= sum( signal( (count-1)* Downsample+ 1: ...
            count* Downsample).^ 2)/ Downsample;   
    end
    %VAD is the power of each 4ms window 
    
    LevelThresh = sum( VAD)/ Nwindows;
    %LevelThresh is set to mean value of VAD
    
    LevelMin= max( VAD);
    if( LevelMin > 0 )
        LevelMin= LevelMin* 1.0e-4;
    else
        LevelMin = 1.0;
    end
    %fprintf( 1, 'LevelMin is %f\n', LevelMin);
    
    VAD( find( VAD< LevelMin))= LevelMin;
    
    for iteration= 1: 12    
        LevelNoise= 0;
        len= 0;
        StDNoise= 0;    
    
        VAD_lessthan_LevelThresh= VAD( find( VAD<= LevelThresh));
        len= length( VAD_lessthan_LevelThresh);
        LevelNoise= sum( VAD_lessthan_LevelThresh);
        if (len> 0)
            LevelNoise= LevelNoise/ len;
            StDNoise= sqrt( sum( ...
            (VAD_lessthan_LevelThresh- LevelNoise).^ 2)/ len);
        end
        LevelThresh= 1.001* (LevelNoise+ 2* StDNoise);  
    end
    %fprintf( 1, 'LevelThresh is %f\n', LevelThresh);
    
    LevelNoise= 0;
    LevelSig= 0;
    len= 0;
    VAD_greaterthan_LevelThresh= VAD( find( VAD> LevelThresh));
    len= length( VAD_greaterthan_LevelThresh);
    LevelSig= sum( VAD_greaterthan_LevelThresh);
    
    VAD_lessorequal_LevelThresh= VAD( find( VAD<= LevelThresh));
    LevelNoise= sum( VAD_lessorequal_LevelThresh);
    
    if (len> 0)
        LevelSig= LevelSig/ len;
    else
        LevelThresh= -1;
    end
    %fprintf( 1, 'LevelSig is %f\n', LevelSig);
    
    if (len< Nwindows)
        LevelNoise= LevelNoise/( Nwindows- len);
    else
        LevelNoise= 1;
    end
    %fprintf( 1, 'LevelNoise is %f\n', LevelNoise);
    
    VAD( find( VAD<= LevelThresh))= -VAD( find( VAD<= LevelThresh));
    VAD(1)= -LevelMin;
    VAD(Nwindows)= -LevelMin;
    
    start= 0;
    finish= 0;
    for count= 2: Nwindows
        if( (VAD(count) > 0.0) && (VAD(count-1) <= 0.0) )
            start = count;
        end
        if( (VAD(count) <= 0.0) && (VAD(count-1) > 0.0) )
            finish = count;
            if( (finish - start)<= MINSPEECHLGTH )
                VAD( start: finish- 1)= -VAD( start: finish- 1);
            end
        end
    end
    %to make sure finish- start is more than 4
    
    if( LevelSig >= (LevelNoise* 1000) )
        for count= 2: Nwindows
            if( (VAD(count)> 0) && (VAD(count-1)<= 0) )
                start= count;
            end
            if( (VAD(count)<= 0) && (VAD(count-1)> 0) )
                finish = count;
                g = sum( VAD( start: finish- 1));
                if( g< 3.0* LevelThresh* (finish - start) )
                    VAD( start: finish- 1)= -VAD( start: finish- 1);
                end
            end
        end
    end
    
    start = 0;
    finish = 0;
    for count= 2: Nwindows
        if( (VAD(count) > 0.0) && (VAD(count-1) <= 0.0) )
            start = count;
            if( (finish > 0) && ((start - finish) <= JOINSPEECHLGTH) )
                VAD( finish: start- 1)= LevelMin;
            end        
        end
        if( (VAD(count) <= 0.0) && (VAD(count-1) > 0.0) )
            finish = count;
        end
    end
    
    start= 0;
    for count= 2: Nwindows
        if( (VAD(count)> 0) && (VAD(count-1)<= 0) )
            start= count;
        end
    end
    if( start== 0 )
        VAD= abs(VAD);
        VAD(1) = -LevelMin;
        VAD(Nwindows) = -LevelMin;
    end
    
    count = 4;
    while( count< (Nwindows-1) )
        if( (VAD(count)> 0) && (VAD(count-2) <= 0) )
            VAD(count-2)= VAD(count)* 0.1;
            VAD(count-1)= VAD(count)* 0.3;
            count= count+ 1;
        end
        if( (VAD(count)<= 0) && (VAD(count-1)> 0) )
            VAD(count)= VAD(count-1)* 0.3;
            VAD(count+ 1)= VAD(count-1)* 0.1;
            count= count+ 3;
        end
        count= count+ 1;
    end
    
    VAD( find( VAD< 0))= 0;
    
    % fid= fopen( 'mat_vad.txt', 'wt');
    % fprintf( fid, '%f\n', VAD);
    % fclose( fid);
    
    if( LevelThresh<= 0 )
        LevelThresh= LevelMin;
    end
    
    logVAD( find( VAD<= LevelThresh))= 0;
    VAD_greaterthan_LevelThresh= find( VAD> LevelThresh);
    logVAD( VAD_greaterthan_LevelThresh)= log( VAD( ...
        VAD_greaterthan_LevelThresh)/ LevelThresh);
    
end

