function modifiedValues = join_short_intervals(binValues,MINJOINTINTERVAL,nSamples)
	
	modifiedValues = binValues(:);
    speech_flag = 0;
	last_end = 1;
	this_start = 1;
    for count= 1:nSamples
        value= binValues(count);
		
        if( (value> 0) && (speech_flag== 0) ) 
            speech_flag= 1;
            this_start= count;
			
			if(this_start - last_end <= MINJOINTINTERVAL)
				modifiedValues(last_end:this_start) = 1;
			end
        end
	
        if( ((value== 0) || (count == (nSamples))) && (speech_flag == 1) ) 
            speech_flag = 0;
			last_end = count;
        end
    end
		
end