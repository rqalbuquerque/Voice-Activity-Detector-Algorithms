
% Zero Cross Rate (ZCR)

function y = zeroCrossRate(x)

	len = length(x);
	%y = sum(sign(abs(x(2:len)-x(1:len-1)))/(2*len));
	y = 0;
	for i=2:len
		y = y + abs(sign(x(i))-sign(x(i-1)));
	end
	y = y/(2*len);

end