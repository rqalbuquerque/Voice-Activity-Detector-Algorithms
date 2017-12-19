
% RMS of energy

function y = rms_energy(x)

	len = length(x);
	y = sqrt(sum(x.*x)/len);

end