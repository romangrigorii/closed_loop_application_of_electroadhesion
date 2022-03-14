function out = filt_noise(length,center,sigma)
x = 1:length;
out = exp(-((x - center).^2)/(2*(sigma^2)));
end