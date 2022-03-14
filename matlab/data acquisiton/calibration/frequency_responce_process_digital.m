%*************************************************************************%
% This piece of code is to be sued to process impulse responses of the
% tribometer and plot the reusulting frequency response curve. Additional
% calculations and plots include the filtered signal that will be read and
% processed by the micro controller.
%*************************************************************************%
%%
load('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\lat freq response\impulse_responses.mat');
addpath('C:\Users\atrox\Desktop\Work\Research\Code\General code\MATLAB code');

sr_orig = 100000;
sr = 20000;
samp = length(imps);
res = 2^14;

ff1 = zeros(samp,res); % original fft response of the tribometer
ff2 = zeros(samp,res); % fft response of the tribometer + low pass filter
ff3 = zeros(samp,res); % fft response of the tribometer + special filter
ff4 = zeros(samp,res); % fft response of the tribometer + special filter adjusted

pha1 = zeros(samp,res);
pha2 = zeros(samp,res);
pha3 = zeros(samp,res);
pha4 = zeros(samp,res);

gd1 = zeros(samp,res);
gd2 = zeros(samp,res);
gd3 = zeros(samp,res);
gd4 = zeros(samp,res);

FFT1 = zeros(samp,res);
FFT2 = zeros(samp,res);
FFT3 = zeros(samp,res);
FFT4 = zeros(samp,res);

f = (0:res)/res*sr;
f = f(1:end-1);
plot_phase = 1;
plot_group = 1;
plot_error_bars = 0;
%% computing filter paramters and the delays associated with filters
cofr1 = 2; % low frequency of the low-pass filter
cofr2 = .2;  % high frequency of the low-pass filter
dcofr1 = 3200; % remove

% low pass filter
ac1 = [1 2*(2*pi*cofr1*sqrt(2)) (2*pi*cofr1)^2];
bc1 = [1 0 0];
[bd1,ad1] = stoz(bc1,ac1,sr); 
% adjusted filter
%/ high pass component
ac2 = [1 2*(2*pi*cofr1*.5) (2*pi*cofr1)^2];
bc2 = [1 2*(2*pi*cofr2*.5) (2*pi*cofr2)^2];
[bd2,ad2] = stoz(bc2,ac2,sr); 

ac3 = [1 2*(2*pi*dcofr1*.35) (2*pi*dcofr1)^2];
bc3 = [1 2*(2*pi*dcofr1*1) (2*pi*dcofr1)^2];
[bd3,ad3] = stoz(bc3,ac3,sr);

[bcf,acf] = tfdata(tf(bc2,ac2)*tf(bc3,ac3));
bcf = cell2mat(bcf);
acf = cell2mat(acf);
[bdf,adf] = stoz(bcf,acf,sr);

%% computing FFTs and saving for post processing
for s = 1:samp
    sig1 = resample(-detrend(imps{s}),sr,sr_orig);
    ff1(s,:) = fft(sig1*2/length(sig1),res);
    sig2 = filter(bd1,ad1,sig1);
    ff2(s,:) = fft(sig2*2/length(sig2),res);
    sig3 = filter(bdf,adf,sig1);
    ff3(s,:) = fft(sig3*2/length(sig3),res);  
%     sig4 = filter(bdf,adf,sig1);
%     ff4(s,:) = fft(sig4*2/length(sig4),res); 
end

ff1 = ff1./abs(vectomat(mean(ff1(:,round((20:1000)/sr*res)),2),res));
ff2 = ff2./abs(vectomat(mean(ff2(:,round((20:1000)/sr*res)),2),res));
ff3 = ff3./abs(vectomat(mean(ff3(:,round((20:1000)/sr*res)),2),res));
ff4 = ff4./abs(vectomat(mean(ff4(:,round((20:1000)/sr*res)),2),res));

FFT1 = 20*log10(abs(ff1));
FFT2 = 20*log10(abs(ff2));
FFT3 = 20*log10(abs(ff3));
FFT4 = 20*log10(abs(ff4));
pha1 = unwrap(angle(ff1))*180/pi;
pha2 = angle(ff2)*180/pi;
pha3 = angle(ff3)*180/pi;
pha4 = angle(ff4)*180/pi;
meanfft = mean(FFT1);
meanfft = interp1(f,meanfft,1:1000);
meanpha = mean(pha1);
meanpha = interp1(f,meanpha,1:1000);

save('C:\Users\atrox\Desktop\Work\Research\projects\Texture player\calibration\frequency responses\fr.mat','meanfft','meanpha');

%% plotting
a = figure;
a.Position(3) = 600;
a.Position(4) = 250;
tsize = 10;
hold on
m1 = mean(10.^(FFT1/20));
m2 = mean(10.^(FFT2/20));
m3 = mean(10.^(FFT3/20));
m4 = mean(10.^(FFT4/20));
s = std(10.^(FFT1/20));
if (plot_error_bars)
    [mi1,mi2] = min(m(round(40*(res)/sr):round(1000*(res)/sr)));
    i = mi2;
    while 20*log10(m(i))<20*log10(mi1)+3
        i = i + 1;
    end
    a = plot([1,3000],[-.4455 -.4455],':');
    a.Color = [.6 .6 .6];
    a.LineWidth = .9;
    a = plot([1,3000],[.4238,.4238],':');
    a.Color = [.6 .6 .6];
    a.LineWidth = .9;
    
    a = plot(f,20*log10(abs(m1-s)));
    a.Color = [.65 .65 .65];
    a.LineWidth = .5;
    a = plot(f,20*log10(m1+s));
    a.Color = [.65 .65 .65];
    a.LineWidth = .5;
end
a = plot(f,20*log10(m1));
a.Color = [.7 .7 .7];
a.LineWidth = 1.2;
a = plot(f,20*log10(m2),':');
a.Color = [0 0 0];
a.LineWidth = 1.2;
a = plot(f,20*log10(m3),'--');
a.Color = [1 0 0];
a.LineWidth = 1.2;

axis([5,1200,-25,25])
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('amplitude (dB)','Interpreter','latex','FontSize',tsize);
set(gca, 'XScale', 'log')
set(gcf,'color','w');
yticks([-40 -30 -20 -10 0 10 20]);
xticks([1 10 100 1000 5000]);

a = plot([.1,10000],[-1.5,-1.5],':');
a.Color = [.6 .6 .6];
a.LineWidth = .9;
a = plot([.1,10000],[1.5,1.5],':');
a.Color = [.6 .6 .6];
a.LineWidth = .9;

a = legend('unfiltered response','hp filter','hp filter + phase adjustment');
a.Location = 'northwest';
a.Position(2) = a.Position(2) + .05;
grid on

a = figure;
hold on
a.Position(3) = 600;
a.Position(4) = 250;
a = plot(f,mean(pha1));
a.Color = [.6 .6 .6];
a.LineWidth = 1.2;
a = plot(f,mean(pha2),':');
a.Color = [0 0 0];
a.LineWidth = 1.2;
a = plot(f,mean(pha3),'--');
a.Color = [1 0 0];
a.LineWidth = 1.2;

axis([5,1200,-100,100])
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('phase (deg)','Interpreter','latex','FontSize',tsize);
set(gca, 'XScale', 'log')
set(gcf,'color','w');
yticks([-90 -60 -30 0 30 60 90])
xticks([1 10 100 1000 5000]);

a = legend('unfiltered response','hp filter','hp filter + phase adjustment');
a.Location = 'southwest';
grid on
