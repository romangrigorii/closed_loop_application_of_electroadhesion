%*************************************************************************%
% This piece of code is to be sued to process impulse responses of the
% tribometer and plot the reusulting frequency response curve. Additional
% calculations and plots include the filtered signal that will be read and
% processed by the micro controller.
%*************************************************************************%
%%
%load('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\lat freq response\impulse_responses.mat');
%addpath('C:\Users\atrox\Desktop\Work\Research\Code\General code\MATLAB code');
samp = length(imps);
res = 2^16;
ff1 = zeros(samp,res); % original fft response of the tribometer
ff2 = zeros(samp,res); % fft response of the tribometer + low pass filter
ff3 = zeros(samp,res); % fft response of the tribometer + special filter
ff4 = zeros(samp,res); % fft response of the tribometer + special filter adjusted
ffx = zeros(samp,res); % fft response of the tribometer + special filter adjusted

pha1 = zeros(samp,res);
pha2 = zeros(samp,res);
pha3 = zeros(samp,res);
pha4 = zeros(samp,res);
phax = zeros(samp,res);

gd1 = zeros(samp,res);
gd2 = zeros(samp,res);
gd3 = zeros(samp,res);
gd4 = zeros(samp,res);

FFT1 = zeros(samp,res);
FFT2 = zeros(samp,res);
FFT3 = zeros(samp,res);
FFT4 = zeros(samp,res);

sr = 100000;
f = (0:(res-1))/(res)*sr;
plot_phase = 1;
plot_group = 1;
plot_error_bars = 0;
%% computing filter paramters and the delays associated with filters
cofr1 = 5;
cofr2 = .05;
dcofr1 = 3300;

ac1 = [1 2*(2*pi*cofr1*.2) (2*pi*cofr1)^2];
bc1 = [1 0 0];
[bd1,ad1] = stoz(bc1,ac1,sr);

ac2 = [1 2*(2*pi*dcofr1*.5) (2*pi*dcofr1)^2];
bc2 = [1 2*(2*pi*dcofr1*1) (2*pi*dcofr1)^2];
[bd2,ad2] = stoz(bc2,ac2,sr);

[bcf,acf] = tfdata(tf(bc1,ac1)*tf(bc2,ac2));
bcf = cell2mat(bcf);
acf = cell2mat(acf);
[bdf,adf] = stoz(bcf,acf,sr);

%% computing FFTs and saving for post processing
for s = 1:samp
    sig1 = -detrend(imps{s});
    ff1(s,:) = fft(sig1*2/length(sig1),res);
    sig2 = filter(bd1,ad1,sig1);
    ff2(s,:) = fft(sig2*2/length(sig2),res);
    sig3 = filter(bd2,ad2,sig2);
    ff3(s,:) = fft(sig3*2/length(sig3),res);  
    sig4 = filter(bdf,adf,sig1);
    ff4(s,:) = fft(sig4*2/length(sig4),res); 
    %ffx(s,:) = ((10.^(interp1(lf,20*log10(mag),f)/20)).*interp1(lf,exp(-sqrt(-1)*ang),f)).';
end

ff1 = ff1./abs(vectomat(mean(ff1(:,round((20:1000)/sr*res)),2),res));
ff2 = ff2./abs(vectomat(mean(ff2(:,round((20:1000)/sr*res)),2),res));
ff3 = ff3./abs(vectomat(mean(ff3(:,round((20:1000)/sr*res)),2),res));
ff4 = ff4./abs(vectomat(mean(ff4(:,round((20:1000)/sr*res)),2),res));
ffx = ffx./abs(vectomat(mean(ffx(:,round((20:1000)/sr*res)),2),res));

FFTx = 20*log10(abs(ffx));
FFT1 = (20*log10(abs(ff1)));
FFT2 = (20*log10(abs(ff2)));
FFT3 = (20*log10(abs(ff3)));
FFT4 = (20*log10(abs(ff4)));

phax = unwrap(angle(ffx.')).'*180/pi;
pha1 = unwrap(angle(ff1.')).'*180/pi;
pha2 = unwrap(angle(ff2.')).'*180/pi;
pha3 = unwrap(angle(ff3.')).'*180/pi+45;
pha4 = unwrap(angle(ff4.')).'*180/pi+270;
%pha4((pha4>0)&(f>2000)) = -180-pha4((pha4>0)&(f>2000));

meanfft = mean(FFT1);
meanfft = interp1(f,meanfft,1:1000);
meanpha = mean(pha1);
meanpha = interp1(f,meanpha,1:1000);

save('C:\Users\atrox\Desktop\Work\Research\projects\Texture player\calibration\frequency responses\fr3.mat','meanfft','meanpha');

%% plotting
a = figure;
a.Position(3) = 600;
a.Position(4) = 175;
tsize = 10;
hold on
m1 = mean(10.^(FFT1/20));
m2 = mean(10.^(FFT2/20));
m3 = mean(10.^(FFT3/20));
m4 = mean(10.^(FFT4/20));
mx = mean(10.^(FFTx/20));
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
a.LineWidth = 1;
a = plot(f,20*log10(m4),'--');
a.Color = [1 0 0];
a.LineWidth = 1;

axis([1.5,5000,-25,25])
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

a = legend({'impulse response (IR)','IR * digital filter L(z)'},'Interpreter','latex','FontSize',tsize);
a.Location = 'northwest';
a.EdgeColor = [1 1 1];
a.Position(2) = a.Position(2) + .05;
grid on

a = figure;
hold on
a.Position(3) = 600;
a.Position(4) = 150;
a = plot(f,mean(pha1));
a.Color = [.6 .6 .6];
a.LineWidth = 1;
a = plot(f,mean(pha4),'--');
a.Color = [1 0 0];
a.LineWidth = 1;

axis([1.5,5000,-100,100])
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('phase (deg)','Interpreter','latex','FontSize',tsize);
set(gca, 'XScale', 'log')
set(gcf,'color','w');
yticks([-90 -60 -30 0 30 60 90])
xticks([1 10 100 1000 5000]);

a = legend({'impulse response (IR)','IR * digital filter L(z)'},'Interpreter','latex','FontSize',tsize);
a.EdgeColor = [1 1 1];
a.Location = 'northeast';
a.Position(2) = a.Position(2) + .05;
grid on
