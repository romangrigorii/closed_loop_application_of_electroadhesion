%*************************************************************************%
% This piece of code is to be sued to process impulse responses of the
% tribometer and plot the reusulting frequency response curve. Additional
% calculations and plots include the filtered signal that will be read and
% processed by the micro controller.
%*************************************************************************%
%%
load('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\nor freq response\impulse_responses.mat');
addpath('C:\Users\atrox\Desktop\Work\Research\Code\General code\MATLAB code');
samp = length(imps);
res = 2^16;
ff1 = zeros(samp,res); % original fft response of the tribometer
ff2 = zeros(samp,res); % fft response of the tribometer + low pass filter

pha1 = zeros(samp,res);
pha2 = zeros(samp,res);

FFT1 = zeros(samp,res);
FFT2 = zeros(samp,res);

sr = 100000;
f = (0:res)/res*sr;
f = f(1:end-1);

%% computing filter paramters and the delays associated with filters
cofr1 = 200;
cofr2 = 1000;

ac = [1 2*(2*pi*cofr1*.5) (2*pi*cofr1)^2];
bc = [0 0 (2*pi*cofr1)^2];
[bd,ad] = stoz(bc,ac,sr); 

%% computing FFTs and saving for post processing
for s = 1:samp
    sig1 = -detrend(imps{s});
    ff1(s,:) = fft(sig1*2/length(sig1),res);
    sig2 = filter(bd,ad,sig1);
    ff2(s,:) = fft(sig2*2/length(sig2),res);
end

ff1 = ff1./abs(vectomat(mean(ff1(:,round((10:50)/sr*res)),2),res));
ff2 = ff2./abs(vectomat(mean(ff2(:,round((10:50)/sr*res)),2),res));

FFT1 = (20*log10(abs(ff1)));
FFT2 = (20*log10(abs(ff2)));

pha1 = unwrap(angle(ff1.')).'*180/pi + 180;
pha2 = unwrap(angle(ff2.')).'*180/pi + 180;

%% plotting
a = figure;
a.Position(3) = 600;
a.Position(4) = 175;
tsize = 10;
hold on
m1 = mean(10.^(FFT1/20));
m2 = mean(10.^(FFT2/20));

a = plot(f,20*log10(m1));
a.Color = [.7 .7 .7];
a.LineWidth = 1;
a = plot(f,20*log10(m2),':');
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

a = legend('impulse response (IR)','digital low-pass IR');
a.Location = 'northwest';
a.EdgeColor = [1 1 1];
a.Position(2) = a.Position(2);
grid on

a = figure;
hold on
a.Position(3) = 600;
a.Position(4) = 150;
a = plot(f,mean(pha1));
a.Color = [.6 .6 .6];
a.LineWidth = 1;
a = plot(f,mean(pha2),':');
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

a = legend('impulse response (IR)','digital low-pass IR');
a.EdgeColor = [1 1 1];
a.Location = 'southwest';
a.Position(2) = a.Position(2)-.05;
grid on
