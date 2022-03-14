dur = 10;
sr = 60000;
t = linspace(0,dur,dur*sr);
st = daq.createSession('ni');
st.Rate = sr;
st.DurationInSeconds = dur;
ch1 = addAnalogInputChannel(st,'Dev2','ai1','Voltage');
ch1.InputType = 'SingleEnded'; %% lateral force
ch2 = addAnalogInputChannel(st,'Dev2','ai9','Voltage');
ch2.InputType = 'SingleEnded'; %% lateral force
ch3 = addAnalogInputChannel(st,'Dev2','ai2','Voltage');
ch3.InputType = 'SingleEnded'; %% lateral force
ch4 = addAnalogOutputChannel(st,'Dev2','ao1','Voltage');
amp = 2.5;
sig = amp*sin(2*pi*50*t);
%sig = sigf*10;
sig(end) = 0;
queueOutputData(st,sig');
%% collect frict vs cur data
freqs = logspace(log10(20),log10(1000),10);
DATA = {};
r = randperm(length(freqs));
dur = 10;
sr = 60000;
t = linspace(0,dur+1,(dur+1)*sr);
st = daq.createSession('ni');
st.Rate = sr;
st.DurationInSeconds = dur;
ch1 = addAnalogInputChannel(st,'Dev2','ai1','Voltage');
ch1.InputType = 'SingleEnded'; %% lateral force
ch2 = addAnalogInputChannel(st,'Dev2','ai9','Voltage');
ch2.InputType = 'SingleEnded'; %% lateral force
ch3 = addAnalogInputChannel(st,'Dev2','ai2','Voltage');
ch3.InputType = 'SingleEnded'; %% lateral force
ch4 = addAnalogOutputChannel(st,'Dev2','ao1','Voltage');
aF = [0 0 (2*pi*3100).^2];
bF = [1 (2*pi*3100*.9) (2*pi*3100).^2];
aH = [1 (2*pi*1100*1.5) (2*pi*1100).^2];
bH = [0 0 (2*pi*1100).^2];
aJ = [1 (2*pi*5000)];
bJ = [0 (2*pi*5000)];


Of = tf(bJ,aJ)/(tf(bF,aF)*tf(bH,aH));
[bz,az] = tfdata(Of);
bz = cell2mat(bz);
az = cell2mat(az);
[bzd,azd] = stoz(bz,az,60000);
for i = 1:length(freqs)
    i
    freqs(r(i))
    sig = 5*sin(2*pi*freqs(7)*t);
    sigf = filter(bzd,azd,sig');
    sigf = sigf(sr/2+1:end-sr/2);
    sig = sig(sr/2+1:end-sr/2);
    sigf(end) = 0;
    queueOutputData(st,sigf);
    x = input('press enter \n');
    out = startForeground(st);
    DATA{r(i)} = out;
end
% save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\cur vs force\data2.mat');
%% processing collect frict vs cur data
sr = 60000;
[bl1,al1] = butter(2,1000/sr,'low'); % used to extract directional changes
[bl2,al2] = butter(2,5*2/sr,'high'); % used to remove low frequency kinematic components
[bl3,al3] = butter(3,2500*2/sr,'low'); % used to remove high frequency force components
[bn,an] = butter(2,10*2/sr,'low'); % used to remove low frequency kinematic components
magg = {};
angg = {};
a_n = {};
cors = {};

for i = 1:length(freqs)
    lat = DATA{i}(:,1)*.25;
    nor = DATA{i}(:,2);
    cur = DATA{i}(:,3);
    ii = 1;
    while abs(lat(ii))<.1
        ii = ii+10;
    end
    lat = lat(ii:end);
    nor = nor(ii:end);
    cur = cur(ii:end);
    latf = filtfilt(bl1,al1,lat);
    norf = filtfilt(bn,an,nor);
    norf = (norf - mean(norf(1:sr/2)))*(-.1526);
    lat1 = latf>0;
    lat1 = abs(diff(lat1)).*(1:length(lat1)-1)';
    lat1 = lat1(lat1~=0);
    
    [bl4,al4] = butter(2,freqs(i)*2/60000*[.7,1.3],'bandpass');
    latff = filtfilt(bl4,al4,lat);
    
    latf = filtfilt(bl2,al2,lat);
    latf = filtfilt(bl3,al3,latf);
    cors{i} = [];
    
    latff = latff.*sign(lat);
    latf = latf.*sign(lat);
    
    sig = 5*sin(2*pi*freqs(i)*t);
    sig = sig(sr/2+1:end-sr/2);
    
    for ii = 1:length(lat1)-1
        s = round((lat1(ii+1)-lat1(ii))*1/2)+lat1(ii);
        e = round((lat1(ii+1)-lat1(ii))*3/4)+lat1(ii);        
        lattf = latf(s:e);
        lattff = latff(s:e);
        norr = norf(s:e);
        sigg = cur(s:e);
        %env = envelope(lattff,round(sr/freqs(i)),'peak')./envelope(sigg,round(sr/freqs(i)),'peak');
        env = findpeaks(lattff)./5*100;       
        angg{i}(ii) = -acos(dot(sigg,lattff)/(norm(sigg)*norm(lattff)))*180/pi;
        magg{i}(ii) = mean(env);
    end        
end

%% plotting
load('C:\Users\atrox\Desktop\Work\Research\projects\Texture player\calibration\frequency responses\fr.mat');
load('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\filter design\low_pass_analog2.mat');

Sm = interp1(1:1000,10.^(meanfft/20),freqs);
Sp = interp1(1:1000,meanpha,freqs);
Lm = interp1(lf,mag,freqs);
Lp = interp1(lf,ang,freqs);


m = []; ms = []; a = []; as = [];
for i = 1:length(freqs)
    m(i) = mean(magg{i});
    ms(i) = std(magg{i});
    a(i) = mean(angg{i});
    as(i) = std(angg{i});
end

m = m./Sm./Lm;
ms =20*.434*m.*ms;
m = 20*log10(m);
a = a - Sp - Lp;
Pm = interp1(freqs,m,1:1000);
Pp = interp1(freqs,a,1:1000);
Sm = 10.^(meanfft/20);
Sp = meanpha;
Lm = interp1(freqs,Lm,1:1000);
Lp = interp1(freqs,Lp,1:1000);

%save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\TFs\P.mat','Pm','Pp');
%save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\TFs\S.mat','Sm','Sp');
%save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\TFs\L.mat','Lm','Lp');


p = figure();
p.Position(3) = 400;
p.Position(4) = 400;

subplot(2,1,1);
hold on
set(gca,'TickLabelInterpreter','latex','FontSize',11);
for i = 1:length(freqs)
    p = plot(freqs(i)*[1 1],[m(i)-ms(i),m(i)+ms(i)]);
    p.Color = [.5 .5 .5];
end
p = plot(freqs,m);
p.Color = [1 0 0];
set(gca, 'XScale', 'log');
% axis([18 1100 .01 .07]);
% yticks([.01 .03 .05 .07]);
axis([18 1100 -10 10]);
yticks([-9 -6 -3 0 3 6 9])
xticks([20 50 100 250 500 1000]);
ylabel('magnitude ($\frac{F}{R}$)','Interpreter','latex','FontSize',10);

subplot(2,1,2);
hold on
set(gca,'TickLabelInterpreter','latex','FontSize',10);
for i = 1:length(freqs)
    p = plot(freqs(i)*[1 1],[a(i)-as(i),a(i)+as(i)]);
    p.Color = [.5 .5 .5];
end
p = plot(freqs,a);
p.Color = [1 0 0];
set(gca, 'XScale', 'log');
% axis([18 1100 -30 0]);
% yticks([-30 -20 -10 0])
axis([18 1100 -210 30]);
yticks([-180 -120 -60 0]);
xticks([20 50 100 250 500 1000]);
xlabel('frequency (Hz)','Interpreter','latex','FontSize',10);
ylabel('phase (deg)','Interpreter','latex','FontSize',10);
set(gcf,'color','w');


%% lat compare
sr = 60000;
[bl1,al1] = butter(4,1000/sr,'low');
[bl2,al2] = butter(2,2*2/sr,'high');
[bl3,al3] = butter(2,2000*2/sr,'low');
[bn,an] = butter(2,2*5/sr,'low');
lat = out(:,1)*.25;
lat1 = filtfilt(bl1,al1,lat);
i = 1;
while abs(lat1(i))<.1
    i = i+10;
end
nor = out(i:end,2);
nor = nor - mean(nor(1:sr/2));
norf = filtfilt(bn,an,nor)*(-.1526);

ref = out(i:end,3);
lat = lat(i:end);
lat1 = lat1(i:end);
lat1 = lat1>0;
lat1 = abs(diff(lat1)).*(1:length(lat1)-1)';
lat1 = lat1(lat1~=0);


%latf = filtfilt(bl2,al2,lat);
latf = filtfilt(bl3,al3,lat);
latf = latf.*sign(lat);

for i = 1:length(lat1)-1
    s = round((lat1(i+1)-lat1(i))/4)+lat1(i);
    e = round((lat1(i+1)-lat1(i))*3/4)+lat1(i);
    latt = latf(s:e);
    reff = ref(s:e);
    norr = norf(s:e);
    st = 1;
    en = length(latt);
    while latt(st)<latt(st+1)
        st = st+1;
    end
    while latt(en)<latt(en-1)
        en = en-1;
    end    
    latt = latt(st:en);
    %latt = latt(end/4:end*3/4);
    reff = reff(st:en);
    %reff = reff(end/4:end*3/4);
    cors(i) = abs(corr(reff',latt));
end
    
%% look up table
st = 61000;
sr = 60000;
[b,a] = butter(1,2/sr*10,'low');
norf = filter(b,a,out(:,2))*(-.1526);
norf = norf(st:end) - mean(norf(1:st/2));
[b,a] = butter(2,2/sr*1300,'low');
lat = filtfilt(b,a,out(st:end,1))*.25;
ac3 = [1 2*(2*pi*cofr1*.5) (2*pi*cofr1)^2];
bc3 = [1 2*(2*pi*cofr2*.5) (2*pi*cofr2)^2];
[bd3,ad3] = stoz(bc3,ac3,sr); 
latf = filtfilt(bd3,ad3,lat);
cur = (out(st:end,3)+5)/10*8;

latf = filtfilt(bd3,ad3,lat);
lat1 = lat>0;
lat1 = abs(diff(lat1)).*(1:length(lat1)-1)';
lat1 = lat1(lat1~=0);
binn = {};
LAT = [];
NOR = [];
CUR = [];
for i = 1:length(lat1)-1
    latt = abs(lat(lat1(i)+100:lat1(i+1)-100));
    norr = norf(lat1(i)+100:lat1(i+1)-100);
    sigg = cur(lat1(i)+100:lat1(i+1)-100);
    st = 1;
    en = length(latt);
    while latt(st)<latt(st+1)
        st = st+1;
    end
    while latt(en)<latt(en-1)
        en = en-1;
    end    
    latt = latt(st:en);
    norr = norr(st:en);
    sigg = sigg(st:en);
    LAT = [LAT;latt];
    NOR = [NOR;norr];
    CUR = [CUR;sigg];
end

bin = {};
mout = []; sout = [];
respts = 51;
res = zeros(respts+1,2);
res(1:respts+1,1) = linspace(0,2,respts+1);
res(1:respts+1,2) = linspace(0,8,respts+1);
resm = (res(1:end-1,:)+res(2:end,:))/2;

for i = 1:respts
    for ii = 1:respts
        bin{i,ii} = LAT((NOR>res(i,1))&(NOR<=res(i+1,1))&(CUR>res(ii,2))&(CUR<=res(ii+1,2)));
        mout(i,ii) = mean(bin{i,ii});
        sout(i,ii) = std(bin{i,ii});
    end
end        

% plotting
data = mout(6:40,7:end)';
[n,m]=size(data);%assumes that d is a n x m matrix
[X,Y]=meshgrid(resm(6:40,1),resm(7:end,2));%your x-y coordinates
x = [];
x(:,1)=X(:); % x= first column
x(:,2)=Y(:); % y= second column
f=data(:); % your data f(x,y) (in column vector)
fun = @(c,x)c(1).*(x(:,1).^c(2)) + c(3).*((x(:,1).^c(4)).*x(:,2));
options=optimset('TolX',1e-6);
c0=[1 1 1 1];
[cc,g,h]=lsqcurvefit(fun,c0,x,f,[],[],options);
Ifit=fun(cc,x);
Ifit=reshape(Ifit,[n m]);%fitting data reshaped as matrix
a = surf(X,flip(Y),flip(Ifit),'FaceAlpha',.5);
a.EdgeColor = 'none';
hold on;
a = plot3(X, flip(Y), flip(data),'.');
for i = 1:35
a(i).Color = [.5 .5 .5];
end
R2 = 1 - sum(sum(h.^2))/(sum(sum((mean(mean(data))-data).^2)));

%% digital lockin
f = linspace(1,60000,120);
mag = zeros(1,200);
t = linspace(0,length(outt)/sr,length(outt));
for i = 1:length(f)
    sig1 = sin(2*pi*t*f(i));
    sig2 = cos(2*pi*t*f(i));
    mag(i) = sum(sqrt(((sig1.*outt(:,1)').^2)+((sig2.*outt(:,1)').^2)));
    i
end

%%
f1 = logspace(1,log10(3000),11);
f2 = logspace(log10(3000),log10(10000),16);
f3 = logspace(log10(10000),log10(60000),10);
%lf = [f1(1:end-1),f2(1:end-1),f3];
lf = logspace(1,3,10);
dur = 2;
sr = 125000;
st = daq.createSession('ni');
st.Rate = sr;
st.DurationInSeconds = dur;
t = linspace(0,dur,sr*dur);
ch1 = addAnalogInputChannel(st,'Dev2','ai1','Voltage');
ch1.InputType = 'SingleEnded'; %% lateral force
ch2 = addAnalogInputChannel(st,'Dev2','ai2','Voltage');
ch2.InputType = 'SingleEnded'; %% lateral force
ang = zeros(1,length(lf));
mag = ang;
sig = {};
for i = 1:length(lf)
    lf(i)
    i
    x = input('');
    out = startForeground(st);
    sig{i} = out;
    ang(i) = acos(dot(out(:,1),out(:,2))/(norm(out(:,1))*norm(out(:,2))))*180/pi;
    mag(i) = max(abs(fft(detrend(out(:,1)))))/max(abs(fft(detrend(out(:,2)))));
end
save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\misc\filter design\low_pass_analog2.mat','sig','mag','ang','sr','dur','lf')

%% plotting mag and phase response of LPF

tsize = 10;
a = figure;
a.Position(3) = 400;
a.Position(4) = 300;
subplot(2,1,1);
a = plot(lf,20*log10(mag),'.-');
a.Color = [0 0 0];
set(gca, 'XScale', 'log');
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('magnitude (dB)','Interpreter','latex','FontSize',tsize);
axis([9 75000 -40 5])
yticks([-40 -20 0])
xticks([10,100,1000,10000,50000])
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
subplot(2,1,2);
a = plot(lf,-ang*180/pi,'.-');
a.Color = [0 0 0];
set(gca, 'XScale', 'log');
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('phase (deg)','Interpreter','latex','FontSize',tsize);
axis([9 75000 -180 10])
yticks([-180 -135 -90 -45 0]);
xticks([10,100,1000,10000,50000])
set(gcf,'color','w');
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);

%% plotting ADC + DAC noise
% load('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\misc\noises\ADC_DAC noise\ADC_DACnoise.mat');

tsize = 11;
a = figure;
a.Position(3) = 250;
a.Position(4) = 200;
f = (0:(length(out)-1))/(length(out))*sr;
a = plot(f,2/length(out)*abs(fft(out)));
a.Color = [.7 .7 .7];
axis([10 1000 10^-10 1])
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('voltage (V)','Interpreter','latex','FontSize',tsize);
title('ADC + DAC noise','Interpreter','latex','FontSize',tsize);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gcf,'color','w');

%% plotting DAC noise
load('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\misc\noises\DAC noise\noise.mat');

tsize = 11;
% a = figure;
% a.Position(3) = 250;
% a.Position(4) = 200;
f = (0:(length(out)-1))/(length(out))*sr;
a = plot(f,abs(fft(2*out/length(out))));
a.Color = [.7 .7 .7];
axis([10 1000 10^-10 1])
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('voltage (V)','Interpreter','latex','FontSize',tsize);
title('DAC noise','Interpreter','latex','FontSize',tsize);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gcf,'color','w');
%% plotting lat sensors and amp noise
convFact = .25;
V = convFact * out(:,2);
tsize = 11;
a = figure;
a.Position(3) = 250;
a.Position(4) = 200;
f = (0:(length(V)-1))/(length(V))*sr;
a = plot(f,abs(fft(2*V/length(V))));
a.Color = [.7 .7 .7];
axis([10 1000 10^-10 1])
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('force (N)','Interpreter','latex','FontSize',tsize);
title('charge amp noise','Interpreter','latex','FontSize',tsize);
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gcf,'color','w');

%% computing and plotting unactuated screen correlations
%%%% save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\misc\noises\unact_screen\data.mat')
%%%% load('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\misc\noises\unact_screen\data.mat')
dur = 30;
sr = 100000;
t = linspace(0,dur,dur*sr);
st = daq.createSession('ni');
st.Rate = sr;
st.DurationInSeconds = dur;
ch1 = addAnalogInputChannel(st,'Dev2','ai1','Voltage');
ch1.InputType = 'SingleEnded'; %% lateral force
ch2 = addAnalogInputChannel(st,'Dev2','ai9','Voltage');
ch2.InputType = 'SingleEnded'; %% lateral force
x = input('press enter \n');
out = startForeground(st);

%[b,a] = butter(2,2/sr*1300,'low');
%out = filtfilt(b,a,out);

lat = detrend(out(:,1))*.25;
nor = ((out(:,2) - mean(out(1:sr/2,2)))*-.1526);

i = 1;
while abs(lat(i))<.05
    i = i + 10;
end
lat = lat(i:end);
nor = nor(i:end);

[bn,an] = butter(2,2*10/sr,'low');
nor = filtfilt(bn,an,nor);

lat1 = lat>0;
lat1 = abs(diff(lat1)).*(1:length(lat1)-1)';
lat1 = lat1(lat1~=0);

FFT = []; g = 1;
for i = 1:length(lat1)-1
    if (lat1(i+1)-lat1(i))>20000
        latt = abs(lat(lat1(i):lat1(i+1)));
        latt = latt(end/8:end*7/8);
        norr = nor(lat1(i):lat1(i+1));    
        norr = norr(end/8:end*7/8);
        FFT(g,1:(2^15)) = log10(abs(fft(latt./norr*.5,2^15))/length(latt)*2);
        g = g + 1;
    end    
end

f = (0:(2^15-1))/(2^15)*sr;
l = 1:length(lat1)-1;


hold on
a = plot(f,FFT(i,:));
a.Color = [.8 .8 .8];
set(gca, 'XScale', 'log');
axis([3 1000 -6 -1])
a = plot(f,mean(FFT));
a.Color = [1 0 0];
a = plot(f,mean(FFT)+std(FFT));
a.Color = [0 0 0];

for i = l(l~=1)
    a = plot(f,FFT(i,:));
    a.Color = [.8 .8 .8];
end
a = plot(f,mean(FFT));
a.Color = [1 0 0];
a = plot(f,mean(FFT)+std(FFT));
a.Color = [0 0 0];
a = legend({'a single swipe','$\mu$ of swipes','$\mu$ + $\sigma$ of swipes'},'Interpreter','latex','FontSize',11);
set(gcf,'color','w');
xlabel('frequency (Hz)','Interpreter','latex','FontSize',11);
ylabel('force (N)','Interpreter','latex','FontSize',11);
%a.Box = 'off';
grid on
%%
fun = @(x,xdata)(x(1)*exp(xdata*x(2))./xdata);
x = lsqcurvefit(fun,[.01,-.01],10:1000,max(FF(:,10:1000)));

a = figure;
a.Position(3) = 400;
a.Position(4) = 250;
hold on
tsize = 11;
a = plot(FF(1,:));
a.Color = [.8 .8 .8];
a = plot([10,2000],10*max(FF(:,10))./[10,2000]);
a.Color = [0 0 0];
a.LineWidth = 1.25;
a = plot(10:1000,x(1)*exp((10:1000)*x(2))./(10:1000));
a.Color = [1 0 0];
a.LineWidth = 1.25;

[j,k] = size(FF);
for i = 1:j
    a = plot(FF(i,:));
    a.Color = [.8 .8 .8];
end
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xticks([10 100 1000]);
yticks = ([.01 .001 .0001 .00001]);
a = plot([10,2000],10*max(FF(:,10))./[10,2000]);
a.Color = [0 0 0];
a.LineWidth = 1.25;
a = plot(10:1000,x(1)*exp((10:1000)*x(2))./(10:1000));
a.Color = [1 0 0];
a.LineWidth = 1.25;
axis([10 1000 .00001 .011])
a = legend({'$\frac{lateral force}{normal force}$','$\frac{1}{f}$ curve','rendering threshold curve'},'Interpreter','latex','FontSize',tsize);
a.Location = 'southwest';
a.Position(2) = a.Position(2) - .04;
a.Position(1) = a.Position(1) - .025;
a.EdgeColor = [1 1 1];
set(gcf,'color','w');
arrow([100.,.001],[200,.01],1)

%% 
tsize = 11;
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
a = figure;
a.Position(3) = 300;
a.Position(4) = 200;
hold on
a = histogram(pars(2,pars(2,:)>0)'./pars(1,pars(2,:)>0)');
a.FaceColor = [14 89 146]/256;
a.EdgeColor = [.5 .5 .5];
a = histogram(-pars(2,pars(2,:)<0)'./pars(1,pars(2,:)<0)');
a.FaceColor = [1 0 0];
a.EdgeColor = [.5 .5 .5];
set(gcf,'color','w');
axis([.2 .6 -25 300])

a = plot(median(pars(2,pars(2,:)>0)'./pars(1,pars(2,:)>0)'),-12.5,'.');
a.Color = [14 89 146]/256;
a.MarkerSize = 25;
a = plot(-median(pars(2,pars(2,:)<0)'./pars(1,pars(2,:)<0)'),-12.5,'.');
a.Color = [1 0 0];
a.MarkerSize = 25;
xlabel('lateral force / normal force','Interpreter','latex','FontSize',tsize);
ylabel('count','Interpreter','latex','FontSize',tsize);
a = legend({'$right \rightarrow left$','$left \rightarrow right$'},'Interpreter','latex','FontSize',tsize-2);
a.Box = 'off';


%% 
fun = @(x,xdata)(x(1)*xdata);
tsize = 11;
a = figure;
a.Position(3) = 300;
a.Position(4) = 300;

hold on
x1 = lsqcurvefit(fun,[0,0],pars(1,pars(2,:)>0)',pars(2,pars(2,:)>0)');
x2 = lsqcurvefit(fun,[0,0],pars(1,pars(2,:)<0)',-pars(2,pars(2,:)<0)');

a = plot(pars(1,pars(2,:)>0)',pars(2,pars(2,:)>0)','o');
a.MarkerSize = 4;
a.LineWidth = 1;
a.Color = [14 89 146]/256;
a = plot(pars(1,pars(2,:)<0)',-pars(2,pars(2,:)<0)','.');
a.Color = [1 0 0];
a = plot([0,2],x1(1)*[0,2]);
a.Color = [0 0 0];
a.LineWidth = 1.2;
a = plot([0,2],x2(1)*[0,2]);
a.Color = [0 0 0];
a.LineWidth = 1.2;

a = legend({'$right \rightarrow left$','$left \rightarrow right$'},'Interpreter','latex','FontSize',tsize-2);
axis([.2 1.3 .05 .5]);
a.Location = 'northwest';
a.Position(1) = a.Position(1);
a.Box = 'off';
xlabel('normal load (N)','Interpreter','latex','FontSize',tsize);
ylabel('lateral force (N)','Interpreter','latex','FontSize',tsize);
text(.75,.15,'$R^{2}_{right \rightarrow left} = .77$','Interpreter','latex','FontSize',tsize-2);
text(.75,.1,'$R^{2}_{right \rightarrow left} = .76$','Interpreter','latex','FontSize',tsize-2);
set(gcf,'color','w');

%% TF formulaiton
P = tf(.03);

aL = [1 2*(2*pi*5500*.25) (2*pi*5500)^2];
bL = [0 0 (2*pi*5500)^2];
L = tf(bL,aL);

aH = [1 2*(2*pi*5*.2) (2*pi*cofr1)^2];
bH = [1 0 0];
H = tf(bH,aH);

aS = [1 2*(2*pi*4400*.1) (2*pi*4400)^2];
bS = [0 0 (2*pi*4400)^2];
S = tf(bS,aS);

aT = [1 2*(2*pi*3000*.1) (2*pi*3000)^2];
bT = [0 0 (2*pi*3000)^2];
T = tf(bT,aT);

a = [1 2*(2*pi*100*.1) (2*pi*100)^2];
b = [0 2*(2*pi*100*.1) 0]*8300;

C = tf(b,a);
C2 = T/(P*(1 - (T*H*S*L)));
O = C*P/((C*P*H*L*S)+1);
O2 = C2*P/((C2*P*H*L*S)+1);

[k,l,p] = bode(C,logspace(1,5,100000));
[k2,l2,p2] = bode(C2,logspace(1,5,100000));

[ko,lo,po] = bode(O,logspace(1,5,100000));
[ko2,lo2,po2] = bode(O2,logspace(1,5,100000));

%% constructing NN
[bl1,al1] = butter(4,1000/sr,'low');
[bl2,al2] = butter(2,5*2/sr,'low');
[bl3,al3] = butter(3,5*2/sr,'low');
[bn,an] = butter(2,2*5/sr,'low');
lat = out(:,1)*.25;
lat1 = filtfilt(bl1,al1,lat);
i = 1;
while abs(lat1(i))<.1
    i = i+10;
end
nor = out(i:end,2);
nor = nor - mean(nor(1:sr/2));
norf = filtfilt(bn,an,nor)*(-.1526);

lat = lat(i:end);
lat1 = lat1(i:end);
lat1 = lat1>0;
lat1 = abs(diff(lat1)).*(1:length(lat1)-1)';
lat1 = lat1(lat1~=0);

latf = filtfilt(bl2,al2,lat);
latf2 = filtfilt(bl3,al3,lat);
%latf = latf.*sign(lat);
lat1 = [1,lat1'];
in = [];
out = [];
len = 100;
for i = 1:length(lat1)-1
    latt = latf(lat1(i):lat1(i+1));
    latt2 = latf2(lat1(i):lat1(i+1));
    
    
