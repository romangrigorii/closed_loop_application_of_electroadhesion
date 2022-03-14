%% DATA COLLECTION
%% sampling force in active touch w open/closed loop
freqsP = logspace(log10(10),log10(1000),20);
DATAP = {};
r = randperm(length(freqsP));
dur = 10;
sr = 10000;
t = linspace(0,dur,dur*sr);
st = daq.createSession('ni');
st.Rate = sr;
st.DurationInSeconds = dur;
ch1 = addAnalogInputChannel(st,'Dev2','ai1','Voltage');
ch1.InputType = 'SingleEnded'; %% lateral force
ch2 = addAnalogInputChannel(st,'Dev2','ai10','Voltage');
ch2.InputType = 'SingleEnded'; %% normal force
ch3 = addAnalogInputChannel(st,'Dev2','ai2','Voltage');
ch3.InputType = 'SingleEnded'; %% reference force
ch4 = addAnalogInputChannel(st,'Dev2','ai9','Voltage');
ch4.InputType = 'SingleEnded'; %% current
ch5 = addAnalogOutputChannel(st,'Dev2','ao0','Voltage');
for i = 1:length(freqsP)
    fprintf('percentage finished: %1.2f working on frequency %3.0f\n',(i-1)/length(freqsP),freqsP(r(i)));
    sig = 1*sin(2*pi*freqsP(r(i))*t);
    sig(end) = 0;
    queueOutputData(st,sig');
    x = input('press enter \n');
    out = startForeground(st);
    DATAP{r(i)} = out;
end
save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\cur vs force\data\data33.mat','DATAP','sr','freqsP')
%% computing low pass filter response
% collecting data
freqsL = round(logspace(log10(10),log10(6000),20));
DATAL = {};
dur = 2;
sr = 60000;
t = linspace(0,dur,dur*sr);
st = daq.createSession('ni');
st.Rate = sr;
st.DurationInSeconds = dur;
ch1 = addAnalogInputChannel(st,'Dev2','ai2','Voltage');
ch1.InputType = 'SingleEnded'; %% lateral force
ch2 = addAnalogInputChannel(st,'Dev2','ai1','Voltage');
ch2.InputType = 'SingleEnded'; %% lateral force
for i = 1:length(freqsL)
    fprintf('percentage finished: %1.4f plug in frequency %3.4f\n',(i-1)/length(freqsL),freqsL(i));
    x = input('press enter \n');
    out = startForeground(st);
    DATAL{i} = out;
end
% processing data
Lmag = []; Lang = [];
for i = 1:length(freqsL)
    %[b,pl] = butter(1,2*freqsL(i)/sr*[.8,1.2],'bandpass');
    in = DATAL{i}(:,1);
    out = DATAL{i}(:,2);
    inf = abs(fft(in));
    outf = abs(fft(out));
    Lmag(i) = max(outf)/max(inf);
    Lang(i) = -acos(dot(in,out)/(norm(in)*norm(out)))*180/pi;
end
save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\filter design\low_pass_analog6.mat','DATAL','freqsL','sr','Lmag','Lang');

%% DATA PROCESSING
%% processing free exploration of the surface w open loop electroadhesion rendering in order to characterize P(s)
sr = 60000;
[bl1,al1] = butter(3,100/sr,'low');        % used for finding changes in direction
[bl2,al2] = butter(1,10*2/sr,'high');       % used for removing low frequency kinematic information
[bl3,al3] = butter(1,1000*2/sr,'low');      % used for removing high frequency noise
[bn,an] = butter(2,50*2/sr,'low');            % used for filtering out high frequency normal force components
lat_const = .25;
nor_const = -.1526;
load('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\cur vs force\data\data33.mat');
PGL = {};
PGLs = {};
PGLc = [];
for i = 1:length(freqsP)
    lat = DATAP{i}(:,1)*lat_const;
    lat = lat - mean(lat(sr/8:sr/2));
    sig = DATAP{i}(:,3)/100;
    cur = (DATAP{i}(:,4) + 5)/2;
    cur =  DATAP{i}(:,3)/2 + 2.5;
    nor = DATAP{i}(:,2)*nor_const;
    lat1 = filtfilt(bl1,al1,lat);
    ii = 1;
    while abs(lat1(ii))<.01
        ii = ii+1;
    end
    lat1 = lat1(ii:end);
    lat = lat(ii:end);
    latd = lat1;
    sig = sig(ii:end);
    cur = cur(ii:end);
    nor = (nor - mean(nor(sr/8:sr/2)));
    nor = nor(ii:end);
    norf = filtfilt(bn,an,nor);
    
    latl = filtfilt(bn,an,lat);
    lat1 = lat1>0;
    lat1 = abs(diff(lat1)).*(1:length(lat1)-1)';
    lat1 = lat1(lat1~=0);
    [blbp,albp] = butter(1,freqsP(i)*2/sr*[.7,1.3],'bandpass');    
    latfbp = filtfilt(blbp,albp,lat);
    curf = filtfilt(blbp,albp,cur);
    envlat = envelope(latfbp,round(sr/freqsP(i)/2/pi),'peak');
    envcur = envelope(curf,round(sr/freqsP(i)/2/pi),'peak');
    
    latf = filtfilt(bl2,al2,lat);
    latf = filtfilt(bl3,al3,latf);
    latfbp = latfbp.*sign(latd);
    latf = latf.*sign(latd);
    lat1 = lat1(2:end);
    k = 1;
    
    for ii = 1:length(lat1)-1
        s = round((lat1(ii+1)-lat1(ii))*1/4)+lat1(ii);
        e = round((lat1(ii+1)-lat1(ii))*3/4)+lat1(ii);
        latll = latl(s:e);
        lattf = latf(s:e);
        lattfbp = latfbp(s:e);
        norr = norf(s:e);
        curr = curf(s:e);
        envlatt = envlat(s:e);
        envcurr = envcur(s:e);
        envrat =  envlatt./envcurr;
        sigg = sig(s:e); 
            
        if length(curr)>6000 && sum(abs((cur(s:e)-2.5))>1.5)==0            
            u = 1;
            PGLs{i,1,k} = [];
            PGLs{i,2,k} = [];
            while u*600 <= length(envrat)
                PGLs{i,1,k} = [PGLs{i,1,k},mean(envrat((u-1)*600+1:u*600))];
                PGLs{i,2,k} = [PGLs{i,2,k},mean(norr((u-1)*600+1:u*600))];
                u = u + 1;
            end
            
            PGLc(i) = k;
            PGL{i,1}(k) = -acos(dot(lattfbp,sigg)/(norm(sigg)*norm(lattfbp)))*180/pi;
            PGL{i,2}(k) = mean(envlatt);
            PGL{i,3}(k) = mean(envlatt(1:round(length(envlatt)/2))./envcurr(1:round(length(envlatt)/2)));
            PGL{i,4}(k) = mean(envlatt(round(length(envlatt)/2):end)./envcurr(round(length(envlatt)/2):end));
            PGL{i,5}(k) = sign(mean(latll));
            k = k +1;
        end
        
    end
    i
end

%% processing free exploration of the surface w closed loop
% extraction process
plotss = 4;
dirs = {};
dirs{1} = 'C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\cur vs force\data\data20.mat';
dirs{2} = 'C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\cur vs force\data\data21.mat';
dirs{3} = 'C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\cur vs force\data\data22.mat';
dirs{4} = 'C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\cur vs force\data\data23.mat';
sr = 60000;
[bl1,al1] = butter(2,10/sr,'low');         % used for finding changes in direction
[bl2,al2] = butter(1,10*2/sr,'high');       % used for removing low frequency kinematic information
[bl3,al3] = butter(1,1000*2/sr,'low');      % used for removing high frequency noise
[blc,alc] = butter(2,40*2/sr,'low');        % used for removing high frequency noise
[bn,an] = butter(2,2*50/sr,'low');            % used for filtering out high frequency normal force components
PGLm = {}; % mangitude of P transfer function
PGL = {};
PGLc = [];
PGLp = {}; % phase of P transfer function
lat_const = .25;
nor_const = -.1526;
Lsig = [];
FF = [];
FFN = [];
for p = 1:plotss
    load(dirs{p});
    for i = 1:length(freqsP)
        lat = DATAP{i}(:,1)*lat_const;
        nor = DATAP{i}(:,2)*nor_const;
        sig = DATAP{i}(:,3);
        cur = DATAP{i}(:,4);
        lat = lat - mean(lat(sr/8:sr/2));
        lat1 = filtfilt(bl1,al1,lat);
        ii = 1;
        while abs(lat1(ii))<.01
            ii = ii+10;
        end
        lat1 = lat1(ii:end);
        lat = lat(ii:end);

        sig = sig(ii:end);
        cur = cur(ii:end);
        
        norf = filtfilt(blc,alc,nor);
        norf = (norf - mean(norf(sr/8:sr/2)));
        norf = norf(ii:end);
 
        latl = filtfilt(bl1,al1,lat);
        lat1 = lat1>0;
        
        lat1 = abs(diff(lat1)).*(1:length(lat1)-1)';
        lat1 = lat1(lat1~=0);
        [blbp,albp] = butter(2,freqsP(i)*2/sr*[.75,1.25],'bandpass');
        latfbp = filtfilt(blbp,albp,lat);
        curf = filtfilt(blbp,albp,cur);
        curf2 = filtfilt(bn,an,cur);
        sigf = filtfilt(blbp,albp,sig);
        envlat = envelope(latfbp,round(sr/freqsP(i)/2),'peak');
        envcur = envelope(curf,round(sr/freqsP(i)/2),'peak');
        
        latf = filtfilt(bl2,al2,lat);
        latf = filtfilt(bl3,al3,latf);
        latfbp = latfbp.*sign(lat);
        latf = latf.*sign(lat);
        lat1 = lat1(2:end);
        k = 1;
        
        for ii = 1:length(lat1)-1
            s = round((lat1(ii+1)-lat1(ii))*1/4)+lat1(ii);
            e = round((lat1(ii+1)-lat1(ii))*3/4)+lat1(ii);
            latll = latl(s:e);
            lattf = latf(s:e);
            lattfbp = latfbp(s:e);
            norr = norf(s:e);
            sigg = sigf(s:e);
            curr = cur(s:e);
            currf = curf2(s:e);
            envlatt = envlat(s:e);
            envcurr = envcur(s:e);
            if (max(abs(curr))<4 && length(sigg)>6001)
                ff = detrend(envlatt./envcurr);
                ff = ff.*hann(length(ff));
                ff = abs(fft(ff,2^16))/length(ff);
                FF = [FF,ff(1:546)];
                %plot(ff(:500));
                %pause(.01);
                ff = detrend(norr);
                ff = ff.*hann(length(ff));
                ff = abs(fft(ff,2^16))/length(ff);
                FFN = [FFN,ff(1:546)];
                Lsig = [Lsig,length(sigg)];
                L = floor(length(sigg)/600)*600;
                V = L/600;
                tn = norr(1:L);
                tnm = []; tns = [];
                for v = 1:V
                    tnm(v) = mean(tn(((v-1)*600+1):v*600));
                    tns(v) = std(tn(((v-1)*600+1):v*600));
                end
                tl = latll(1:L).*sign(mean(lat(s:e)));
                tlm = [];
                for v = 1:V
                    tlm(v) = mean(tl(((v-1)*600+1):v*600));
                end             
                tc = envlatt(1:L)./envcurr(1:L)*2;
                tcm = []; tcs = [];
                for v = 1:V
                    tcm(v) = mean(tc(((v-1)*600+1):v*600));
                    tcs(v) = std(tc(((v-1)*600+1):v*600));
                end
                tcu = currf(1:L)/2;
                tcum = [];
                for v = 1:V
                    tcum(v) = mean(tcu(((v-1)*600+1):v*600));
                end
                PGL{p,i,k,1} = tnm;
                PGL{p,i,k,2} = tlm;
                PGL{p,i,k,3} = tcm;
                PGL{p,i,k,4} = tcum;
                PGL{p,i,k,5} = sign(mean(latll));
                PGL{p,i,k,6} = tns;
                PGL{p,i,k,7} = tcs;             
                
                PGLc(p,i) = k;
                PGLp{p,i}(k) = -acos(dot(curr,lattfbp)/(norm(lattfbp)*norm(curr)))*180/pi;
                PGLm{p,i,1}(k) = mean(envlatt);
                PGLm{p,i,2}(k) = mean(envlatt(1:round(length(envlatt)/2)));
                PGLm{p,i,3}(k) = mean(envlatt(round(length(envlatt)/2):end));
                PGLm{p,i,4}(k) = mean(envlatt)*sign(mean(latbb));
                k = k +1;
            end
        end
    end
    p
end
%% plotting frequency relationship of normal load and effect variation
hold on
pl = plot(linspace(0,500,546),mean(FFN.')/max(mean(FFN.')))
pl.Color = [0 0 1];
pl.LineWidth = .75;
pl = plot(linspace(0,500,546),mean(FF.')/max(mean(FF.')))
pl.Color = [1 0 0];
pl.LineWidth = .75;
% pl = plot(linspace(0,500,546),(std(FFN.')+mean(FFN.'))/max(mean(FFN.')),':')
% pl.Color = [0 0 1];
% pl.LineWidth = .75;
for i = 1:103
    pl = plot(linspace(0,500,546),FF(:,i)/max(mean(FF.')))
    pl.Color = [.9 .9 .9];
    pl.LineWidth = .05;
end
pl = plot(linspace(0,500,546),mean(FF.')/max(mean(FF.')))
pl.Color = [1 0 0];
pl.LineWidth = .75;
% pl = plot(linspace(0,500,546),(std(FF.')+mean(FF.'))/max(mean(FF.')),':')
% pl.Color = [1 0 0];
% pl.LineWidth = .75;
%set(gca, 'YScale', 'log');
axis([0 250 0 2]);
yticks([0 .5 1]);
xticks([20 100 200]);
xlabel('frequency (Hz)','Interpreter','latex','FontSize',10);
pl = plot(linspace(0,500,546),mean(FFN.')/max(mean(FFN.')))
pl.Color = [0 0 1];
pl.LineWidth = .75;
pl = plot(linspace(0,500,546),mean(FF.')/max(mean(FF.')))
pl.Color = [1 0 0];
pl.LineWidth = .75;
pl = legend({'$W$','$P$'},'Interpreter','latex');
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex','FontSize',10);
%% plotting histograk for open loop effect strength

pl = figure;
hold on
pl.Position(3) = 350;
pl.Position(4) = 200;
vecn = [];
vecr = [];
cols = spring();
for i = 1:length(freqsP)
    for k = 1:PGLc(i)
        vecr = [vecr,PGLs{i,1,k}];
        vecn = [vecn,PGLs{i,2,k}];
    end
end

histogram(vecr,.05:.0001:.175);
axis([.05 .175 0 75]);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex','FontSize',10);
xlabel('P','Interpreter','latex','FontSize',10);
ylabel('count','Interpreter','latex','FontSize',10);
%% plotting normal force vs current envelope for closed loop case
pl = figure;
hold on
pl.Position(3) = 350;
pl.Position(4) = 200;
vecn = [];
vecl = [];
vecc = [];
vecu = [];
vecns = [];
veccs = [];
vecl = [];
cols = spring();
for p = 1:plotss
    for i = 1:length(freqsP)
        for k = 1:PGLc(p,i)
            vecn = [vecn,PGL{p,i,k,1}];
            vecl = [vecl,PGL{p,i,k,2}];
            vecc = [vecc,PGL{p,i,k,3}];
            vecu = [vecu,PGL{p,i,k,4}];
            vecns = [vecns,PGL{p,i,k,6}];
            veccs = [veccs,PGL{p,i,k,7}];
        end
    end
end

figure
hold on
NN = 80;
CC = 80;
MAT = zeros(NN,CC);
vecnn = [];
veccn = [];
for n = 1:NN
    for c = 1:CC
        MAT(n,c) = sum((vecl>=(.4*c/(CC+1))).*(vecl<(.4*((c+1)/(CC+1)))).*(vecn>=(2*n/(NN+1))).*(vecn<(2*((n+1)/(NN+1)))));
        if MAT(n,c)>25
            pl = plot((2*((n+.5)/(NN+1))),.2*(c+.5)/(CC+1),'.');
            pl.MarkerSize = (MAT(n,c)^.5)/2 + 1;
            pl.Color = [0 0 0];
            veccn = [veccn,vecl(logical((vecl>=(.4*c/(CC+1))).*(vecl<(.4*((c+1)/(CC+1)))).*(vecn>=(2*n/(NN+1))).*(vecn<(2*((n+1)/(NN+1))))))];
            vecnn = [vecnn,vecn(logical((vecl>=(.4*c/(CC+1))).*(vecl<(.4*((c+1)/(CC+1)))).*(vecn>=(2*n/(NN+1))).*(vecn<(2*((n+1)/(NN+1))))))];
        elseif MAT(n,c)>0
            pl = plot((2*((n+.5)/(NN+1))),.2*(c+.5)/(CC+1),'o');
            pl.MarkerSize = 2;
            pl.Color = [.8 .8 .8];
        end
    end
end
l = MAT>0;
l = MAT(:).*l(:);
l = l(l>0);

fun = @(c,x)(c(1).*(x.^c(2)));
[f,r,t] = lsqcurvefit(fun,[.04,1],vecnn,veccn);

% pl = figure;
% hold on
% pl.Position(3) = 300;
% pl.Position(4) = 300;
nvec = .01:.01:2;
% 
% pl = plot(vecn,vecc,'.');

pl = plot(nvec,f(1).*(nvec.^f(2)));
%pl = plot(nvec,f(1).*(exp(nvec*f(2))));
pl.Color = [1 0 0];
pl.LineWidth = 1.5;

pl = plot(median(vecn),.010,'x');
pl.MarkerSize = 6;
pl.LineWidth = 1.5;
pl.Color = [1 0 0 .5];

pl = plot(.15,median(vecc),'x');
pl.MarkerSize = 6;
pl.LineWidth = 1.5;
pl.Color = [1 0 0 .5];

axis([.15 1.25 .010 .14]);
set(gcf,'color','w');
xlabel('W (N)','Interpreter','latex','FontSize',10);
ylabel('$P$ (N/mA)','Interpreter','latex','FontSize',10);
xticks([.25 .5 .75 1 1.25])
yticks([.03 .06 .09 .12])
set(gca,'TickLabelInterpreter','latex','FontSize',10);
%text(.8,.045,'r = -.68','Interpreter','latex','FontSize',13);
%text(.8,.04,'$\rho $ = -.72','Interpreter','latex','FontSize',13);
text(.85,.09,'$R^{2}= .62$','Interpreter','latex','FontSize',10);
text(.85,.105,'$P = .041W^{-.63}$','Interpreter','latex','FontSize',10);
set(gcf,'renderer','Painters')


figure
hold on
NN = 80;
CC = 40;
MAT = zeros(NN,CC);
vecnn = [];
veccn = [];
for n = 1:NN
    for c = 1:CC
        MAT(n,c) = sum((vecc>=(.2*c/(CC+1))).*(vecc<(.2*((c+1)/(CC+1)))).*(vecn>=(2*n/(NN+1))).*(vecn<(2*((n+1)/(NN+1)))));
        if MAT(n,c)>11
            pl = plot((2*((n+.5)/(NN+1))),.2*(c+.5)/(CC+1),'.');
            pl.MarkerSize = (MAT(n,c)^.5)/2 + 1;
            pl.Color = [0 0 0];
            veccn = [veccn,vecc(logical((vecc>=(.2*c/(CC+1))).*(vecc<(.2*((c+1)/(CC+1)))).*(vecn>=(2*n/(NN+1))).*(vecn<(2*((n+1)/(NN+1))))))];
            vecnn = [vecnn,vecn(logical((vecc>=(.2*c/(CC+1))).*(vecc<(.2*((c+1)/(CC+1)))).*(vecn>=(2*n/(NN+1))).*(vecn<(2*((n+1)/(NN+1))))))];
        elseif MAT(n,c)>0
            pl = plot((2*((n+.5)/(NN+1))),.2*(c+.5)/(CC+1),'o');
            pl.MarkerSize = 2;
            pl.Color = [.8 .8 .8];
        end
    end
end
fun = @(c,x)(c(1).*(x.^c(2)));
%fun = @(c,x)(c(1).*(exp(x*c(2))));
[f,r,t] = lsqcurvefit(fun,[.04,-2],vecnn,veccn);

% pl = figure;
% hold on
% pl.Position(3) = 300;
% pl.Position(4) = 300;
nvec = .01:.01:2;
% 
% pl = plot(vecn,vecc,'.');

pl = plot(nvec,f(1).*(nvec.^f(2)));
%pl = plot(nvec,f(1).*(exp(nvec*f(2))));
pl.Color = [1 0 0];
pl.LineWidth = 1.5;

pl = plot(median(vecn),.010,'x');
pl.MarkerSize = 6;
pl.LineWidth = 1.5;
pl.Color = [1 0 0 .5];

pl = plot(.15,median(vecc),'x');
pl.MarkerSize = 6;
pl.LineWidth = 1.5;
pl.Color = [1 0 0 .5];

axis([.15 1.25 .010 .14]);
set(gcf,'color','w');
xlabel('W (N)','Interpreter','latex','FontSize',10);
ylabel('$P$ (N/mA)','Interpreter','latex','FontSize',10);
xticks([.25 .5 .75 1 1.25])
yticks([.03 .06 .09 .12])
set(gca,'TickLabelInterpreter','latex','FontSize',10);
%text(.8,.045,'r = -.68','Interpreter','latex','FontSize',13);
%text(.8,.04,'$\rho $ = -.72','Interpreter','latex','FontSize',13);
text(.85,.09,'$R^{2}= .62$','Interpreter','latex','FontSize',10);
text(.85,.105,'$P = .041W^{-.63}$','Interpreter','latex','FontSize',10);
set(gcf,'renderer','Painters')
%% effects of time
pglt1 = [];
pglt2 = [];

for p = 1:plotss
    for i = 1:length(freqsP)
        pglt1 = [pglt1,mean(PGL{p,i,1,3})];
        pglt2 = [pglt2,mean(PGL{p,i,PGLc(p,i),3})];
    end
end
[h,p] = ttest2(pglt1,pglt2);
%% computing swipe direction and swipe location effects on strength of the effect
meanF = []; meanS = []; meanL = []; meanR = []; TS = []; TE = [];
for i = 1:length(freqs)
    flag1 = 1;
    flag2 = 1;
    meanF(i) = mean(PGLm{i,2});
    meanS(i) = mean(PGLm{i,3});
    iL = 0;
    iR = 0;
    meanR(i) = 0;
    meanL(i) = 0;
    for ii = 1:length(PGLm{i,4})
        if (PGLm{i,4}(ii)<0)
            iR = iR + 1;
            meanR(i) = meanR(i) - PGLm{i,4}(ii);
            if flag1
                TS(1,i) = -PGLm{i,4}(ii);
                flag1 = 0;
            else
                TE(1,i) = -PGLm{i,4}(ii);
            end
        else
            iL = iL + 1;
            meanL(i) = meanL(i) + PGLm{i,4}(ii);
            if flag2
                TS(2,i) = PGLm{i,4}(ii);
                flag2 = 0;
            else
                TE(2,i) = PGLm{i,4}(ii);
            end
        end
    end
    meanR(i) = meanR(i)/iR;
    meanL(i) = meanL(i)/iL;
end

%% computing P, L, G
% loading G and L transfer function characteristics in order to compute P
load('C:\Users\atrox\Desktop\Work\Research\projects\Texture player\calibration\frequency responses\fr.mat');
load('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\filter design\low_pass_analog4.mat');
Sm = interp1(1:1000,10.^(meanfft/20),freqsP);
Sp = interp1(1:1000,meanpha,freqsP);
Lm = interp1(freqsL,Lmag,freqsP);
Lp = interp1(freqsL,Lang,freqsP);
m = []; ms = []; a = []; as = [];
plotss = 1;
for p = 1:plotss
    for i = 1:length(freqsP)
        m(p,i) = mean(PGL{i,3}(1:end-1));
        ms(p,i) = std(PGL{i,3});
        a(p,i) = mean(PGL{i,1});
        as(p,i) = std(PGL{i,1});
    end
end
%m = m./Sm./Lm;
a = a - Sp - Lp;
% Pm = interp1(freqsP,m,20:1000,'linear','extrap');
% Pp = interp1(freqsP,a,20:1000,'linear','extrap');
% Sm = 10.^(interp1(1:1000,meanfft,20:1000,'linear','extrap')/20);
% Sp = interp1(1:1000,meanpha,20:1000,'linear','extrap');
% Lm = interp1(freqsP,Lm,20:1000,'linear','extrap');
% Lp = interp1(freqsP,Lp,20:1000,'linear','extrap');
% save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\TFs\P.mat','Pm','Pp');
% save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\TFs\S.mat','Sm','Sp');
% save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\TFs\L.mat','Lm','Lp');

%% plotting P response
%pl = figure();
%pl.Position(3) = 400;
%pl.Position(4) = 300;
%subplot(2,1,1);
hold on
set(gca,'TickLabelInterpreter','latex','FontSize',10);
plotss = 1;
cols = winter();
shaps = {'s','o','d','*'};
for p = 1:plotss
    for i = 1:length(freqsP)
        for ii = 1:PGLc(p,i)
            if PGL{i,5}(ii)<0
                pl = plot(freqsP(i)*1.03,PGL{i,3}(ii),'.');
                pl.Color = [.6 .6 .6];
                %pl.MarkerSize = 3;
                %pl.LineWidth = .5;
                %pl.Color = [cols(round(p/plotss*50),:)];
                %p.MarkerSize = 3;
            else
                pl = plot(freqsP(i)*.97,PGL{i,3}(ii),'.');
                pl.Color = [.6 .6 .6];
                %pl.MarkerSize = 3;
                %pl.LineWidth = .5;
                %pl.Color = [cols(round(p/plotss*50),:)];
                %p.MarkerSize = 3;
            end
        end
    end
    pl = plot(freqsP,m(p,:));
    pl.Color = [1 0 0];
    %pl.Color = [cols(round(p/plotss*50),:),.7];
    pl.LineWidth = 1;
    %     pl = plot(1200,.01*p,shaps{p});
    %     pl.Color = [.6 .6 .6];
    %     pl.MarkerSize = 6;
    %     pl.LineWidth = 1.1;
end
set(gca, 'XScale', 'log');
%axis([18 1200 .05 .15]);
yticks([.07 .1 .13 .16]);
axis([9 1200 .06 .16]);
xticks([10 100 1000]);
ylabel('magnitude (N)','Interpreter','latex','FontSize',10);
xlabel('frequency (Hz)','Interpreter','latex','FontSize',10);
% subplot(2,1,2);
% hold on
% set(gca,'TickLabelInterpreter','latex','FontSize',10);
% for p = 1:plotss
%     pl = plot([1,10000],-[1,1]*(p-1)*45,'--');
%     pl.Color = [0,0,0];
%     for i = 1:length(freqsP)
%         pl = plot(freqsP(i)*[1 1],[pl(p,i)-as(p,i),pl(p,i)+as(p,i)]-((plotss-p)*45));
%         pl.Color = [.5 .5 .5];
%         %pl.Color = [cols(round(p/plotss*64),:)];
%         pl.LineWidth = .8;
%     end
%     pl = plot(freqsP,pl(p,:)-((plotss-p)*45));
%     %pl.Color = [1 0 0];
%     pl.Color = [cols(round(p/plotss*50),:),.7];
%     pl.LineWidth = 1;
%     %     pl = plot(1200,-(plotss-p)*45,shaps{p});
%     %     pl.Color = [.6 .6 .6];
%     %     pl.MarkerSize = 6;
%     %     pl.LineWidth = 1.1;
% end
% set(gca, 'XScale', 'log');
% axis([18 1200 -285 15]);
% yticks([-270 -180 -90 0])
% xticks([20 75 250 1000]);
% ylabel('phase (deg)','Interpreter','latex','FontSize',10);
% set(gcf,'color','w');

%% plotting L response and esimate
load('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\filter design\low_pass_analog6.mat');

aL = [1 2*(2*pi*6600*.68) (2*pi*6600)^2];
bL = [0 0 (2*pi*6600)^2];
L = tf(bL,aL);
[k,l,~] = bode(L,freqsL*2*pi);

tsize = 10;
pl = figure;
pl.Position(3) = 400;
pl.Position(4) = 300;
subplot(2,1,1);
hold on
pl = plot(freqsL,20*log10(Lmag),'.-');
pl.Color = [.8 .8 .8];
pl = plot(freqsL,20*log10(squeeze(k)),'.-');
pl.Color = [1 0 0];
set(gca, 'XScale', 'log');
ylabel('magnitude (dB)','Interpreter','latex','FontSize',tsize);
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
subplot(2,1,2);
hold on
pl = plot(freqsL,Lang,'.-');
pl.Color = [.8 .8 .8];
pl = plot(freqsL,squeeze(l),'.-');
pl.Color = [1 0 0];
set(gca, 'XScale', 'log');
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('phase (deg)','Interpreter','latex','FontSize',tsize);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);

%% plotting G response and estimate
aG = [1 2*(2*pi*4400*.6) (2*pi*4400)^2];
bG = [0 0 (2*pi*4400)^2];
G = tf(bG,aG);

[k,l,p] = bode(G,logspace(0,5,10000));
tsize = 10;
pl = figure;
pl.Position(3) = 400;
pl.Position(4) = 300;
subplot(2,1,1);
hold on
pl = plot(squeeze(p)/2/pi,20*log10(squeeze(k)));
pl.Color = [1 0 0];
set(gca, 'XScale', 'log');
ylabel('magnitude (dB)','Interpreter','latex','FontSize',tsize);
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
axis([10 1000 -5 5])
yticks([-5 0 5])
subplot(2,1,2);
hold on
pl = plot(squeeze(p)/2/pi,squeeze(l));
pl.Color = [1 0 0];
set(gca, 'XScale', 'log');
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('phase (deg)','Interpreter','latex','FontSize',tsize);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
axis([10 1000 -10 1])
yticks([-10 -5 0])


%% computing DC - 10kHz response of L,G,C,T
load('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\lat freq response\impulse_responses.mat');
% computing G
L = length(imps);
FFT = [];
sr = 100000;
pad = 2^16;
f = (0:pad)/pad*sr;
f = f(1:end-1);
mm = [];
r1 = round(10*pad/sr);
r2 = round(1000*pad/sr);
for i = 1:L
    FFT(i,1:pad) = fft(imps{i},pad);
    mm(i) = mean(abs(FFT(i,r1:r2)));
    FFT(i,1:pad) = FFT(i,1:pad);%/mm(i);
end
Gmag = 20*log10(abs(FFT));
Gmeanmag = mean(Gmag);
Gstdmag = .434*std(abs(FFT))/mean(abs(FFT));
Gang = angle(FFT);
for i = 1:L
    Gang(i,:) = unwrap(Gang(i,:));
    if mean(Gang(i,r1))>pi/2
        Gang(i,:) = Gang(i,:) - pi;
    end
    if mean(Gang(i,r1))<-pi/2
        Gang(i,:) = Gang(i,:) + pi;
    end
end
Gang = Gang*180/pi;
Gmeanang = mean(Gang);
Gstdang = std(Gang);
aG = [1 2*(2*pi*4400*.44) (2*pi*4400)^2];
bG = [0 0 (2*pi*4400)^2];
G = tf(bG,aG);
[mG,pG,fG] = bode(G,f*2*pi);
% computing L
aL = [1 2*(2*pi*6200*.68) (2*pi*6200)^2];
bL = [0 0 (2*pi*6200)^2];
L = tf(bL,aL);
[mL,pL,fL] = bode(L,f*2*pi);

% P
P = tf(.03);
[mP,pP,fP] = bode(P,f*2*pi);
%% plotting DC - 10kHz response of P,L,G
% G
subplot(2,3,1);
hold on
% for i = 1:L
%     a = plot(f,Gmag(i,:));
%     a.Color = [.8 .8 .8];
%     a.LineWidth = .5;
% end
p = plot(f,Gmeanmag - mean(Gmeanmag(20:800)));
p.Color = [1,.8,.8];
p.LineWidth = .8;
p = plot(f,20*log10(squeeze(mG)),'--');
p.Color = [0 0 0];
p.LineWidth = .8;
p = plot(f(r1:r2),Gmeanmag(r1:r2) - mean(mean(Gmeanmag(20:800))));
p.Color = [1,0,0];
p.LineWidth = 1.2;
set(gca, 'XScale', 'log');
axis([8 6000 -11 11])
yticks([-10 0 10])
xticks([10 100 1000]);
set(gca,'TickLabelInterpreter','latex','FontSize',10);
text(20, 5,'$|$G$|$','Interpreter','latex','FontSize',10);
ylabel('magnitude (dB)','Interpreter','latex','FontSize',10);
subplot(2,3,4);
hold on
% for i = 1:L
%     a = plot(f,Gang(i,:));
%     a.Color = [.8 .8 .8];
%     a.LineWidth = .5;
% end
p = plot(f,Gang(12,:));
p.Color = [1,.8,.8];
p.LineWidth = .8;
p = plot(f,squeeze(pG),'--');
p.Color = [0 0 0];
p.LineWidth = .8;
p = plot(f(r1:r2),Gang(12,r1:r2));
p.Color = [1,0,0];
p.LineWidth = 1.2;
set(gca, 'XScale', 'log');
axis([8 6000 -100 100])
yticks([-90 0 90]);
xticks([10 100 1000]);
set(gca,'TickLabelInterpreter','latex','FontSize',10);
ylabel('phase (deg)','Interpreter','latex','FontSize',10);
text(20, 45,'$\angle$G','Interpreter','latex','FontSize',10);
% L
subplot(2,3,2);
hold on
p = plot(freqsL,20*log10(Lmag));
p.Color = [1,.8,.8];
p.LineWidth = .8;
p = plot(f,20*log10(squeeze(mL)),'--');
p.Color = [0 0 0];
p.LineWidth = .8;
p = plot(interp1(freqsL,freqsL,10:1000),20*log10(interp1(freqsL,Lmag,10:1000)));
p.Color = [1,0,0];
p.LineWidth = 1.2;
axis([8 6000 -11 11])
yticks([-10 0 10]);
xticks([10 100 1000]);
set(gca, 'XScale', 'log');
set(gca,'TickLabelInterpreter','latex','FontSize',10);
text(20, 5,'$|$L$|$','Interpreter','latex','FontSize',10);
subplot(2,3,5);
hold on
p = plot(freqsL,Lang);
p.Color = [1,.8,.8];
p.LineWidth = .8;
p = plot(f,squeeze(pL),'--');
p.Color = [0 0 0];
p.LineWidth = .8;
p = plot(interp1(freqsL,freqsL,10:1000),interp1(freqsL,Lang,10:1000));
p.Color = [1,0,0];
p.LineWidth = 1.2;
axis([8 6000 -100 100])
yticks([-90 0 90]);
xticks([10 100 1000]);
set(gca, 'XScale', 'log');
set(gca,'TickLabelInterpreter','latex','FontSize',10);
text(20, 45,'$\angle$L','Interpreter','latex','FontSize',10);
xlabel('frequency (Hz)','Interpreter','latex','FontSize',10);
% P
subplot(2,3,3);
hold on
p = plot(f,20*log10(squeeze(mP)),'--');
p.Color = [0 0 0];
p.LineWidth = .8;
p = plot(freqsP,20*log10(m));
p.Color = [1 0 0];
p.LineWidth = 1;
set(gca, 'XScale', 'log');
axis([8 6000 -40 -20])
yticks([-40 -30 -20]);
xticks([10 100 1000]);
set(gca,'TickLabelInterpreter','latex','FontSize',10);
text(20, -25,'$|$P$|$','Interpreter','latex','FontSize',10);
subplot(2,3,6);
hold on
p = plot(f,squeeze(mL),'--');
p.Color = [0 0 0];
p.LineWidth = .8;
p = plot(freqsP,pl);
p.Color = [1 0 0];
p.LineWidth = 1;
set(gca, 'XScale', 'log');
axis([8 6000 -100 100])
yticks([-90 0 90]);
xticks([10 100 1000]);
set(gca,'TickLabelInterpreter','latex','FontSize',10);
text(20, 45,'$\angle$P','Interpreter','latex','FontSize',10);
set(gcf,'color','w');
%% estimating transfer functions of P,G,L

P = 1;

aG = [1 2*(2*pi*4400*.6) (2*pi*4400)^2];
bG = [0 0 (2*pi*4400)^2];
G = tf(bG,aG);

aL = [1 2*(2*pi*6200*.68) (2*pi*6200)^2];
bL = [0 0 (2*pi*6200)^2];
L = tf(bL,aL);

[bb,aa] = butter(2,35*2/10000,'high');
[bb,aa] = ztos(bb,aa,10000);
H = tf(bb,aa);

%% calculating C fit
Co = 1/(P*(1 - (G*L))); % reference tracking
J = Co*P*L*G;
D = tf([.01*2*pi],[1 1]);

Oo = Co*P/((Co*P*L*G)+1);
[k2,l2,p2] = bode(Co,logspace(-2,5,100000));
[ko2,lo2,po2] = bode(Oo,logspace(-2,5,100000));
% 
% a21 = [0 2*(2*pi*4400*.68) 0];
% b21 = [1 2*(2*pi*4400*.68) ((2*pi*4400)^2)]*18.5/4;
% 
% a22 = [1 2*(2*pi*7100*.66) ((2*pi*7100)^2)];
% b22 = [0 0 (2*pi*7100)^2];

a21 = [0 2*(2*pi*3500*.8) 0];
b21 = [1 2*(2*pi*3500*.8) ((2*pi*3500)^2)]*10.5*.05;

a22 = [1 2*(2*pi*6000*.75) ((2*pi*6000)^2)];
b22 = [0 0 (2*pi*6000)^2];

C = tf(b21,a21)*tf(b22,a22);
%C = tf(bc,ac);
[b2,a2] = tfdata(C);
b2 = cell2mat(b2);
a2 = cell2mat(a2);
[b2d,a2d] = stoz(b2,a2,10000);

O = C*P/((C*P*L*G)+1);
S = 1 - O;
[bo,ao] = tfdata(O);
bo = cell2mat(bo);
ao = cell2mat(ao);
[bdo,ado] = stoz(bo,ao,10000);

[k,l,p] = bode(C,logspace(-2,5,100000));
[kd,d,~] = freqz(b2d,a2d,5000,10000);
ld = unwrap(angle(kd))*180/pi;
kd = 20*log10(abs(kd));

[ko,lo,po] = bode(O,logspace(-2,5,100000));
[kdo,do,~] = freqz(bdo,ado,5000,10000);
ldo = unwrap(angle(kdo))*180/pi;
kdo = 20*log10(abs(kdo));
%% plotting C fit
tsize = 10;
a = figure;
%a.Position(3) = 500;
%a.Position(4) = 250;
subplot(2,2,1);
hold on
pl = plot(squeeze(p2)/2/pi,20*log10(squeeze(k2)));
pl.Color = [0 0 0];
pl.LineWidth = .7;
pl = plot(d,kd);
pl.Color = [1 0 0];
pl.LineWidth = .7;
set(gca, 'XScale', 'log');
ylabel('magnitude (dB)','Interpreter','latex','FontSize',tsize);
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
axis([8 5000 0 95])
xticks([10 100 1000]);
yticks([0 20 40 60 80])
pl = legend({'ideal C','designed C'},'Interpreter','latex','FontSize',tsize);
pl.Position(2) = pl.Position(2)+.05;
pl.Box = 'off';
subplot(2,2,3);
hold on
pl = plot(squeeze(p2)/2/pi,squeeze(l2));
pl.Color = [0 0 0];
pl.LineWidth = .7;
pl = plot(d,ld-360);
pl.Color = [1 0 0];
pl.LineWidth = .7;
set(gca, 'XScale', 'log');
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('phase (deg)','Interpreter','latex','FontSize',tsize);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
axis([8 5000 -180 0])
yticks([-180 -120 -60 0])
xticks([10 100 1000]);

%% plotting O fit
tsize = 10;
% a = figure;
% a.Position(3) = 400;
% a.Position(4) = 300;
subplot(2,2,2);
hold on
%a = plot(10:1000,20*log10(interp1(squeeze(p2)/2/pi,squeeze(ko2),10:1000)));
%a.Color = [0 0 0];
%a.LineWidth = 1;
%a = plot(20:1000,20*log10(interp1(squeeze(p2)/2/pi,squeeze(ko),20:1000)));
%a = plot(d(11:1000),kdo(11:1000));
%a.Color = [1 0 0];
%a.LineWidth = 1;

pl = plot(squeeze(p2)/2/pi,20*log10(squeeze(ko2)));
pl.Color = [0 0 0];
pl.LineWidth = .7;

%a = plot(10:1000,20*log10(interp1(squeeze(p2)/2/pi,squeeze(ko2),10:1000)));
%a.Color = [0 0 0];
%a.LineWidth = 1;
%a = plot(squeeze(p2)/2/pi,20*log10(squeeze(ko)));
pl = plot(d,kdo);
pl.Color = [1 0 0];
pl.LineWidth = .7;
%a = plot(20:1000,20*log10(interp1(squeeze(p2)/2/pi,squeeze(ko),20:1000)));
%a = plot(d(11:1000),kdo(11:1000));
%a.Color = [1 0 0];
%a.LineWidth = 1;
pl = legend({'ideal T','designed T'},'Interpreter','latex','FontSize',tsize);
pl.Location = 'northeast';
pl.Box = 'off';
pl.Position(2) = pl.Position(2)+.05;

set(gca, 'XScale', 'log');
ylabel('magnitude (dB)','Interpreter','latex','FontSize',tsize);
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
axis([8 5000 -13 13])
yticks([-12 -6 0 6 12])
% xticks([10,20,50,100,200,500,1000,3000])
xticks([10 100 1000]);
subplot(2,2,4);
hold on

pl = plot(squeeze(p2)/2/pi,squeeze(lo2));
pl.Color = [0 0 0];
pl.LineWidth = .7;
%a = plot(10:1000,interp1(squeeze(p2)/2/pi,squeeze(lo2),10:1000))
%a.Color = [0 0 0];
%a.LineWidth = 1;
%a = plot(squeeze(p2)/2/pi,squeeze(lo)+360);
pl = plot(d,ldo-360);
pl.Color = [1 0 0];
pl.LineWidth = .7;
%a = plot(20:1000,interp1(squeeze(p2)/2/pi,squeeze(lo)+360,20:1000));
%a = plot(d(11:1000),ldo(11:1000)-360);
%a.Color = [1 0 0];
%a.LineWidth = 1;

set(gca, 'XScale', 'log');
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
ylabel('phase (deg)','Interpreter','latex','FontSize',tsize);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
axis([8 5000 -100 45])
yticks([-90 -45 0 45])
% xticks([10,20,50,100,200,500,1000,3000])
xticks([10 100 1000]);
%% creating force segmnets to cose the loop over
directory = 'C:\Users\atrox\Desktop\Work\Research\projects\z Finished\High-bandwidth tribometry as means of capturing natural texture\NU tribometry\processed data\';
dirr = dir(directory);
dirr = dirr(3:end);
speeds = {'80','120','180','270'};
[bl,al] = butter(1,10*2/125000,'high');
[bh,ah] = butter(1,1000*2/125000,'low');
SIG = {};
T = {};
for t = 1:7
    for s = 1:1
        load(strcat(directory,dirr(t).name,'\Finger\',speeds{s},'\Features\Features.mat'));
        sig = [];
%         T{t,s} = length(LATT{1})/125000;
%         for k = 1:length(LATT)
%             if (k==1)
%                 sig = filtfilt(bl,al,LATT{k}.');
%                 sig = filtfilt(bh,ah,sig);
%             else
%                 sigtemp = filtfilt(bl,al,LATT{k}.');
%                 sigtemp = filtfilt(bh,ah,sigtemp);
%                 sig(end-6249:end) = (sigtemp(1:6250).*(1:6250) + sig(end-6249:end).*flip(1:6250))/6250;
%                 sig = [sig,sigtemp(6251:end)];
%             end
%         end       
%         sigtemp = sig;
%         while length(sig)<(10*1250000)
%             sig(end-6249:end) = (sigtemp(1:6250).*(1:6250) + sig(end-6249:end).*flip(1:6250))/6250;
%             sig = [sig,sigtemp(6251:end)];
%         end
%         sig = filtfilt(bh,ah,sig);
%         sig = resample(sig,60000,125000);
%         SIG{t,s} = sig(1:10*60000);
        k = 1;
        while length(LATT{k})<62500
            k = k + 1;
        end
        sig = filtfilt(bl,al,LATT{k}.');
        
        sigtemp = sig;
        while length(sig)<(20*125000)
            sig(end-6249:end) = (sigtemp(1:6250).*(1:6250) + sig(end-6249:end).*flip(1:6250))/6250;
            sig = [sig,sigtemp(6251:end)];
        end
        sig = filtfilt(bh,ah,sig);
        sig = resample(sig,60000,125000);
        SIG{t,s} = sig(1:20*60000);        
    end
    t
end
save('C:\Users\atrox\Desktop\Work\Research\projects\z Finished\High-bandwidth tribometry as means of capturing natural texture\NU tribometry\combined recordings\data.mat','SIG');

%% generating white noise sig
[bl,al] = butter(1,10*2/60000,'high');
[bh,ah] = butter(1,1000*2/60000,'low');
sig = randn(1,60000*(dur+2));
sig = filtfilt(bl,al,sig);
sig = filtfilt(bh,ah,sig);
sig = sig(60001:end-60000);
sig = sig/std(sig)/10*2;

SIG = {};
for t = 1:10
    SIG{t,1} = sig*t;
end

%% playing back texture signal
tex = randperm(7);
dur = 20;
sr = 60000;
t = linspace(0,dur,dur*sr);
st = daq.createSession('ni');
st.Rate = sr;
st.DurationInSeconds = dur;
ch1 = addAnalogInputChannel(st,'Dev2','ai1','Voltage');
ch1.InputType = 'SingleEnded'; %% lateral force
ch2 = addAnalogInputChannel(st,'Dev2','ai10','Voltage');
ch2.InputType = 'SingleEnded'; %% normal force
ch3 = addAnalogInputChannel(st,'Dev2','ai2','Voltage');
ch3.InputType = 'SingleEnded'; %% reference force
ch4 = addAnalogInputChannel(st,'Dev2','ai9','Voltage');
ch4.InputType = 'SingleEnded'; %% current
ch5 = addAnalogOutputChannel(st,'Dev2','ao0','Voltage');
DATAP = {};
for te = 1:7
    for sp = 1     
        sig = SIG{tex(te),sp}(1:dur*sr)*100;
        fprintf('working on textue%3.0f\n working on speed%3.0f\n',te,sp),
        sig(end) = 0;
        queueOutputData(st,sig');
        x = input('press enter \n');
        out = startForeground(st);
        DATAP{tex(te),1} = out;
    end
end
save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\closed_loop_results\textures\data11.mat','DATAP','sr');

%% processing texture playback under open and closed loop conditions
dirs = {};
dirs{1} = 'C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\closed_loop_results\textures\data5.mat';
dirs{2} = 'C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\closed_loop_results\textures\data9.mat';
sr = 60000;
[bl1,al1] = butter(3,10/sr,'low');          % used for finding changes in direction
[bl2,al2] = butter(1,10*2/sr,'high');       % used for removing low frequency kinematic information
[bl3,al3] = butter(2,1000*2/sr,'low');      % used for removing high frequency noise
[bn,an] = butter(2,2*50/sr,'low');            % used for filtering out high frequency normal force components
lat_const = .25;
nor_const = -.1526;
resultsmat = {};
fftmat = {};
shiftmat = {};
kvals = [];
ct = {};
S = {};
E = {};
sh = [];
for p = 1:2
    load(dirs{p});
    for te = 1:7
        for sp = 1:1
            lat = DATAP{te,sp}(:,1)*lat_const;
            nor = DATAP{te,sp}(:,2)*nor_const;
            sig = DATAP{te,sp}(:,3);
            cur = DATAP{te,sp}(:,4);
            lat1 = detrend(filtfilt(bl1,al1,detrend(lat)));
            ii = 1;
            while abs(lat1(ii))<.1
                ii = ii+10;
            end
            sh(te,sp,p) = ii;
            lat = lat(ii:end);
            latd = detrend(lat);
            sig = sig(ii:end);
            lat1 = lat1(ii:end);
            cur = cur(ii:end);
            norf = filtfilt(bn,an,nor);
            norf = (norf - mean(norf(sr/8:sr/2)))*nor_const;
            norf = norf(ii:end);
            lat1 = lat1>0;
            lat1 = abs(diff(lat1)).*(1:length(lat1)-1)';
            lat1 = lat1(lat1~=0);
            
            latfs = filtfilt(bl2,al2,lat);
            latf = filtfilt(bl3,al3,latfs);

            sigf = filtfilt(bl2,al2,sig);
            sigf = filtfilt(bl3,al3,sigf);
            
            latf = latf.*sign(latd);
            latfs = latfs.*sign(latd);
            lat1 = lat1(2:end);
            k = 1;
            resultsmat{te,sp,p,1} = [];
            resultsmat{te,sp,p,8} = [];
            resultsmat{te,sp,p,9} = [];
            for ii = 1:length(lat1)-1
                s = round((lat1(ii+1)-lat1(ii))*1/4)+lat1(ii);
                e = round((lat1(ii+1)-lat1(ii))*3/4)+lat1(ii);
                lattf = latf(s:e);
                lattfs = latfs(s:e);
                norr = norf(s:e);
                sigg = sigf(s:e);
                curr = cur(s:e);
                
                if (max(abs(curr))<4 && length(sigg)>3000)
                    S{te,k,p} = s;
                    E{te,k,p} = e;
                    cormat = [];
                    for c = 1:100
                        cormat(c) = corr(sigg(1:end-c)/100,lattf(c+1:end));
                    end
                    
                    [~,c] = max(cormat);
                    shiftmat{te,sp,p,k} = c;
                    ss = sigg(1:end-c)/100;
                    ll = lattf(c+1:end);
                    sres = sum((ss - ll).^2);
                    stot = sum((ss- mean(ss)).^2);
                    
                    resultsmat{te,sp,p,1} = [1 - (sres/stot),resultsmat{te,sp,p,1}];
                    resultsmat{te,sp,p,9} = [resultsmat{te,sp,p,9},corr(ss,ll)];
                    resultsmat{te,sp,p,8} = [mean(abs((sigg(1:end-c)/100)-lattf(c+1:end)))/mean(abs((sigg(1:end-c)/100))),resultsmat{te,sp,p,8}];
                    if k == 6
                        resultsmat{te,sp,p,2} = lattf;
                        resultsmat{te,sp,p,3} = sigg;
                        resultsmat{te,sp,p,4} = curr;
                        resultsmat{te,sp,p,5} = lattfs;
                        resultsmat{te,sp,p,6} = norr;
                        resultsmat{te,sp,p,7} = c;
                        %fftmat{te,sp,p,1} = fftromanw(sigg/100,12000,6000,hann(12000));
                        %fftmat{te,sp,p,2} = fftromanw(lattf,12000,6000,hann(12000));
                    end
                    ct{te,sp,p} = k;
                    k = k + 1;                    
                end
            end
        end
        te
    end
end
%% creating noise signals to close the loop over
sr = 60000;
nsignals = 5;
[bl,al] = butter(1,10*2/sr,'high');
[bh,ah] = butter(1,1000*2/sr,'low');
dur = 10;
SIG = {};
for n = 1:nsignals
    sig = randn(1,sr*dur);
    sigf = detrend(filtfilt(bh,ah,sig));
    sigf = sigf/std(sigf);
    sigf = n/nsignals*sigf;
    SIG{n,1} = sigf;        
    sig = randn(1,sr*dur);
    sigf = detrend(pink_filter(10,10000,sr,sig));
    sigf = sigf/std(sigf);
    sigf = n/nsignals*sigf;
    SIG{n,2} = sigf;
end
save('C:\Users\atrox\Desktop\Work\Research\projects\z Finished\High-bandwidth tribometry as means of capturing natural texture\NU tribometry\combined recordings\noise.mat','SIG');

%% playback of nosie signals
tex = randperm(5);
dur = 10;
sr = 60000;
st = daq.createSession('ni');
st.Rate = sr;
st.DurationInSeconds = dur;
ch1 = addAnalogInputChannel(st,'Dev2','ai1','Voltage');
ch1.InputType = 'SingleEnded'; %% lateral force
ch2 = addAnalogInputChannel(st,'Dev2','ai10','Voltage');
ch2.InputType = 'SingleEnded'; %% normal force
ch3 = addAnalogInputChannel(st,'Dev2','ai2','Voltage');
ch3.InputType = 'SingleEnded'; %% reference force
ch4 = addAnalogInputChannel(st,'Dev2','ai9','Voltage');
ch4.InputType = 'SingleEnded'; %% current
ch5 = addAnalogOutputChannel(st,'Dev2','ao0','Voltage');
DATAP = {};
for te = 1:5
    for sp = 1     
        sig = SIG{tex(te),1};
        fprintf('working on textue%3.0f\n working on speed%3.0f\n',te,1),
        sig(end) = 0;
        queueOutputData(st,sig');
        x = input('press enter \n');
        out = startForeground(st);
        DATAP{tex(te),1} = out;
        sig = SIG{tex(te),2};
        fprintf('working on textue%3.0f\n working on speed%3.0f\n',te,2),
        sig(end) = 0;
        queueOutputData(st,sig');
        x = input('press enter \n');
        out = startForeground(st);
        DATAP{tex(te),2} = out;
    end
end
save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\closed_loop_results\textures\datan4.mat','DATAP','sr');

%% processing texture playback of noise signal
dirs = {};
dirs = 'C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\closed_loop_results\textures\datan2.mat';
sr = 60000;
[bl1,al1] = butter(1,10/sr,'low');          % used for finding changes in direction
[bl2,al2] = butter(1,10*2/sr,'high');       % used for removing low frequency kinematic information
[bl3,al3] = butter(1,1000*2/sr,'low');      % used for removing high frequency noise
[bn,an] = butter(2,2*50/sr,'low');            % used for filtering out high frequency normal force components
lat_const = .25;
nor_const = -.1526;
resultsmat = {};
fftmat = {};
shiftmat = {};
kvals = [];
ct = {};
S = {};
E = {};
sh = [];
load(dirs);
for p = 1:2
    for te = 1:5
        for sp = 1:1
            lat = DATAP{te,p}(:,1)*lat_const;
            nor = DATAP{te,p}(:,2)*nor_const;
            sig = DATAP{te,p}(:,3);
            cur = DATAP{te,p}(:,4);
            lat1 = detrend(filtfilt(bl1,al1,detrend(lat)));
            ii = 1;
            while abs(lat1(ii))<.1
                ii = ii+10;
            end
            sh(te,p) = ii;
            lat = lat(ii:end);
            latd = detrend(lat);
            sig = sig(ii:end);
            lat1 = lat1(ii:end);
            cur = cur(ii:end);
            norf = filtfilt(bn,an,nor);
            norf = (norf - mean(norf(sr/8:sr/2)))*nor_const;
            norf = norf(ii:end);
            lat1 = lat1>0;
            lat1 = abs(diff(lat1)).*(1:length(lat1)-1)';
            lat1 = lat1(lat1~=0);
            
            latfs = filtfilt(bl2,al2,lat);
            latf = filtfilt(bl3,al3,latfs);

            sigf = filtfilt(bl2,al2,sig);
            sigf = filtfilt(bl3,al3,sigf);
            
            latf = latf.*sign(latd);
            latfs = latfs.*sign(latd);
            lat1 = lat1(2:end);
            k = 1;
            resultsmat{te,p,1} = [];
            resultsmat{te,p,8} = [];
            resultsmat{te,p,9} = [];
            for ii = 1:length(lat1)-1
                s = round((lat1(ii+1)-lat1(ii))*1/4)+lat1(ii);
                e = round((lat1(ii+1)-lat1(ii))*3/4)+lat1(ii);
                lattf = latf(s:e);
                lattfs = latfs(s:e);
                norr = norf(s:e);
                sigg = sigf(s:e);
                curr = cur(s:e);
                
                if (max(abs(curr))<3 && length(sigg)>3000)
                    S{te,k,p} = s;
                    E{te,k,p} = e;
                    cormat = [];
                    for c = 1:100
                        cormat(c) = corr(sigg(1:end-c)/100,lattf(c+1:end));
                    end
                    
                    [~,c] = max(cormat);
                    shiftmat{te,p,k} = c;
                    ss = sigg(1:end-c)/100;
                    ll = lattf(c+1:end);
                    sres = sum((ss - ll).^2);
                    stot = sum((ss- mean(ss)).^2);
                    
                    resultsmat{te,p,1} = [1 - (sres/stot),resultsmat{te,p,1}];
                    resultsmat{te,p,9} = [resultsmat{te,p,9},corr(ss,ll)];
                    resultsmat{te,p,8} = [mean(abs((sigg(1:end-c)/100)-lattf(c+1:end)))/mean(abs((sigg(1:end-c)/100))),resultsmat{te,p,8}];
                    if k == 1
                        resultsmat{te,p,2} = lattf;
                        resultsmat{te,p,3} = sigg;
                        resultsmat{te,p,4} = curr;
                        resultsmat{te,p,5} = lattfs;
                        resultsmat{te,p,6} = norr;
                        resultsmat{te,p,7} = c;
                    end
                    ct{te,p} = k;
                    k = k + 1;                    
                end
            end
        end
        te
    end
end
%% extracting points of interseaction
CPTS = {};
for te = 1:9
    for sp = 1:1
        k = 1;
        CPTS{te,sp} = [];
        s1 = []; e1 = []; s2 = []; e2 = [];
        for c = 1:ct{te,sp,1}
            s1(c) = S{te,c,1}+sh(te,sp,1);
            e1(c) = E{te,c,1}+sh(te,sp,1);
        end
        for c = 1:ct{te,sp,2}
            s2(c) = S{te,c,2}+sh(te,sp,2);
            e2(c) = E{te,c,2}+sh(te,sp,2);
        end
        for s = 1:length(s1)
            for ss = 1:length(s2)
                if (e1(s) - s2(ss)) > 3000 && (e1(s) - s2(ss)) < 10000
                    CPTS{te,sp}(k,1:4) = [s2(ss)-sh(te,sp,1),e1(s)-sh(te,sp,1),s2(ss)-sh(te,sp,2),e1(s)-sh(te,sp,2)];
                    k = k + 1;
                end
            end
        end
    end
end
%%
texts = [6,7,5,1,3,2,8,9];
textnames = {'EV','FL','HT','HT','DM','MS','SW','WN','PN'};
for te = 1:8
    for sp = 1:1
        subplot(1,8,te)
        hold on
        pl = plot([-100 100],[-100 100],':');
        pl.Color = [.6 .6 .6];
        pl.LineWidth = 1;
        pl = plot(resultsmat{texts(te),sp,p,3}(1:end-shiftmat{texts(te),sp,p,1})*10,resultsmat{texts(te),sp,p,2}(shiftmat{texts(te),sp,p,1}+1:end)*1000);
        if p == 2
            pl.Color = [.4 .4 .4 .5];
        else
            pl.Color = [1 0 0 .5];
        end
        
        if te == 1
            ylabel('rendered friction (mN)','Interpreter','latex','FontSize',10);
            xlabel('reference friction (mN)','Interpreter','latex','FontSize',10);
        end
        %title(textnames{texts(te)},'Interpreter','latex','FontSize',5);
        set(gca,'TickLabelInterpreter','latex','FontSize',7);
    if (te == 1 || te == 2)
        axis([-12.5 12.5 -12.5 12.5]);
        %yticks([-10 0 10]);
        yticks([]);
        xticks([-10 0 10]);        
    elseif (te == 5 || te == 4 || te == 3)
        axis([-17.5 17.5 -17.5 17.5]);
        %yticks([-15 0 15]);    
        yticks([])
        xticks([-15 0 15]);  
    elseif (te == 7 || te == 8)
        axis([-32.5 32.5 -32.5 32.5]);
        %yticks([-30 0 30]);
        yticks([])
        xticks([-30 0 30]); 
    else
        axis([-42.5 42.5 -42.5 42.5]);
        %yticks([-40 0 40]);
        yticks([])
        xticks([-40 0 40]);
    end
    end
end
set(gcf,'color','w');
set(gcf,'renderer','Painters')
%% plotting tracking comparison
pl = figure;
pl.Position(3) = 400;
pl.Position(4) = 200;
te = 3;
hold on
for p = 1:2
    load(dirs{p});
    lat = DATAP{te,sp}(:,1)*lat_const;
    nor = DATAP{te,sp}(:,2);
    sig = DATAP{te,sp}(:,3);
    cur = DATAP{te,sp}(:,4);
    lat1 = detrend(filtfilt(bl1,al1,detrend(lat)));
    ii = 1;
    while abs(lat(ii))<.1
        ii = ii+10;
    end
    lat = lat(ii:end);
    latd = detrend(lat);
    sig = sig(ii:end);
    lat1 = lat1(ii:end);
    cur = cur(ii:end);
    norf = filtfilt(bn,an,nor);
    norf = (norf - mean(norf(sr/8:sr/2)))*nor_const;
    norf = norf(ii:end);
    lat1 = lat1>0;
    lat1 = abs(diff(lat1)).*(1:length(lat1)-1)';
    lat1 = lat1(lat1~=0);
    
    latfs = filtfilt(bl2,al2,lat);
    latf = filtfilt(bl3,al3,latfs);
    sigf = sig;
    sigf = filtfilt(bl2,al2,sig);
    sigf = filtfilt(bl3,al3,sigf);
    
    latf = latf.*sign(latd);
    if p == 1
         pl = plot(sigf(CPTS{te,sp}(1,1):CPTS{te,sp}(1,2))/100);
         pl.Color = [0 0 0];
         pl.LineWidth = 1;
         pl = plot(latf(CPTS{te,sp}(1,1)+22:CPTS{te,sp}(1,2)+22));
         pl.Color = [1 0 0];
         pl.LineWidth = 1;        
    else
         pl = plot(latf(CPTS{te,sp}(1,3):CPTS{te,sp}(1,4)));
         pl.Color = [0 0 1];
         pl.LineWidth = 1;  
         pl = plot(sigf(CPTS{te,sp}(1,3):CPTS{te,sp}(1,4))/100);
         pl.Color = [0 1 0];
         pl.LineWidth = 1;    
    end
end
%% plotting R squared for nat texture

texts = [7,6,1,3,5,2];
textnames = {'EV','FL','HT','HT','DM','MS','SW'};
%subplot(1,2,1]]);
hold on
for i = 1:6
    st = std(resultsmat{texts(i),1,p,1});
    mn = mean(resultsmat{texts(i),1,p,1});
    pl = plot([i,i]+.05,[mn,mn]+[-st,st]);
    if p == 1
        pl.Color = [1 0 0];
        pl.LineWidth = 1.5;
        pl = plot(i+.05,mn,'.');
        pl.MarkerSize = 16;
        pl.Color = [1 0 0];
    else
        pl.Color = [0 0 0];
        pl.LineWidth = 1.5;
        pl = plot(i+.05,mn,'.');
        pl.MarkerSize = 16;
        pl.Color = [0 0 0];
    end
end
xticklabels(textnames(texts));
xticks([1,2,3,4,5,6,7,8]);
set(gca,'TickLabelInterpreter','latex','FontSize',10);
axis([.5 8.5 -.5 1]);
yticks([0 .2 .4 .6 .8 1.0])
set(gcf,'color','w');
ylabel('R$^{2}$','Interpreter','latex','FontSize',11);
xlabel('texture','Interpreter','latex','FontSize',11);
pl = plot([-1 10],[.9 .9],':');
pl.Color = [0 0 0];
pl.LineWidth = .5;
% 
% subplot(1,2,2);
% hold on
% for i = 1:8
%     st = std(resultsmat{texts(i),1,9});
%     mn = mean(resultsmat{texts(i),1,9});
%     pl = plot([i,i]+.05,[mn,mn]+[-st,st]);
%     pl.Color = [.3 .3 .3];
%     pl = plot(i+.05,mn,'.');
%     pl.MarkerSize = 11;
%     pl.Color = [0 0 0];
% end
% xticklabels(textnames(texts));
% xticks([1,2,3,4,5,6,7,8]);
% set(gca,'TickLabelInterpreter','latex','FontSize',8);
% axis([.5 8.5 0 1]);
% yticks([0 .2 .4 .6 .8 .9 1.0])
% set(gcf,'color','w');
% ylabel('Pearson r','Interpreter','latex','FontSize',10);
% xlabel('texture','Interpreter','latex','FontSize',10);
% 
% pl = plot([-1 10],[.95 .95],':');
% pl.Color = [0 0 0];
% pl.LineWidth = .5;

%% plotting R squared for noise

hold on
for i = 1:5
    st = std(resultsmat{i,p,1});
    mn = mean(resultsmat{i,p,1});
    pl = plot([i,i]+.05,[mn,mn]+[-st,st]);
    if p == 1
        pl.Color = [0 0 0];
        pl.LineWidth = 1.5;
        pl = plot(i+.05,mn,'.');
        pl.MarkerSize = 16;
        pl.Color = [0 0 0];
    else
        pl.Color = [0 0 0];
        pl.LineWidth = 1.5;
        pl = plot(i+.05,mn,'.');
        pl.MarkerSize = 16;
        pl.Color = [0 0 0];
    end
end
xticks([1,2,3,4,5,6,7,8]);
set(gca,'TickLabelInterpreter','latex','FontSize',10);
axis([.5 8.5 -.5 1]);
yticks([0 .2 .4 .6 .8 1.0])
set(gcf,'color','w');
ylabel('R$^{2}$','Interpreter','latex','FontSize',11);
xlabel('texture','Interpreter','latex','FontSize',11);
pl = plot([-1 10],[.9 .9],':');
pl.Color = [0 0 0];
pl.LineWidth = .5;
% 
% subplot(1,2,2);
% hold on
% for i = 1:8
%     st = std(resultsmat{texts(i),1,9});
%     mn = mean(resultsmat{texts(i),1,9});
%     pl = plot([i,i]+.05,[mn,mn]+[-st,st]);
%     pl.Color = [.3 .3 .3];
%     pl = plot(i+.05,mn,'.');
%     pl.MarkerSize = 11;
%     pl.Color = [0 0 0];
% end
% xticklabels(textnames(texts));
% xticks([1,2,3,4,5,6,7,8]);
% set(gca,'TickLabelInterpreter','latex','FontSize',8);
% axis([.5 8.5 0 1]);
% yticks([0 .2 .4 .6 .8 .9 1.0])
% set(gcf,'color','w');
% ylabel('Pearson r','Interpreter','latex','FontSize',10);
% xlabel('texture','Interpreter','latex','FontSize',10);
% 
% pl = plot([-1 10],[.95 .95],':');
% pl.Color = [0 0 0];
% pl.LineWidth = .5;


%% plotting temp and freq signals
pl = figure;
pl.Position(3) = 800;
pl.Position(4) = 300;
t = (0:15000-1)/60000*1000;
f = (0:12000-1)/12000*60000;
texts = [6,3,5,2,1,7,8,9];
textnames = {'EV','DM','HT','FL','MS','SW','WN','PN'};
f = sr*(1:6000)/12000;
for te = 1:length(texts)
    subplot(2,4,te)
    hold on
    pl = plot(t-80,resultsmat{texts(te),sp,1,3}(1:15000)/100*1000);
    pl.Color = [.4 .4 .4];
    pl.LineWidth = .8;
    
    pl = plot(t(1:end-resultsmat{texts(te),sp,1,7}+1)-80,resultsmat{texts(te),sp,1,2}(resultsmat{texts(te),sp,1,7}:15000)*1000);
    pl.Color = [1 0 0 1];
    pl.LineWidth = .6;
    if (te == 5 || te == 1)
        axis([0 50 -7.5 6.5]);
        yticks([-5 0 5]);
    elseif (te == 4)
        axis([0 50 -32.5 32.5]);
        yticks([-30 -15 0 15 30]);
    elseif (te == 6 || te == 2)
        axis([0 50 -12.5 12.5]);
        yticks([-10 0 10]);    
    elseif (te == 3)
        axis([0 50 -17.5 17.5]);
        yticks([-15 0 15]);            
    else
        axis([0 50 -32.5 32.5]);
        yticks([-30 -15 0 15 30]);
    end
    xticks([0 25 50])
    set(gca,'TickLabelInterpreter','latex','FontSize',9);
    if te==1
        ylabel('force (mN)','Interpreter','latex','FontSize',10);
    end
    title(textnames{te},'Interpreter','latex','FontSize',8);
 
%     subplot(2,8,te+8)
%     hold on
%     pl = plot(f,fftmat{texts(te),sp,1}/max(fftmat{texts(te),sp,1}));
%     pl.Color = [.7 .7 .7];
%     pl.LineWidth = .8;
%     pl = plot(f,fftmat{texts(te),sp,2}/max(fftmat{texts(te),sp,1}));
%     pl.Color = [1 0 0 1];
%     pl.LineWidth = .8;   
%     set(gca, 'XScale', 'log');
%     %set(gca, 'YScale', 'log');
%     axis([5 600 0 1.05]);
%     xticks([10 100]);
end
xlabel('time (ms)','Interpreter','latex','FontSize',10);
set(gcf,'color','w');
set(gcf,'renderer','Painters')


%%
plot(resultsmat{te,1,1,2},resultsmat{te,1,1,3}/100);
%% plotting a single swipe
lat_const = .25;
nor_const = -.1526;
lat = lat_const * out(:,1);
cur = (out(:,4)+5)/2;
sig = out(:,3)/100;
nor = nor_const * out(:,2);
[b,a] = butter(1,2/sr*1000,'low');
[bl1,al1] = butter(1,10*2/sr,'high'); 
[bn,an] = butter(2,2/sr*50,'low');
latf = lat;
latf = filtfilt(b,a,lat);
norf = filtfilt(bn,an,nor);
norf = norf - mean(norf(10000:30000));
latf = latf - mean(latf(10000:30000));

x1 = 320000:440000;
norf = norf(x1);
latf = latf(x1);
%curf = filtfilt(b,a,cur);
curf = cur(x1);
sig = sig(x1);
x = (x1 - x1(1))/sr;

%x1 = (.06113:.07113)*sr;
set(gca,'TickLabelInterpreter','latex','FontSize',10);
subplot(3,1,1);
hold on
s = 24000;
e = 28500;

pl = plot((x(s:e)-x(s))*8 + x(7500),8*(sig(s:e) + mean(latf(s:e) - sig(s:e)) - .18));
pl.Color = [0 0 1];
pl = plot(x,latf);
pl.Color = [0 0 0];
pl = plot((x(s:e)-x(s))*8 + x(7500),8*(latf(s:e)-.18));
pl.Color = [0 0 0];
axis([0 (length(x1)/sr) -.35 .35]);
ylabel('friction force (N)','Interpreter','latex','FontSize',10);
set(gca,'TickLabelInterpreter','latex','FontSize',10);

subplot(3,1,2);
hold on
pl = plot(x,curf);
pl.Color = [0 0 0];
yticks([0 2.5 5.0]);
axis([0 (length(x1)/sr) 0 5]);

ylabel('modulation current (mA)','Interpreter','latex','FontSize',10);
set(gca,'TickLabelInterpreter','latex','FontSize',10);
subplot(3,1,3);
pl = plot(x,norf);
pl.Color = [0 0 0];
axis([0 (length(x1)/sr) 0 .5]);
ylabel('normal load (N)','Interpreter','latex','FontSize',10);
xlabel('time (s)','Interpreter','latex','FontSize',10);
set(gca,'TickLabelInterpreter','latex','FontSize',10);
set(gcf,'color','w');