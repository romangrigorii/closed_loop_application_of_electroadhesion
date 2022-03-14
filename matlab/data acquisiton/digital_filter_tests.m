%% digital filter tests
addpath(genpath('C:\Users\atrox\Desktop\Work\Research\Code\General code\MATLAB code'));
% init 

cofr1 = .5; % low frequency of the low-pass filter
cofr2 = 5;  % high frequency of the low-pass filter
dcofr1 = 7700; % remove
dcofr2 = 3200; % boost

dur = 10;
sr = 20000;
cors = zeros(2,1000);

% filt param compute

ac1 = [1 2*(2*pi*cofr2*sqrt(2)) (2*pi*cofr2)^2];
bc1 = [1 0 0];
[bd1,ad1] = stoz(bc1,ac1,sr);

ac2 = [1 2*(2*pi*cofr2*.5) (2*pi*cofr2)^2];
bc2 = [1 2*(2*pi*cofr1*.5) (2*pi*cofr1)^2];
[bd2,ad2] = stoz(bc2,ac2,sr); 

ac3 = [1 2*(2*pi*dcofr1*.01) (2*pi*dcofr1)^2];
bc3 = [0 0 1];
[bd3,ad3] = stoz(bc3,ac3,sr);

ac4 = [1 2*(2*pi*dcofr2*.2) (2*pi*dcofr2)^2];
bc4 = [1 2*(2*pi*dcofr2*1) (2*pi*dcofr2)^2];
[bd4,ad4] = stoz(bc4,ac4,sr);

[bcf,acf] = tfdata(tf(bc2,ac2)*tf(bc4,ac4));
bcf = cell2mat(bcf);
acf = cell2mat(acf);
[bdf,adf] = stoz(bcf,acf,sr);

% running correlation tests

for n = 1:1000
    n
    sig = randn(1,dur*sr);
    sig = pink_filter(.1,10000,sr,sig);
    [bb,aa] = butter(3,2/25000*[20,1000],'bandpass');
    sigg = filter(bb,aa,sig);
    sig_recorded = ifft(fft(sigg).*interp1((0:(length(ff1)-1))/(length(ff1))*100000,mean(ff1).',linspace(0,sr,dur*sr)));
    sig_recorded = 2*real(sig_recorded(1:dur*sr));
    
    out =  filteriir(bd1,ad1,sig_recorded);    
    out1 = filteriir(bd2,ad2,out);      
    out2 = filteriir(bd3,ad3,out1);  
    out3 = filteriir(bdf,adf,sig_recorded);
    % out4 = filteriir2([bd1;bd2;bd3],[ad1;ad2;ad3],sig_recorded);
    
    cors(1,n) = corr(out',sigg');
    cors(2,n) = corr(out3',sigg');
    cors(3,n) = corr(sig_recorded',sigg');  
    % cors(4,n) = corr(out4',sigg');  
end

% plotting

tsize = 12;
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
a = figure;
a.Position(3) = 600;
a.Position(4) = 300;
hold on

a = histogram(cors(3,:),.90:.0001:1);
a.EdgeColor = [.6 .6 .6];
a.FaceColor = [.6 .6 .6];

a = histogram(cors(1,:),.90:.0001:1);
a.EdgeColor = [0 0 0];
a.FaceColor = [0 0 0];

a = histogram(cors(2,:),.90:.0001:1);
a.EdgeColor = [1 0 0];
a.FaceColor = [1 0 0];

a = plot(median(cors(3,:)),-20,'.');
a.Color = [.6 .6 .6];
a.MarkerSize = 20;

a = plot(median(cors(1,:)),-20,'.');
a.Color = [0 0 0];
a.MarkerSize = 20;

a = plot(median(cors(2,:)),-20,'.');
a.Color = [1 0 0];
a.MarkerSize = 20;

a = plot([.951 1],[0 0]);
a.Color = [1 1 1];
a.LineWidth = .5;

a = legend('original vs recorded signal','original vs recorded signal high passed','original vs recorded signal special filtered');
a.Box = 'off';
a.Location = 'northwest';
title('white noise','Interpreter','latex','FontSize',tsize);

xlabel('Pearson''s r','Interpreter','latex','FontSize',tsize);
ylabel('count','Interpreter','latex','FontSize',tsize);
set(gcf,'color','w');
axis([.95,1,-20,475])
xticks([.95 .96 .97 .98 .99 1])