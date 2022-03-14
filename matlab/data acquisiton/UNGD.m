%% UNDG = universal negative group delay
%  this code tests out various paramters for a n order high pass filter
%  cascaded with a second order transfer function utilized to vary phase
%  and group delays to minimize ignal distortion
sr = 4000;
sn = logspace(-1,1,10);
sd = logspace(-3,1,10);
delfreq = 5:15;
cofreq = 5:15;
coord = 1:3;
trials = 50;
rr = zeros(length(sn),length(sd),length(delfreq),length(cofreq),length(coord),trials);
dbstopband = zeros(length(sn),length(sd),length(delfreq),length(cofreq),length(coord),trials);
dbpassband = zeros(length(sn),length(sd),length(delfreq),length(cofreq),length(coord),trials);
varofdel = zeros(length(sn),length(sd),length(delfreq),length(cofreq),length(coord),trials);
maxofdel = zeros(length(sn),length(sd),length(delfreq),length(cofreq),length(coord),trials);
maxpha = zeros(length(sn),length(sd),length(delfreq),length(cofreq),length(coord),trials);

[b,a] = butter(2,15*2/sr,'high');
for t = 1:trials
    t
    sigw = randn(sr,1);
    sigp = pink_filter(1,1000,sr,sigw);
    sig1 = filtfilt(b,a,sigp);
    for i = 1:length(sn)
        for j = 1:length(sd)
            for k = 1:length(delfreq)
                for l = 1:length(cofreq)
                    for p = 1:length(coord)                                           
                        [b1,a1] = butter(coord(p),2*cofreq(l)/sr,'high');
                        [h1,j1] = grpdelay(b1,a1,sr);
                        freq1 = j1/(2*pi)*sr;
                        del1 = h1/sr.*freq1*360;
                        
                        b2 = [1 2*(2*pi*delfreq(k)*sn(i)) (2*pi*delfreq(k))^2];
                        a2 = [1 2*(2*pi*delfreq(k)*sd(j)) (2*pi*delfreq(k))^2];
                        [b2,a2] = tfdata(c2d(tf(b2,a2),1/sr));
                        b2 = cell2mat(b2);
                        a2 = cell2mat(a2);    
                        
                        sig2 = filter(b1,a1,sigp);
                        sig2 = filter(b2,a2,sig2);
                        
                        rr(i,j,k,l,p,t) = corr(sig1',sig2');
                        
                        [h2,j2] = grpdelay(b2,a2,sr);
                        freq2 = j2/(2*pi)*sr;
                        del2 = h2/sr.*freq2*360;
                        
                        varofdel(i,j,k,l,p,t) = var(del1(41:2001)+del2(41:2001));
                        maxofdel(i,j,k,l,p,t) = max(abs(del1(41:2001)+del2(41:2001)));
                        
                        [j1,~,~] = freqz(b1,a1,sr);
                        [j2,~,~] = freqz(b2,a2,sr);
                        jj = j1.*j2;
                        dbstopband(i,j,k,l,p,t) = 20*log10(abs(jj(2)));
                        dbpassband(i,j,k,l,p,t) = 20*log10(abs(jj(41)));
                        
                        maxpha(i,j,k,l,p,t) = max(abs(unwrap(atan2(imag(jj),real(jj)))*180/pi));
                        
                    end
                end
            end
        end
    end
end
save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\misc\filter design\filter_res.mat','varofdel','maxpha','rr','dbstopband','dbpassband','sr','sn','sd','delfreq','cofreq','trials');