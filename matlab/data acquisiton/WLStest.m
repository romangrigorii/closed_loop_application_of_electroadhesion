sr = 10000;
sig = randn(1,11*sr);
[b1,a1] = butter(1,2*1000/10000,'low');
[b2,a2] = butter(1,2*20/10000,'high');
sig = filtfilt(b1,a1,filtfilt(b2,a2,sig));
sig = sig(5001:end-5000);
[bl3,al3] = butter(2,10*2/sr,'low');
sig1 = filtfilt(bl3,al3,sige);
t = linspace(0,10,10*sr);

ff2 = zeros(25,100);
for j = 1:25
    [b,a] = fir1(10*j*2,2*1/sr,'low');
    for ii = 1:100
        sige = sin(2*pi*.5*t+.3) + sig*ii/100;        
        sig2 = filter(b,a,sige);
        sig3 = sig2;
        for i = (j*10+1):length(sig2)-(j*10)
            for k = 0:(j*10)-1
                P1 = sig3(i-(j*10)+k);
                P2 = sig3(i+k);
                P3 = P2 + (P2-P1)/j/10;
                sig3(i+1) = P3;
            end
        end
        ff2(j,ii) = corr(sig',(sige-sig3)').^2;
    end
    j
end

sr = 10000;
[bl1,al1] = butter(2,5*2/10000,'low');
[bl2,al2] = butter(3,5*2/10000,'low');
for q = 1:20
    lat = DATA{q}(:,1)*.25;
    ii = 1;
    while abs(lat(ii))<.1
        ii = ii+10;
    end
    lat = lat(ii:end);    
    lat = resample(lat,sr,60000);    
    sig1 = filtfilt(bl1,al1,lat);
    sig2 = filter(bl2,al2,lat);
    sig1r = resample(sig1,100,sr);
    sig2r = resample(sig2,100,sr);
    
    lat1 = sig1>0;
    lat1 = abs(diff(lat1)).*(1:length(lat1)-1)';
    lat1 = lat1(lat1~=0);
    
    for j = 1:50
        [b,a] = fir1(5*j*2,2*.1/sr,'low');
        sig2 = filter(b,a,lat);
        sig2r = resample(sig2,100,5000);
        sig3 = sig2;
        for i = (j*5+1):length(sig2)-(j*5)
            for k = 1:(j*5)-1
                P1 = sig2(i-(j*5)+k);
                P2 = sig3(i);
                P3 = P2 + (P2-P1)/j/5;
                sig3(i) = P3;
            end
        end
        for l = 1:length(lat1)-1
            ff2(j,q) = ff2(j,q)+corr(lat(lat1(l):lat1(l+1))-sig1(lat1(l):lat1(l+1)),lat(lat1(l):lat1(l+1))-sig3(lat1(l):lat1(l+1)));
        end
        ff2(j,q) = ff2(j,q)/(length(lat1)-1);
    end
    q
end

[k,l,j] = freqz(b,a,sr);

sigf = filter(b,a,sig);
dsigf = diff(sigf);
dsigff = dsigf;
for i = 3:length(dsigf)
    if (abs(dsigf(i-1)-dsigf(i-2))*5)<abs(dsigf(i)-dsigf(i-1))
        dsigff(i) = dsigf(i-1);
    end
end


h = 2;
X = [[1 -h h^2 -(h^3)];[1 -(h-1) (h-1).^2 (-(h-1)^3)];[1 0 0 0];[1 h-1 (h-1)^2 (h-1)^3];[1 h h^2 h^3]];    
K = ((h+1).^2-(j.^2)).*((h+2)^2 - (j.^2)).*((h+3)^2-(j.^2));
K = K/sum(K);

Bs = inv((X.'*(K)*X))*X.'*K;


h = 20;
j = -h:h;
K = ((h+1).^2-(j.^2)).*((h+2)^2 - (j.^2)).*((h+3)^2-(j.^2));
K = K/sum(K);
X = [[1 -h h^2 -(h^3)];[1 -(h-1) (h-1).^2 (-(h-1)^3)];[1 0 0 0];[1 h-1 (h-1)^2 (h-1)^3];[1 h h^2 h^3]];    

S = zeros(size(j));
for s = 1:length(j)
    S(s) = sum(K.*(j.^(s-1)));
end
w = K.*((S(5)-(S(3).*(j.^2)))/(S(1)*S(5) - (S(3).^2)));

