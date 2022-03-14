%%dataloc = 'C:\Users\Roman\Desktop\Work\Research\Publications\Munich Conference paper\';
%%run('C:\Users\Roman\Desktop\Work\Research\General code\MATLAB code\WHC2017\params.m')
%%run('C:\Users\Roman\Desktop\Work\Research\General code\MATLAB code\HEADER.m');
addpath('C:\Users\atrox\Desktop\Work\Research\Code\General code\MATLAB code');
speds = [80 120 180 270];
dataloc = 'C:\Users\atrox\Desktop\Work\Research\projects\z Finished\High-bandwidth tribometry as means of capturing natural texture\NU tribometry\updated_data\';
materials = dir(dataloc);
materials = materials(3:end-1);
srNU = 125000;
[bn,an] = butter(3,20*2/srNU,'low');
[bp,ap] = butter(3,20*2/srNU,'low');

dL = round(srNU/10);
fun = @(x,xdata)x(1)*xdata + x(2);

params = {};
FF = {};

latfarr = {}; 
k = 1;
for m = 1:length(materials)
    material = materials(m).name;
    fprintf('Working on: %s\n', material(1:end-2));
    for s = 1:4
        speed = speds(s);
        fprintf('     speed: %d\n', speed);
        trialz = dir(strcat(dataloc,material,'\Finger\', num2str(speed)));
        trialz = trialz(3:end);
        q = 1; qq = 1;
        arr = []; arr2 = [];
        for t = 1:length(trialz)-1
            triall = trialz(t).name;
            % dispstat(strcat('Percent done:  ',num2str(round(100*t/(length(trialz)-1)))));
            content = dir(strcat(dataloc,material,'\Finger\', num2str(speed),'\',triall));
            content = content(3:end);
            dirr = strcat(strcat(dataloc,material,'\Finger\', num2str(speed),'\',triall,'\',content(length(content)).name));
            load(dirr);
            counterNBits = 32;
            signedThreshold = 2^(counterNBits-1); %Correcting Position data to have negative values
            signedData = data(:,2);
            signedData(signedData > signedThreshold) = signedData(signedData > signedThreshold) - 2^counterNBits;
            position = signedData*0.003837; %Converting Position data to mm
            position = filtfilt(bp,ap,position);
            
            vel = derivR(position,1,srNU);
            normalForce = data(:,1);
            normalForce = filtfilt(bn,an,normalForce);
            
            nFcpeace = detrend(normalForce(1:srNU*.4));
            nFc = mean(normalForce(1:srNU*.4) - nFcpeace); %Zeroing normal force at beginning of data collection for 1s
            normalForce = (normalForce - nFc)*3.58207; %Converting normal force to N
            lateralForce = data(:,3);
            lFc = mean(lateralForce(1:1000)); %Zeroing lateral force at beginning of data collection for 1s
            lateralForce = (lateralForce - lFc)*(-0.943169*2.33084); %Converting lateral force to N
            vellim = abs(vel)>50;
            vellim2 = derivR(vellim,1,srNU)>0;
            vellim2 = vellim2.*(1:length(vellim))';
            vellim2 = vellim2(vellim2~=0);
            vellim3 = derivR(vellim2,1,2);
            vellim2 = vellim2(vellim3>mean(vellim3));
            vellim3 = [diff(vellim2)>mean(diff(vellim2));0==1];
            k = 1;
            for i = 1:length(vellim2)-1
                if vellim3(i)==1
                    velc = vel(vellim2(i):vellim2(i+1));
                    norc = normalForce(vellim2(i):vellim2(i+1));
                    latc = lateralForce(vellim2(i):vellim2(i+1));
                    if (t==2)
                        if (k==1)
                            latfa{m,s} = latc;
                        end
                    end
                    k = k + 1;
                    L = length(latc);
                    N = floor(L/dL);
                    ff = abs(fft(latc)/L*2);
                    f = (0:(length(ff)-1))/(length(ff))*srNU;
                    ff = interp1(f,ff,1:2000);
                    FF{m,s,2}(qq,1:2) = [mean(norc),mean(velc)];
                    FF{m,s,1}(qq,1:2000) = ff;
                    qq = qq + 1;
                    for n = 1:N
                        latcp = latc((n-1)*dL+1:n*dL);
                        norcp = norc((n-1)*dL+1:n*dL);
                        velcp = velc((n-1)*dL+1:n*dL);
                        arr(1:3,q) = [mean(latcp),mean(norcp),mean(velcp)];
                        q = q + 1;
                    end
                end
            end
            params{m,s} = arr;
        end
    end
end

save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\nat texture\tex_params.mat','params','materials','FF');

%% extracting correlation paramters
cor = [];
P = [];
tsize = 11;
siz = [];
for m = 1:length(materials)
    for s = 1:4
        for c = 1:2000
            [cor(m,s,c,1),P(m,s,c,1)] = corr(FF{m,s,1}(FF{m,s,2}(:,1)>.05,c),FF{m,s,2}(FF{m,s,2}(:,1)>.05,1));
            [cor(m,s,c,2),P(m,s,c,2)] = corr(FF{m,s,1}(FF{m,s,2}(:,1)>.05,c),FF{m,s,2}(FF{m,s,2}(:,1)>.05,2));
        end
        siz(m,s) = length(FF{m,s,2});
    end
end

v1 = 1:2000;
v2 = logspace(0,log10(2000),2000);
cc = [];
hold on
for m = 1:length(materials)
    for s = 1:4
        cc(m,s,1:2000) = smooth(interp1(v1',squeeze(cor(m,s,:,1)),v2'),20);
    end
end

%% plotting correlation stuff
cmap = colormap('hot');
tsize = 11;
a = figure;
materialz = {'Empire Velveteen','Faux Leather','Hucktowel parallel','Hucktowel perpendicular','Denim','Microsuede','Swimwear'};
a.Position(3) = 1500;
a.Position(4) = 200;
for m = 1:length(materials)
    subplot(1,7,m);
    if (m == 1)
        ylabel('Pearson''s   r','Interpreter','latex','FontSize',tsize);
    end
    if (m == 4)
        xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);
    end
    title(materialz{m},'Interpreter','latex','FontSize',tsize-1);
    hold on
    a = plot([.01 1000],[.5 .5],'--');
    a.Color = [0 0 0];
    a.LineWidth = .5;
     a = plot([.01 1000],[0 0],'--');
    a.Color = [0 0 0];
    a.LineWidth = .5;   
    for s = 1:4
        a = plot(v2((P(m,s,:,1)<=.05)),squeeze(cc(m,s,(P(m,s,:,1)<=.05),1)),'.');
        % a.Color = 1 - [1 1 1]/5*s;
        a.Color = cmap(round((6-s)/7*63),:);
        %a.LineWidth = 1;
        a.MarkerSize = 5;
        if ~isempty(squeeze(cor((P(m,s,:,1)>.05))))
            f = 1:2000;
            a = plot(v2((P(m,s,:,1)>.05)),squeeze(cc(m,s,(P(m,s,:,1)>.05),1)),'.');
            a.Color = [.9 .9 .9];
            a.MarkerSize = 5;
        end
    end
    set(gca, 'XScale', 'log')
    axis([1 1000 -.25 1])
    xticks([1 10 100 1000]);
    set(gca,'TickLabelInterpreter','latex','FontSize',tsize-1);
end
set(gcf,'color','w');
%%

hold on
for m = 1:length(materials)-1
    for s = 1:4
        if ~isempty(squeeze(cor((P(m,s,:,1)>.01))))
            f = 1:2000;
            a = plot(f((P(m,s,:,1)>.01)),squeeze(cc(m,s,(P(m,s,:,1)>.01),1)),'.');
            a.Color = [1 0 0];
            a.MarkerSize = 3;
        end
    end
end
set(gca, 'XScale', 'log')
axis([1 1000 -.1 1])
axis([1 1000 -.25 1])
set(gcf,'color','w');
xlabel('$log_{10}$(P value)','Interpreter','latex','FontSize',tsize);
ylabel('count','Interpreter','latex','FontSize',tsize);
xlabel('frequency (Hz)','Interpreter','latex','FontSize',tsize);

%%
a = figure;
a.Position(1) = a.Position(1) - 500;
a.Position(2) = a.Position(2) - 400;
a.Position(4) = 800;
a.Position(3) = 1600;
tsize = 11;
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
materialz = struct2cell(materials);
for i = 1:7
    subplot(2,7,i);
    hold on
    for s = 1:4
        a = plot(params{i,s}(2,:),params{i,s}(1,:),'.');
        a.Color = [.6 .6 .6];
        pause(1);
    end
    xlabel('normal force (N)','Interpreter','latex','FontSize',tsize);
    ylabel('lateral force (N)','Interpreter','latex','FontSize',tsize);
    title(materialz{1,i}(1:end-2),'Interpreter','latex','FontSize',tsize);    
    axis([0 1.5 -1 1])
    subplot(2,7,i+7);
    hold on
    for s = 1:4
        a = plot(params{i,s}(1,:)./params{i,s}(2,:),params{i,s}(3,:),'.');
        a.Color = [.6 .6 .6];
        pause(1);
    end   
    xlabel('lateral force / normal force','Interpreter','latex','FontSize',tsize);
    ylabel('swipe velocity (mm/s)','Interpreter','latex','FontSize',tsize);
    axis([-2 2 -375 375])    
end
