% ************************************************ %
% Experimental procedure is set-up in this program %
% ************************************************ %
addpath(genpath('C:\Users\atrox\Desktop\Work\Research\projects\Texture player\code\matlab'));

%% DISCONNECT AND CLEAR %%
delete(instrfind);
pic_disconnect();
clear;
%% CONNECT %%
port = connect();
%%
dur = 12;
sr = 62500;
%% %% finding No
st = daq_connect(1,sr,zeros(1,sr)');
out = startForeground(st);
mN = mean(out(:,2));
%% passing filtered noise data to the daq
VtoN = -.1815;
VtoL = .3;
Ptomm = 1/10;
out = randn(1,(dur+2)*sr);
[b,a] = butter(3,[20,1100]*2/sr);
out = filtfilt(b,a,out);
out = out(sr+1:end-sr);
sig = out/std(out);
st = daq_connect(dur+2,sr,sig');
x = input('ready?\n');
port_engage(port,'r');
out = startForeground(st);

lat = out(:,1);
nor = out(:,2);
[b,a] = butter(5,50*2/sr,'low');
nor = (filtfilt(b,a,nor-mN))*VtoN;
[b,a] = butter(3,2000*2/sr,'low');
lat = filtfilt(b,a,lat);
lat = lat*VtoL;

pos = zeros(1,2000);
vel = pos;
port_engage(port,'s');
for i = 1:2000
    temp = fscanf(port);
    temp = str2num(temp);
    pos(i) = temp(1);
    vel(i) = temp(2);
end
pos = pos*Ptomm;
vel = vel*Ptomm;

[b,a] = butter(3,15*2/200,'low');
velf = filtfilt(b,a,vel);

%% syncing data

veli = interp1(tim(1:end-1),vel,linspace(tim(1),tim(end),sr*10));
dmat = abs(diff(out(:,4)));
dmat = [dmat',0];
loc = (dmat>2.5).*(1:length(dmat));
loc = loc(loc~=0);
shiftt = 10*sr - loc(1);
veli = veli(shiftt:end);

%% processing
velw = [zeros(1,69000-1),veli(69000:end)];
vel_up = velw>50;
vel_down = velw<-50;
f_up = lat(vel_up)./nor(vel_up);
s_up = sig(vel_up);
f_down = lat(vel_down)./nor(vel_down);
s_down = sig(vel_down);

L = length(f_up);
f_up = interp1(linspace(0,L/sr,L),f_up,linspace(0,L/sr,L/sr*2500));
s_up = interp1(linspace(0,L/sr,L),s_up,linspace(0,L/sr,L/sr*2500));
L = length(f_up);

hold on
for i = 1:L
    a = plot((s_up(i)+5)/10*8,-f_up(i),'.');
    a.Color = (1-i/L)*[1,1,1];
end

L = length(f_down);
f_down = interp1(linspace(0,L/sr,L),f_down,linspace(0,L/sr,L/sr*2500));
s_down = interp1(linspace(0,L/sr,L),s_down,linspace(0,L/sr,L/sr*2500));
L = length(f_down);

hold on
for i = 1:L
    a = plot((s_down(i)+5)/10*8,f_down(i),'.');
    a.Color = (1-i/L)*[1,1,1];
end
