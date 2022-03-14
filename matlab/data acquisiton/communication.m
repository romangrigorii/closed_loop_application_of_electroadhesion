%% DISCONNECT AND CLEAR %%
pic_disconnect();
clear;
%% CONNECT %%
port = connect();
%% setting up DAQ
t = 6;
sr = 50000;
freq = 250;
car = sin(linspace(0,t,t*sr)*20000);
mod = 2.5*(sin(linspace(0,t,t*sr)*freq)+1)/2;
sig = (mod.*car)';
st = daq_connect(t,sr,sig);
%% zeroing microcontroller and matlab variables
vel = {};
pos = {};
port_engage(port,'z');
cfreq = 200;
res = 5/10000;
%% initialize motion paramters %% 
velocity = -50; %mm/s
transitionT = .2;
displacement = 75;
ttotal = motion_init(port, velocity, transitionT, displacement);
%pause(3);
%% starting motion %%
port_engage(port,'a');
port_engage(port,'z');
port_engage(port,'b');

port_engage(port,'r');
port_engage(port,'e');

%data = startForeground(st);
array_size = str2num(port_read(port,'s'));
t = (0:array_size-1)/cfreq;
vel_rec = zeros(array_size,1);
vel_des = zeros(array_size,1);
for i= 1:array_size
    out = str2num(fscanf(port));
    vel_rec(i) = out(1);
    vel_des(i) = out(2);
end

%% initiating automatic motion and data collection
vel_vec = [25,75,50,100];
freq_vec = [10,20,50,100];
transitionT = .2;
displacement = 75;
port_engage(port,'a');
port_engage(port,'z');
port_engage(port,'b');
send_and_collect(port, amp_vec, freq_vec, vel_vec,transitionT,displacement);


%% miscellenous

g = 6;
vel{g,1} = vel_rec;
vel{g,2} = vel_des;
pos{g,1} = integralr(vel_rec,0,6,0);
pos{g,2} = integralr(vel_des,0,6,0);





hold on
for g = 1:6
    a = plot((pos{g,1}-(max(pos{g,1}) + min(pos{g,1}))/2)*res ,vel{g,1}*res);
    a.Color = [1 0 0];
    a = plot((pos{g,2}-(max(pos{g,2}) + min(pos{g,2}))/2)*res,vel{g,2}*res);
    a.Color = [0 0 0];
end


hold on
a = plot(t,vel_rec*res);
a.Color = [1 0 0];
a = plot(t,vel_des*res);
a.Color = [0 0 0];
pos_des = integralr(vel_des,0,6,0);



fprintf(port,'k');
v = fscanf(port);
v = fscanf(port);
fprintf(port,'p');
fprintf(port, .15);
prop_constants = str2num(fscanf(port))


port_read(port,'d')

pic_disconnect();