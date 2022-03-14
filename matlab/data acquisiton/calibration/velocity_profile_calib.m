addpath(genpath('C:\Users\atrox\Desktop\Work\Research\projects\Texture player\code\matlab'));
L = 800;
vel_rec = zeros(L,1);
vel_des = zeros(L,1);
srfreq = 62500;
all_data = {};
%% finding baseline zero
%input = ('press enter if ready');

direction = 'r';
port = connect();
vels = [1,2,3];
port_write(port,'d',vels);
port_engage(port,'i');
port_write(port,'v',direction);

data = startForeground(st);
release(st);
for i = 1:length(vels)
    st = daq_connect(5,srfreq);
    port_write(port,'v','r');
    port_engage(port,'r');
    port_engage(port,'e');
    out1 = startForeground(st);
    array_size = str2num(port_read(port,'s'));
    for ii= 1:array_size
        out = str2num(fscanf(port));
        vel_rec(ii) = out(1);
    end
    all_data{i,1} = vel_rec;
end
end