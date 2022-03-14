function send_and_collect_1(port,amp_vec,freq_vec,vel_vec,transitionT,displacement,direction)
L = length(vel_vec);
vel_rec = zeros(800,1);
vel_des = zeros(800,1);
velIerr = zeros(800,1);
all_data = {};
carfreq = 20000;
srfreq = 62500;
%% finding baseline zero
st = daq_connect(2,srfreq,zeros(2*srfreq,1));
data = startForeground(st);
mn = mean(data(:,1));
ml = mean(data(:,2));
release(st);
x = input('ready to start\n');
for i = 1:L    
    %% initialize
    velocity = direction*vel_vec(i); %mm/s
    ttotal = motion_init(port, velocity, transitionT, displacement)+.25;
    sig = create_sig(amp_vec(i),freq_vec(i),carfreq,srfreq,ttotal,1);
    st = daq_connect(ttotal,srfreq,sig);
    
    %% start actuation and data aquisition  
    
    port_engage(port,'r');
    pause(.01);
    port_engage(port,'e');
    data = startForeground(st);
    pause(.01);
    array_size = str2num(port_read(port,'s'));
    pause(.01);
    for ii= 1:array_size
        out = str2num(fscanf(port));
        vel_rec(ii) = out(1);
        vel_des(ii) = out(2);
        velIerr(ii) = out(3);
    end
    all_data{i,1} = data(:,1)-mn; % normal force
    all_data{i,2} = data(:,2)-ml; % lateral force
    all_data{i,3} = data(:,3); % vibrometry
    all_data{i,4} = vel_des;
    all_data{i,5} = vel_rec;
    all_data{i,6} = data(:,4); % trigger
    all_data{i,7} = velIerr; % trigger
    direction = -1*direction;
    release(st);
    fprintf(strcat(num2str(i/L*100),'\n'));
end
save('C:\Users\atrox\Desktop\Work\Research\projects\Texture player\data\unprocessed data\datax.mat','all_data')
end