function send_and_collect_2(port,amp_vec,freq_vec,vel_vec,transitionT,displacement,direction,num_brst,dirr)
L = length(vel_vec);
vel_rec = zeros(800,1);
vel_des = zeros(800,1);
all_data = {};
carfreq = 20000;
srfreq = 62500;
%% finding baseline zero
% st = daq_connect(2,srfreq,zeros(2*srfreq,1));
% data = startForeground(st);
% mn = mean(data(:,1));
% ml = mean(data(:,2));
% release(st);
ml = 0;
mn = 0;
x = input('ready to start\n');
L2 = L/num_brst;
for d = 1:2
    for i = 1:L2
        %% initialize
        velocity = direction*vel_vec(i); %mm/s
        ttotal = motion_init(port, velocity, transitionT, displacement)+.25;
        sig = create_sig(amp_vec((i-1)*num_brst+1:i*num_brst),freq_vec((i-1)*num_brst+1:i*num_brst),carfreq,srfreq,ttotal,2,15625);
        st = daq_connect(ttotal,srfreq,sig);
        
        %% start actuation and data aquisition
        
        port_engage(port,'r');
        port_engage(port,'e');
        data = startForeground(st);
        array_size = str2num(port_read(port,'s'));
        for ii= 1:array_size
            out = str2num(fscanf(port));
            vel_rec(ii) = out(1);
            vel_des(ii) = out(2);
        end
        all_data{i,1,d} = data(:,1)-mn; % normal force
        all_data{i,2,d} = data(:,2)-ml; % lateral force
        all_data{i,3,d} = data(:,3); % vibrometry
        all_data{i,4,d} = vel_des;
        all_data{i,5,d} = vel_rec;
        all_data{i,6,d} = data(:,4); % trigger
        direction = -1*direction;
        release(st);
        fprintf(strcat(num2str(i/L2*100),'\n'));
    end
end
save(strcat(dirr,'\output_data.mat'),'all_data');
end