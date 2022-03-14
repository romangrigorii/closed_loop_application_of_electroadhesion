%% 9203 kistler piezo force sensor
sensor_param = 40; % units = pC/N
sensitivity = 100; % units = 100pC/mu
scale = .1;        % units = .1 mu/volt

conver = sensitivity*scale/sensor_param; % units = N/V