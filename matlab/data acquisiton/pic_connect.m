function out = pic_connect(channel)
pic_disconnect();
port = serial(channel, 'BaudRate', 230400, 'FlowControl', 'hardware');
fopen(port);
out = port;
end