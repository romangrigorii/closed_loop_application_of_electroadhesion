function ttotal = motion_init(port,velocity,transitionT,displacement)
ttotal = (transitionT*2) + abs(displacement/velocity);
res = 5/10000;
cfreq = 200;
velocity = velocity/res;
displacement = displacement/res;

port_write(port,'v',velocity + .000000001);
port_write(port,'t',transitionT + .000000001);
port_write(port,'l',displacement + .000000001);
fprintf(strcat('motion params set to\n',num2str(str2num(port_read(port,'x')).*[res,1,res]),'\n'));
end