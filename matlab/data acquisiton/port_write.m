function port_write(port,option,data)
fprintf(port,option);
v = fscanf(port);
switch option
    case 'd'
        fprintf(port,length(data)+.0000001);
        for i = 1:length(data)
            fprintf(port,data(i)+.0000001);
            %fscanf(port)
        end
    case 'v'
        fprintf(port,data);
end
end
