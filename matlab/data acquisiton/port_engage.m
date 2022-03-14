function port_engage(port,option)
fprintf(port,option);
out = fscanf(port);
end
