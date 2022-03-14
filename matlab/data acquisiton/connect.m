%% connects to the PIC
function port = connect()
port = pic_connect('COM6');
fprintf(port,'');
fprintf(fscanf(port));
end