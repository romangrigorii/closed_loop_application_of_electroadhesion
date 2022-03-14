function pic_disconnect()
if ~isempty(instrfind)
    fclose(instrfind);
end
end