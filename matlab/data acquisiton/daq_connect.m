%% set-up DAQ for collecting data
function st = daq_connect(duration,sr,sig)
st = daq.createSession('ni');
st.Rate = sr;
st.DurationInSeconds = duration;
ch1 = addAnalogInputChannel(st,'Dev2','ai1','Voltage');
ch1.InputType = 'SingleEnded'; %% lateral force
ch2 = addAnalogInputChannel(st,'Dev2','ai9','Voltage');
ch2.InputType = 'SingleEnded'; %% normal force
ch3 = addAnalogInputChannel(st,'Dev2','ai2','Voltage');
ch3.InputType = 'SingleEnded'; %% DAC reading
ch4 = addAnalogInputChannel(st,'Dev2','ai10','Voltage');
ch4.InputType = 'SingleEnded'; %% SYN reading
ch5 = addAnalogOutputChannel(st,'Dev2','ao1','Voltage');
queueOutputData(st,sig);
end