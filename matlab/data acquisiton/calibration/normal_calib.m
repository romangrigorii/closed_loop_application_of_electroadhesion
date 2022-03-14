%% collects and gives calibration paramters from the normal force readings
homedir = 'C:\Users\atrox\Desktop\Work\Research\projects\Friction control\';
sr = 100000;
weights = [0,10,20,30,40,50,70,100];
st = daq.createSession('ni');
st.Rate = sr;
st.DurationInSeconds = 5;
ch1 = addAnalogInputChannel(st,'Dev2','ai9','Voltage');
ch1.InputType = 'SingleEnded'; %% lateral force
m = weights;
[b,a] = butter(2,2/sr*50,'low');
for w = 1:length(weights)
     x = input('');
     out = startForeground(st);
     out = filtfilt(b,a,out);
     m(w) = mean(out);    
end
%save(strcat(homedir,'code\matlab\data acquisiton\calibration\calib_info\NF.mat'));

%% finding stop bands to filter out the force signal 

out1 = out;
stop_band_c = [12,21.8,38.3];
for i = 1:length(stop_band_c)
   [b,a] = butter(1,2/sr*[stop_band_c(i),stop_band_c(i)].*[.98,1.02],'stop');
    out1 = filtfilt(b,a,out1);
end
[b,a] = butter(6,2/sr*50,'low');
out1 = filtfilt(b,a,out1);

%% plotting norforce calib info
tsize = 10;
a = figure;
a.Position(3) = 200;
a.Position(4) = 150;
a = plot(weights/1000*9.8,-m'+5.449,'.-');
a.MarkerSize = 10;
a.Color = [0 0 0];
axis([0 1 0 7]);
xlabel('normal load (N)','Interpreter','latex','FontSize',tsize);
ylabel('output voltage (V)','Interpreter','latex','FontSize',tsize);
set(gca,'TickLabelInterpreter','latex','FontSize',tsize);
set(gcf,'color','w');
text(.1,5,'$R^{2} >.99$','Interpreter','latex','FontSize',tsize);
box off