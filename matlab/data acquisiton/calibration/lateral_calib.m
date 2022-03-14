% sr = 100000;
% st = daq.createSession('ni');
% st.Rate = sr;
% st.DurationInSeconds = 20;
% ch1 = addAnalogInputChannel(st,'Dev2','ai1','Voltage');
% ch1.InputType = 'SingleEnded'; %% lateral force
% x = input('');
% out = startForeground(st);
i = 1;
tr = .5;
out1 = detrend(out);
imps = {};
p = 0;
while i<= length(out)
    while abs(out1(i))<tr
        i = i + 1;
    end
    while (abs(out1(i))>abs(out1(i-1)))
        i = i + 1;
    end
    i = i-1;
    if abs(out1(i))>tr
        while (abs(out1(i))>abs(out1(i-1)))
            i = i -1;
        end
        p = p + 1;
        imps{p} = out1(i:i+.75*sr);
    end
    i = i + .25*sr;
end
% save('C:\Users\atrox\Desktop\Work\Research\projects\Friction control\publishable material\lat freq response\impulse_responses.mat','imps')
