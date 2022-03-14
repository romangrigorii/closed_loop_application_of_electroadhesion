function out = pf(f1,f2,n,sr,t,o)
outpf = zeros(sr*(t+2),1);
outwf = zeros(sr*(t+2),1);
f = linspace(f1,f2,n);
tt = linspace(0,(t+2),(t+2)*sr);
for i = 1:n
    r = randn(1,1);
    outpf = outpf + sin(2*pi*tt*f(i)+ r)'/(f(i)^o);
    outwf = outwf + sin(2*pi*tt*f(i)+ r)';
end
outpf = outpf/std(outpf);
outpf = outpf(sr+1:end-sr);
outwf = outwf/std(outwf);
outwf = outwf(sr+1:end-sr);
out = [outwf,outpf];
end