%% studying resonance behavior of tribometer

k1 = 2000000; % stiffness value of normal force sensors (both)
k2 = 3200; % stiffness value of leaf springs (both)
m1 = .7; % mass of bottom block
m2 = .027; % mass of mount
B = 100;
m0 = .3;
k0 = 10000;
hold on

[b,a] = ss2tf([[0 1/m1 0 0 ];[-k1 0 -k2 0];[0 1/m1 0 -1/m2];[0 0 k2 0]],[[0];[1];[0];[0]],[0 0 0 1],[0]);
[j,v,l] = bode(b,a);
plot(l/(2*pi),log2(j))
set(gca, 'XScale', 'log')

[b,a] = ss2tf([[-k1/B -1/m1 0 0 ];[k1 0 -k2 0];[0 1/m1 0 -1/m2];[0 0 k2 0]],[[0];[0];[0];[1]],[0 0 0 1],[0]);
[j,v,l] = bode(b,a);
plot(l/(2*pi),log2(j))
set(gca, 'XScale', 'log')

[b,a] = ss2tf([[-B/m0 -k1 0 0 0 -k0];[1/m0 0 -1/m1 -0 0 0];[0 k1 0 -k2 0 0];[0 0 1/m1 0 -1/m2 0];[ 0 0 0 k2 0 0];[1/m0 0 0 0 0 0]],[[0];[0];[0];[0];[1];[0]],[0 0 0 0 1 0],[0]);
[j,v,l] = bode(b,a);
plot(l/(2*pi),log2(j))
set(gca, 'XScale', 'log')