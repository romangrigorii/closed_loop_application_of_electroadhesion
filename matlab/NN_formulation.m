function out = NN_formulation(input,labels,levels,plott)
[Li,~] = size(input);
in = abs(max(labels)-min(labels)); % calculating the range of the labels
levels = [levels,1]; %% force a single output topology
a = {}; a{1} = input(1,:)'; b = {}; W = {}; s = {}; z = {};
L = length(levels);
%Forward propogating initial values
for i = 1:L
    b{i} = in*(rand(length(levels(i)),1));
    W{i} = in*(rand(levels(i),length(a{i})));
    z{i} = W{i}*a{i} + b{i};
    a{i+1} = sigmoidd(z{i},0);
end
% Back propagation
gradCw = {}; gradCb = {}; C = []; nu = Li/L; lim = 100; % number of epochs

for p = 1:lim
    Ccur = []; % current cost
    for j = randperm(Li)
        a{1} = input(j,:)';
        
        % backpropagation
        s{L} = (a{L+1}-labels(j)).*sigmoidd(z{L},1);
        for i = flip(1:(L-1))
            s{i} = (W{i+1}.'*s{i+1}).*sigmoidd(z{i},1);
        end
        for i = 1:L
            gradCb{i} = s{i};
            gradCw{i} = s{i}*a{i}';
        end
        for i = 1:L
            W{i} = W{i} - .0005*nu*gradCw{i};
            b{i} = b{i} - .0005*nu*gradCb{i};
        end
        for i = 1:L % forward propagation the newly applied weights
            z{i} = W{i}*a{i} + b{i};
            a{i+1} = sigmoidd(z{i},0);
        end
        err = [];
        for j = 1:Li
            err(j) = classifyNN(W,b,input(j,:)') - labels(j);
        end
        Ccur = [Ccur,mean(abs(err))];
    end
    fprintf('percent: %d error: %f\n',round(100*p/lim),Ccur(end));
    %fprintf('percent: %d \n',round(100*p/lim));
    for i = 1:L
        W{i} = .0001*nu*(rand(size(W{i}))) + W{i}; % applying stochasticity
        b{i} = .0001*nu*(rand(size(b{i}))) + b{i};
    end
    C = [C,mean(Ccur)];
    close
    hold on
    plot(labels);
    plot(err + labels);
    pause(.5);
end
if plott == 1
    figure
    plot(C);
    xlabel('epoch');
    ylabel('error');
    title('Classification error as function of epochs on data set 2');
end
out{1} = W;
out{2} = b;
end