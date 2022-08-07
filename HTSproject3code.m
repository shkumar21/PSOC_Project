clear all;

%% Simulation Parameters

T = 5;          % Time in Hrs
N = 21;         % Samples
dt = T/(N-1);   % Length of samples
t = 0:dt:T;     % vector of time
M = 6;          % Optimization variables = number of constraints = pg1,pg2,pg3,pgh1,pgh2,pgh3 ==6
%cost functions for thermal power plants
c(1,:) = [ 0.005 11 200 ];
c(2,:) = [ 0.009 10 180 ];
c(3,:) = [ 0.007 10 230 ];

%minimum and maximum power constraints for Thermal and hydro power plants
Pmin = [50*ones(N,1), 37*ones(N,1), 45*ones(N,1), 43*ones(N,1), 43*ones(N,1), 43*ones(N,1)]; 
% ones(21,1)* 50 it means 50 repeated in 21 rows and one coloumn for all i
Pmax = [200*ones(N,1), 150*ones(N,1), 180*ones(N,1), 280*ones(N,1), 280*ones(N,1), 280*ones(N,1)]; 
										
Vmin = 7e6;% Min water volume of reservoir
Vmax = 8e6; % Max water volume of reservoir

Pload0 = 530;         % Nominal load in MW
Pload = Pload0*[0.50  0.53  0.55  0.53  0.5  0.54  0.70  0.90  0.95  1.10 1.20  1.40  1.70  1.65  1.50  1.30  1.00  0.90  0.80  0.50 0.54];

pmin = reshape(Pmin',M*N,1); % became (126 row,1 column) dimesnsion [ 50;37;45;43;43;43;.....;43]
pmax = reshape(Pmax',M*N,1); % became (126 row,1 column) dimesnsion [ 200;150;180;280;280;.....;280]

% example below q(opt) u can understand the concept. p(6,:) represent optimum values for Hydroelectric generation 
q = @(P) 320e3 + 6e3*P(6,:)+ exp(0.035*P(6,:));

% volume formula --- reshape(p,M,N) means that from 126*1 dimension make each row represent generation unit
V = @(p) sum(q(reshape (p,M,N)))*dt;

% nonlinear constraints water volume for each power 
noncon = @(p) deal([V(p)-Vmin; Vmax-V(p)], []);

Aload = kron(eye(N),ones(1,M));% linear equality constraints pg1+pg2+pg3+ph=pload ---- [1 1 1 1 1 1] dimension 21X126

%% Optimization

p0 = pmin;      %intial guess
% PERFORMING THE OPTIMIZATION OPERATION TO MINIMIZE COST
opt = optimset ('Algorithm', 'Interior-point', 'Display', 'iter', 'MaxIter', 1e5, 'MaxFunEvals', 1e5 );

[popt,fopt] = fmincon(@(p) costFun(p,c,M), p0, [], [], Aload, Pload', pmin, pmax,noncon, opt);


%% Plotting

% reshaping of vector in matrix [M,N]

popt1=popt; % lenght of popt1 126*1 
Popt = reshape(popt, M, N); % after reshape each row represent a generation unit and hydrounit 6*21 dimension
Pmint = reshape(pmin, M, N); % same as above
Pmaxt = reshape(pmax, M, N); % same as above
Pgen = Aload*popt; % Aload dimension = (N*126) X popt Dimension =(126*1) --> pgen dimension = (N*1)
qt = q(Popt); % popt= reshape(popt, M, N) this illustrate the above
Vt =cumsum(qt*dt);% multiply each value by dt 21 values

%% subplot

figure(1)

for i = 1:M
    subplot (3,2,i)
    plot(t,Pmint(i,:), 'r--',t,Pmaxt(i,:),'r--',t,Popt(i,:),'b'); % Popt after reshape
    title(['P',num2str(i)]);
    xlabel('dt')
    ylabel('Power P')

end
subplot (3, 2, 5)
plot(t,Vmin*ones(size(t)), 'r--',t,Vmax*ones(size(t)),'r--',t,Vt,'b'); % make vmin and Vmax in the length of T
title('V')
xlabel('dt')
ylabel('Volume V')
subplot(3,2,6)
plot(t,Pload,'rx',t,Pgen,'b')
legend('Pload', ' Pgen')
xlabel('dt')
ylabel('Power P')






