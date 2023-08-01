%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binomial Method for Calculating Employee Stock Option%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

tic
%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=0.05; % interest rate
sigma=0.3; % Volatility
T=10; % expiration time
K=50; % Strike Price
s0=50; % Stock Price
N=500; % Number of Time Step
v=3; % Vesting Period
M=2; % The parameter M is the ratio of the stock price to the strike price necessary to trigger voluntary early exercise
exr=0.06; % Exit Rate
D=0.025; % Dividen yield
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate Other Parameters
delta=T/N;
u=exp(sigma*sqrt(delta));
d=exp(-sigma*sqrt(delta));
p=(exp((r-D)*delta)-d)/(u-d);
fd = exp(-r*delta);
pe = 1-exp(-exr*delta);
S(1,1)=s0;
id = N;

% Build Binomial Tree
for j=2:N+1
    S(1,j)=u*S(1,j-1);
end
for i=2:N+1
    for j=i:N+1
        S(i,j)=d*S(i-1,j-1);
    end
end

% ESO Calculation
O = zeros(N+1,1);
OB = zeros(N+1,1);

for i=1:N+1
    O(i,N+1)=max(S(i,N+1)-K,0);
end

for j=1:N
    for i=1:N+1-j
        if (N-j)*delta>v
            if S(i,N-j+1)< M*K
                O(i,N-j+1)=(exp(-exr*delta)*exp(-r*delta)*(p*O(i,N-j+2)+(1-p)*O(i+1,N-j+2)))+((1-exp(-exr*delta))*max(0,S(i,N-j+1)-K));
            else
                O(i,N-j+1)= S(i,N-j+1)-K;
            end
        else
            O(i,N-j+1)=exp(-exr*delta)*exp(-r*delta)*(p*O(i,N-j+2)+(1-p)*O(i+1,N-j+2));    
        end
    end
end
disp(['ESO Price = ',num2str(O(1))])

for j=1:N+1
    OB(j,N+1)=max(S(j,N+1)-K,0);
end

for j=1:N
    for i=1:N+1-j
        OB(i,N-j+1)=exp(-r*delta)*(p*OB(i,N-j+2)+(1-p)*OB(i+1,N-j+2));    
    end
end

disp(['European Call Option Price = ',num2str(OB(1))])

% Plot Binomial Tree
k=1;
for i=1:N+1
for j=1:i
t(k) = i-1 ;
S1(k)=S(j,i);
k=k+1;
end
end
figure (1)
scatter(t,S1,9,'bd','filled')
hold on

xlabel('t');
ylabel('S','rotation',0);
title ({
    ['Binomial CRR Tree']
    ['r = ' num2str(r) ', \sigma = ' num2str(sigma) ', T = ' num2str(T) ', S_0 = ' num2str(s0) ', K = ' num2str(K) ', D = ' num2str(D) ', M = ' num2str(M) ', N = ' num2str(N)]
    ['']})
grid on
toc