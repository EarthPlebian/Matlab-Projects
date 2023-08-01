%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trinomial Kamrad-Ritchken Method for Calculating Employee Stock Option%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

tic

%%%%%%%%%%%%% Data-data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
B=M*K; % Barrier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate stretch parameter
dt=T/N;
NU=log(s0/B)/(sigma*sqrt(dt));
NO=floor(NU);

if NU == NO || NU < 1
    lambda = 1;
else
    lambda=NU/NO;
end

% Calculate Other Parameters
u=exp(lambda*sigma*sqrt(dt));
d=1/u;
miu = r - D - (sigma^2/2);
pu = (1/(2*lambda^2)) + ((miu*sqrt(dt))/(2*lambda*sigma));
pd = (1/(2*lambda^2)) - ((miu*sqrt(dt))/(2*lambda*sigma));
pm = 1 - (1/lambda^2);

% Build trinomial Kamrad-Ritchken Tree
S(1,1) = s0;
for j = 2:N+1
    S(1,j) = S(1,j-1) * u ;
    mm = 2*j-1 ;
    for i = 2:mm-1
        S(i,j) = S(i-1,j-1) ;
    end
    S(mm,j) = S(mm-2,j-1) * d ;
end

% ESO Calculation
for i=1:2*(N+1)-1
    O(i,N+1)=max(S(i,N+1)-K,0);
end

for j=N:-1:1
    for i=1:2*j-1
        if j*dt>v
            if S(i,j)< M*K
                O(i,j)=(exp(-exr*dt)*exp(-r*dt)*(pu*O(i,j+1)+pm*O(i+1,j+1)+pd*O(i+2,j+1)))+((1-exp(-exr*dt))*max(0,S(i,j)-K));
            else
                O(i,j)= S(i,j)-K;
            end
        else
            O(i,j)=exp(-exr*dt)*exp(-r*dt)*(pu*O(i,j+1)+pm*O(i+1,j+1)+pd*O(i+2,j+1));    
        end
    end
end
disp(['ESO Price = ',num2str(O(1))])

for j=1:2*(N+1)-1
    OB(j,N+1)=max(S(j,N+1)-K,0);
end

for j=N:-1:1
    for i=1:2*j-1
        OB(i,j)=exp(-r*dt)*(pu*OB(i,j+1)+pm*OB(i+1,j+1)+pd*OB(i+2,j+1));   
    end
end

disp(['European Call Option Price = ',num2str(OB(1))])

% Plot Trinomial K-R Tree
k = 1 ;
for j = 1:N+1
    mm = 2*j-1 ;
    for i = 1:mm
        t(k) = j-1 ;
        S1(k) = S(i,j) ;
        k = k + 1 ;
    end
end

plot(t,S1,'b.')
hold on
x = 0:0.005:N ;
y = B*ones(size(x)) ;
plot (x,y,'r-')
hold on
y = K *ones(size(x));
plot (x,y,'g-')
xlabel ('t') ;
ylabel ('S') ;
title ({
    ['Trinomial Kamrad-Ritchken Tree']
    ['r = ' num2str(r) ', \sigma = ' num2str(sigma) ', T = ' num2str(T) ', S_0 = ' num2str(s0) ', K = ' num2str(K) ', D = ' num2str(D) ', M = ' num2str(M) ', N = ' num2str(N)]
    ['']})
text (5,B+1,'Barrier')
grid on
toc