% Monte Carlo Simulation Method for European Call and Put Option

clc
tic
%format short
format long
rng('default')

%%%%%%%%%% Parameters Value %%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$
S0 = 295; % Stock Price
K =300; % Strike Price
r = 0.019; %Interest Rate
sigma = 0.1361; % Volatility
T =24/360; % Expiration Time
M = 1e6; % number of simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$

% Monte Carlo Simulations
ST = S0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*randn(M,1));

Call = exp(-r*T)*max(ST-K,0);
Put = exp(-r*T)*max(K-ST,0);

aM_Call = mean(Call);
bM_Call = std(Call);

aM_Put = mean(Put);
bM_Put = std(Put);

CI_Call = [aM_Call - 1.96*bM_Call/sqrt(M) , aM_Call + 1.96*bM_Call/sqrt(M)];
CI_Put = [aM_Put - 1.96*bM_Put/sqrt(M), aM_Put + 1.96*bM_Put/sqrt(M)];

disp(['Monte Carlo Simulation Results with ',num2str(M), ' samples'])

disp([char(10) 'European Call Option Price = ',num2str(aM_Call)])
disp(['95% Confidence Interval for European Call Options = ',num2str(CI_Call)])

disp([char(10) 'European Put Option Price = ',num2str(aM_Put)])
disp(['95% Confidence Interval for European Put Options = ',num2str(CI_Put), char(10)])

% Black-Scholes
d1 = (log(S0/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
C_BS = S0*normcdf(d1) - K*exp(-r*(T))*normcdf(d2) ;
P_BS = C_BS + K*exp(-r*T) - S0 ;


disp([char(10) 'Black-Scholes European Call Option Price = ',num2str(C_BS)])

disp([char(10) 'Black-Scholes European Put Option Price = ',num2str(P_BS)])

% Plot Simulation 
histfit(ST,[],'normal'); % distribution of Stock Price
figure
plot(Call); % history of simulated call option prices
figure
plot(Put); % history of simulated put option prices
figure
plot(ST); % history of simulated stock prices

toc