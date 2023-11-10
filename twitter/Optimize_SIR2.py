import numpy as np
import pandas as pd
import numba 
from scipy.optimize import minimize
from Datenum import days_from_zero
import datetime

def SIR_rhs(del_t, SIR, pars, capacity):
    dSIR = np.zeros(3)
    alpha = pars[0]
    beta = pars[1]
    mu = pars[2]
    dSIR[0] = (-(alpha * SIR[0] * SIR[1]) + (mu * ((capacity - SIR[3] / capacity) * SIR[3]))) * del_t
    dSIR[1] = (alpha * SIR[0] * SIR[1] - beta * SIR[1]) * del_t
    dSIR[2] = (beta * SIR[1]) * del_t
    return dSIR

def SIR_Model(del_t, SIR0, params, grid_size, time, capacity):
    SIR = np.zeros((grid_size + 1, 5))
    SIR[0, :] = SIR0
    for i in range(1, grid_size + 1):
        dSIR = SIR_rhs(del_t, SIR[i - 1, :], params, capacity)
        SIR[i, 0] = SIR[i - 1, 0] + dSIR[0] #1st column
        SIR[i, 1] = SIR[i - 1, 1] + dSIR[1] #2nd column
        SIR[i, 2] = SIR[i - 1, 2] + dSIR[2] #3rd column
        SIR[i, 3] = SIR[i, 0] + SIR[i, 1] + SIR[i, 2] #4th column
        SIR[i, 4] = SIR[i - 1, 4] + del_t #5th column
        SIR[i, :] = np.where(SIR[i, :] < 0, 0, SIR[i, :])

    #time_temp = pd.to_datetime(time) - pd.to_datetime(time[0])
    #time_temp = time_temp / np.timedelta64(1, 'D')  # convert time difference to days
    time = pd.Series(time)
    time_temp = days_from_zero(time)
    time_temp = [item - time_temp[0] for item in time_temp]

   # I = np.zeros((len(time_temp), 5))
    I = np.zeros((len(time_temp), SIR.shape[1])) #Initialize I with zeroes
    for j in range(len(time_temp)):
        k = time_temp[j]
        row = np.where(SIR[:, 4] <= k)[0]
        if row.size >0: # Only proceed if 'row' is not empty
            I[j, :] = SIR[row[-1], :]
    return SIR, I

def error_sum_of_squares(dataI, I):
    ESS = np.sum((dataI - I[:, 1])**2)
    Rsquared = 1 - ESS / np.sum((dataI - np.mean(dataI))**2)
    return ESS, Rsquared
#time needs to go back to datetime 
def Optimize_SIR(time, dataI, dataS, tspan):
    A = 1E-10
    B = 1E4
    n = 500
    LA = np.log10(A); LB = np.log10(B)
    alpha = 10**(LA + (LB-LA) * np.random.rand(n, 1))               # Infectious rate
    beta = 10**(LA + (LB-LA) * np.random.rand(n, 1))               # recovery rate
    mu = -10**(LA + (LB-LA) * np.random.rand(n, 1)) # growth rate
    params = np.column_stack((alpha, beta, mu))
    dataN = np.zeros(len(time))  # Initialize dataN as a zero array
    dataN[0] = np.array(dataS[0] + dataI[0])  # Assign the first element
    capacity = max(dataS)
    del_t = .001
    grid_size = round((1/del_t)*tspan)
    time_model = np.zeros(grid_size + 1)
    for j in range(1, grid_size + 1):
        time_model[j] = time_model[j-1] + del_t
    SIR0 = [dataS[0], dataI[0], 0, dataN[0], 0]  # initial conditions
    ESS = np.zeros(n)
    Rsquared = np.zeros(n)
    for i in range(n):
        print(f"Current iteration: {i}")  # Display current iteration
        SIR, I = SIR_Model(del_t, SIR0, params[i, :], grid_size, time, capacity)
        ESS[i], Rsquared[i] = error_sum_of_squares(dataI, I)
    idx = np.argmin(ESS)
    SIR_values = [alpha[idx], beta[idx], mu[idx], Rsquared[idx]]
    SIR, I = SIR_Model(del_t, SIR0, params[idx, :], grid_size, time, capacity)
    SIR_time = np.zeros(grid_size + 1)
    for j in range(1, grid_size + 1):
        SIR_time[j] = SIR_time[j-1] + del_t
    return SIR, ESS, SIR_values, SIR_time, Rsquared



'''
Optimize_SIR MATLAB
%% This script runs the SIR Model and matches it to the Data
function [SIR,ESS,SIR_values,SIR_time,Rsquared] = Optimize_SIR(time,dataI,dataS,tspan)
A = 1E-10;
B = 1E4;
n = 200;
LA = log10(A); LB = log10(B);
alpha           = 10.^(LA + (LB-LA) * rand(n,1));               % Infectious rate
beta            = 10.^(LA + (LB-LA) * rand(n,1));               % recovery rate 
% mu              = rand(n,1)*0;                  % growth rate
mu            = -10.^(LA + (LB-LA) * rand(n,1));
params          = [alpha,beta,mu];
dataN(1:length(time)) = zeros();
dataN(1)        = dataS(1) + dataI(1);
capacity        = max(dataS);
del_t           = .001;
grid_size       = round((1/del_t)*tspan,0);
time_model(grid_size+1) = zeros();
for j = 2:grid_size+1;
    time_model(j) = time_model(j-1)+del_t;
end
SIR0            = [dataS(1),dataI(1),0,dataN(1),0];  % initial conditions
for i = 1:n
    disp(['Current iteration: ' num2str(i)]);
    [SIR I]         = SIR_Model(del_t,SIR0,params(i,:),grid_size,time,capacity);
    [ESS(i),Rsquared(i)]  = error_sum_of_squares(dataI,I);
end
[M idx] = min(ESS);
% [M idx] = max(Rsquared);
SIR_values = [alpha(idx), beta(idx),mu(idx),Rsquared(idx)];
[SIR I]         = SIR_Model(del_t,SIR0,params(idx,:),grid_size,time,capacity);
SIR_time = time(1);
for j = 2:grid_size+1;
    SIR_time(j) = SIR_time(j-1)+del_t;
end
end


SIR_Model MATLAB
function [SIR I] = SIR_Model(del_t,SIR0,params,grid_size,time,capacity)
    SIR(grid_size+1,5)= zeros();
    SIR(1,:) = SIR0;
    for i = 2:grid_size+1
        dSIR=SIR_rhs(del_t,SIR(i-1,:),params,capacity);
        SIR(i,1)=SIR(i-1,1)+dSIR(1);
        SIR(i,2)=SIR(i-1,2)+dSIR(2);
        SIR(i,3)=SIR(i-1,3)+dSIR(3);
        SIR(i,4)=SIR(i,1)+SIR(i,2)+SIR(i,3);
        SIR(i,5)=SIR(i-1,5)+del_t;
        if SIR(i,:) < 0
            SIR(i,:) = 0;
        end
    end

    %% Find the values of the SIR model that corresponds to the time of the data
   time_temp = datenum(time);
   time_temp = time_temp-time_temp(1);
    for j = 1:length(time_temp)
        k = time_temp(j);
        [row column] = find(SIR(:,5)<=k);
        I(j,:) = SIR(row(end),:);
    end
% I = 1;
end


SIR_rhs MATLAB
% This function will solve the right hand side of the differential equations
% and provide the next timestep for the number of susc, inf, rec
% Nick Duncan
% 
% January 2018

function dSIR = SIR_rhs(del_t,SIR,pars,capacity)
dSIR = zeros(3,1);
alpha=pars(1);
beta=pars(2);
mu=pars(3);
% mu = 0;
%  dSIR(1) = (-(alpha*SIR(1)*SIR(2))+mu*SIR(1))*del_t; 
dSIR(1) = (-(alpha*SIR(1)*SIR(2))+(mu*((capacity-SIR(4)/capacity))*SIR(4)))*del_t; 
dSIR(2) = (alpha*SIR(1)*SIR(2) - beta*SIR(2))*del_t;
dSIR(3) = (beta*SIR(2))*del_t;
end






'''