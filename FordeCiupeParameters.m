function [p] = FordeCiupeParameters()

% Model rate parameters
p.mu = 4e-5;      % background per capita death rate
p.b = p.mu;       % constant birth rate of susceptibles.
p.beta = 0.25;    % infection rate (constant across age of infection)
p.m_a = 0.0001;   % excess per capita death due to infection (asymptomatic)
p.m_s = 0.0001;   % excess per capita death due to infection (symptomatic)
p.gamma = 0.7;    % relative infectiousness of asymptomatic patient 
                  % (compared to symptomatic)
p.f = 0.7;        % portion of infections resulting in symptomatic cases

% Viral profile and testing parameters
%PCR   LOD=10^2                      % a.o.i = age of infection

% %Asymptomatic patients
% p.tau_1a = 0.544;    % a.o.i at which virus becomes detectable
% p.tau_2a = 2.5;    % a.o.i at which patient starts being infectious
% p.tau_3a = 10.5;    % a.o.i on which patient stops being infectious
% p.tau_4a = 10.95;   % a.o.i at which virus stops being detectable
 
% % Symptomatic patients
% p.tau_1s = 0.544;    % a.o.i at which virus becomes detectable
% p.tau_2s = 2.5;    % a.o.i at which patient starts being infectious
% p.tau_3s = 10.5;    % a.o.i on which patient stops being infectious
% p.tau_4s = 10.95;   % a.o.i at which virus stops being detectable

 
% Test parameters

% PCR
% p.testingDelay = 5;      % time required to get test results back
% p.capacity = 0.1;        % number of tests available per day

% RDT, Antigen
p.testingDelay = 0.5;      % time required to get test results back
p.capacity = 0.1429;       % number of tests available per day

% Paper-strip
% p.testingDelay = 0.1;    % time required to get test results back
% p.capacity = 0.1;        % number of tests available per day

%Antigen test LOD=10^5
%Asymptomatic patients
p.tau_1a = 2.77;    % a.o.i at which virus becomes detectable
p.tau_2a = 2.5;    % a.o.i at which patient starts being infectious
p.tau_3a = 10.5;    % a.o.i on which patient stops being infectious
p.tau_4a = 7.37;   % a.o.i at which virus stops being detectable

% Symptomatic patients
p.tau_1s = 2.77;    % a.o.i at which virus becomes detectable
p.tau_2s = 2.5;    % a.o.i at which patient starts being infectious
p.tau_3s = 10.5;    % a.o.i on which patient stops being infectious
p.tau_4s = 7.37;   % a.o.i at which virus stops being detectable
 
%Paper test LOD=10^6
%Asymptomatic patients
% p.tau_1a = 3.48;    % a.o.i at which virus becomes detectable
% p.tau_2a = 2.5;    % a.o.i at which patient starts being infectious
% p.tau_3a = 10.5;    % a.o.i on which patient stops being infectious
% p.tau_4a = 6.14;   % a.o.i at which virus stops being detectable
 
% % Symptomatic patients
% p.tau_1s = 3.48;    % a.o.i at which virus becomes detectable
% p.tau_2s = 2.5;    % a.o.i at which patient starts being infectious
% p.tau_3s = 10.5;    % a.o.i on which patient stops being infectious
% p.tau_4s = 6.14;   % a.o.i at which virus stops being detectable
 

% Initial conditions

p.I0 = 0.01;            % intiial infected population (a.o.i. = 0)
p.S0 = p.b/p.mu-p.I0;   % initial susceptible population
p.Pos0 = 0;             % initial number of positive tests

% Maximum time
p.finalTime = 150;% 150; %200;%5*365;

% Maximum age of infection
p.finalTau = 150;% 150; %200;

% Size of grid step
p.dt = .1; % 
p.dtau = p.dt;

% Number of time steps
p.Q = floor(p.finalTime / p.dt);

% Number of tau steps
p.K = floor(p.finalTau / p.dtau);

% Age and Age Step for the end of infected class
p.time_in_class=14;  % Length of time patients are considered to be infected 
p.ClassEnd=floor(p.time_in_class/p.dtau); % age step for end of infected time
