% Numerical Integration of an age-structured model for testing during an
% epidemic, including a processing time delay. 

% The focus of the model is SARS-CoV-2, for which there is (1) a highly
% sensitive PCR test, for which test results are often delayed by days, and
% (2) the potential to manufacture a significantly less sensitive test that
% gives essentially instant results.

% The model allows us to examine the trade-off between sensitivity and
% frequency in a public health setting.

% In this code, the model is solved, as in FordeCiupeMain.m for several
% combinations of test sensitivity and testing capacity based on a fixed
% testing budget.


% Set global parameters for values that may change in loop

global loop_marker tau_1a tau_2a tau_3a tau_4a tau_1s tau_2s tau_3s tau_4s; 
global capacity testingDelay;

% Set global parameters used in the in-host model 
global beta1 c k1 beta2 k2 pi1 pi2 Gamma delta10 delta20;


%Initializing the program timer and loop counter
tic;
loop_marker=0;

% Calling parameter values
p=FordeCiupeParameters();

% Patient Data

% Generating a patient viral profile based on the within-host model of 
% Ke, et al. (BioRxiv)

% mean patient parameter values
c=10;
k1=4;
k2=4;
Gamma=0.01;
beta1=51.35e-8;
delta10=1.98;
pi1=49.8;
beta2=70.67e-8;
delta20=0.53;
pi2=0.34;

% patient initial values
T10 = 4.8e+6;
E10=0;
I10 = 10;
V10=0;
T20=4.8e+8;
E20=0;
I20=1;
V20=0;
y0 = [T10 E10 I10 V10 T20 E20 I20 V20];

% Final time
Tf = 20;

% Solving the ODE for a general patient 
options2=odeset('NonNegative',1);
[t y] = ode15s('covid_two_Ke', [0 Tf], y0);

% Determining the initial age of infection at which virus is detectable and
% the final age after which it becomes undetectable, dependent on the
% sensitivity level of the test. (Based on the viral profile generated
% above.)

% Defining the sensitivity levels to be tested
sensitivity_levels=linspace(2,5,7);

% Pre-allocation
tau1_values=NaN(length(sensitivity_levels),1);
tau4_values=NaN(length(sensitivity_levels),1);

for counter=1:length(sensitivity_levels)

    % Searching forward from 0 until the sensitivity level is reached
    for k=1:length(y(:,4))
        if log10(y(k,4))>=sensitivity_levels(counter)
            tau1_values(counter)=t(k);
            break
        end
    end
    % Searching backward from end until the sensitivity level is reached
    for k=1:length(y(:,4))
        ktwo=length(y(:,4))-k+1;
        if log10(y(ktwo,4))>=sensitivity_levels(counter)
            tau4_values(counter)=t(ktwo);
            break
        end
    end
end

% The parameters below may be changed in a loop of the solver code

% Viral Timing parameters
tau_1a = p.tau_1a;    % a.o.i at which virus becomes detectable
tau_2a = p.tau_2a;    % a.o.i at which patient starts being infectious
tau_3a = p.tau_3a;    % a.o.i on which patient stops being infectious
tau_4a = p.tau_4a;    % a.o.i at which virus stops being detectable

tau_1s = p.tau_1s;    % a.o.i at which virus becomes detectable
tau_2s = p.tau_2s;    % a.o.i at which patient starts being infectious
tau_3s = p.tau_3s;    % a.o.i on which patient stops being infectious
tau_4s = p.tau_4s;    % a.o.i at which virus stops being detectable


% Defining the lengths of test return delay to be simulated
testingDelayArray=0.5:0.5:5;

% Pre-allocation
CaseData=NaN(length(testingDelayArray),length(tau1_values));

% Calculated testing capacity associated with each sensitivity levels
% for a fixed testing budget
capacity_array=[0.1 0.1459 0.2336 0.3962 0.5858 0.8321 1];

% Nested loops to determine the final epidemic size at Day 150 for each
% combination of test sensitivity (and associated capacity) with test
% return delay
for testingCounter=1:length(testingDelayArray) 
    for tauCounter=1:length(tau1_values)        
    
    % Setting the delay, test sensitivity and capacity for the current run
    testingDelay=testingDelayArray(testingCounter);
    tau_1s=tau1_values(tauCounter);
    tau_1a=tau_1s;
    tau_4s=tau4_values(tauCounter);
    tau_4a=tau_4s;
    capacity=capacity_array(tauCounter);

    [S,I_a,I_s,Pos,Cases]=FordeCiupeSolver_loop(p);
    
    % Recording the outcome of the current run
    CaseData(testingCounter,tauCounter)=Cases(end);
    
    % Simulation progress counter. Outputs every 20 runs.
    simulationNumber=(testingCounter-1)*length(tau1_values)+tauCounter;
    if mod(simulationNumber,20)==0
        simulationNumber
    end

%hold on
    end
end

figure(1)
HeatmapName=heatmap(sensitivity_levels,testingDelayArray,CaseData,  'XLabel', 'Test sensitivity (log10 virus titer)', 'YLabel', 'Test return delay (days)', 'FontSize', 12);
colormap(flipud(hot))
caxis(HeatmapName,[0,0.7]);
set(gca,'FontSize', 12)
grid off

% End of timer
toc


