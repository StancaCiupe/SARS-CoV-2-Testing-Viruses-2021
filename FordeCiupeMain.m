% Numerical Integration of an age-structured model for testing during an
% epidemic, including a processing time delay. 

% The focus of the model is SARS-CoV-2, for which there is (1) a highly
% sensitive PCR test, for which test results are often delayed by days, and
% (2) the potential to manufacture a significantly less sensitive test that
% gives essentially instant results.

% The model allows us to examine the trade-off between sensitivity and
% frequency in a public health setting.

%Initializing the program timer
tic;

% Calling parameter values
p=FordeCiupeParameters();

% Solving the system and storing tracking statistics
% S = Susceptible population over time
% I_a, I_s = Infected asymptomatic and symptomatic populations, over time
%             and age of infection
% Pos = Cumulative number of positive test results returned to date
% Cases = Cumulative number of infections (detected or not) to date
% DailyCounts = Tracks the number of new infections and new positive test
%                results each day
% CurrentI_a, CurrentI_s = sum of I_a and I_s over age classes that are
%                           considered to be infected (day 0-14).

[S,I_a,I_s,Pos,Cases,DailyCounts,CurrentI_a,CurrentI_s]=FordeCiupeSolver(p);

% independent variable arrays for plotting

times=linspace(0,p.finalTime,p.Q+1);
days=linspace(1,p.finalTime,p.finalTime);
ages=linspace(0,p.finalTau,p.K+1);

% Generation of the figure

% Subplot 1: Currently infected populations over time
figure(2)
subplot(2,2,1);
plot(times, CurrentI_a,'b',times,CurrentI_s,'r', 'LineWidth',1);
legend('Asymptomatic','Symptomatic');
ylabel('Populations')
set(gca,'fontsize', 12)
axis([0 150 0 0.15])

% Subplot 2: Cumulative number of cases and of case detections over time
subplot(2,2,2);
plot(times, Cases,'m', times,Pos,'g', 'LineWidth',1);
legend('Total Cases', 'Detected Cases');
set(gca,'fontsize', 12);
axis([0 150 0 1])

% Subplot 3: New cases and new case detections reported each day
subplot(2,2,3:4);
bar(days,DailyCounts(2,:),'y','FaceAlpha',1);
hold on;
bar(days,DailyCounts(1,:),'b','FaceAlpha',0.5);
hold on
legend('Cases','Detections');
hold on
set(gca,'fontsize', 12);
xlabel('Days');
ylabel('Daily cases');
hold on;
axis([0 150 0 0.015])

% End of the program timer
toc


