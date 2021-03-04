function [S,I_a,I_s,Pos,Cases] = FordeCiupeSolver_loop(p)

global tau_1a tau_2a tau_3a tau_4a tau_1s tau_2s tau_3s tau_4s capacity testingDelay;

% Pre-allocation of solution arrays
S =   NaN(1,p.Q+1);
Pos= NaN(1,p.Q+1);
Cases=NaN(1,p.Q+1);
                                
I_a   = NaN(p.K+1,p.Q+1);
I_s   = NaN(p.K+1,p.Q+1);

% Pre-allocation of age and time dependent function arrays
rho_a   = NaN(p.K+1,p.Q);
rho_s   = NaN(p.K+1,p.Q);

lambda_a = NaN(p.K+1);
lambda_s = NaN(p.K+1);

% Initial conditions
S(1) = p.S0;
Pos(1)=p.Pos0;
Cases(1)=p.I0;

% Initial conditions for Infection class, age-dependent
I_a(1,1) = (1-p.f)*p.I0/p.dtau;
I_s(1,1) = p.f*p.I0/p.dtau;
for j = 2:p.K+1         
    I_a(j,1)=0;
    I_s(j,1)=0;
end

% Define force of infection multiplier functions (lambda_a and lambda_s)

for k=1:(p.K+1)
    if ((k >= floor(tau_2a/p.dtau)) && (k<= floor(tau_3a/p.dtau)))
        lambda_a(k)=p.gamma;
    else lambda_a(k)=0;
    end 
    if ((k >= floor(tau_2s/p.dtau)) && (k<= floor(tau_3s/p.dtau)))
        lambda_s(k)=1;
    else lambda_s(k)=0;
    end 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Numerical integration of the iterated system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Multiplicative factors for loss to death during either a time step,
%    or during the testing delay

testloss_a = exp(-(p.mu+p.m_a)*testingDelay);  % death of infected patients
                                               % during testing delay
testloss_s = exp(-(p.mu+p.m_s)*testingDelay);  % death of infected patients
                                               % during testing delay
steploss_a = exp(-(p.mu+p.m_a)*p.dtau);        % death of infected patients
                                               % during a time step
steploss_s = exp(-(p.mu+p.m_s)*p.dtau);        % death of infected patients
                                               % during a time step
                                               
for q=1:p.Q         % time steps
   
    InfInt=0;         % integrates new infections over age classes
    DetInt=0;         % integrates new quarantines over age classes
    
    % If the all infection classes are empty, no further infections can
    % occur, and the remainder of the simulation can be skipped
    
    dont_skip_flag=1;
    if (sum(I_a(:,q))==0) && (sum(I_s(:,q))==0)
        dont_skip_flag=0;
    end
    
    
    % Calculate the testing rates for this time point, all ages
      
    % Assuming a random testing rate over all populations
    current_pop=S(q)+p.dtau*(sum(I_a(:,q))+sum(I_s(:,q)));
    rho_current = -log(1-capacity/current_pop);
        if (~isreal(rho_current))|(rho_current==-Inf)
            rho_current=NaN;
        end
        
    % restricting the testing to be successful only for ages of infection
    % within the range determined by the test sensitivity
    for k=1:(p.K+1)
        if (k >= floor(tau_1a/p.dtau)) && (k<= floor(tau_4a/p.dtau))
            rho_a(k,q)=rho_current;
        else rho_a(k,q)=0;
        end    
        if (k >= floor(tau_1s/p.dtau)) && (k<= floor(tau_4s/p.dtau))
                rho_s(k,q)=rho_current;
        else rho_s(k,q)=0;
        end
        
    end
    
    if (dont_skip_flag==1)
    
    % Updating infect ages classes other than the first,
    % both symptomatic and asymptomatic patients
    for k=1:p.K
        
        % carrying forward empty age classes     
        if (I_a(k,q)==0) && (I_s(k,q)==0)
            I_a(k+1,q+1)=0;
            I_s(k+1,q+1)=0;
            continue;
        end
        
        % Define the number of time steps for testing delay
        L=floor(testingDelay/p.dtau);
        
        % Calculate the number removed from the current age class by tests
        % conducted L times steps prior

        if (k <= L) | (q <= L)
            omega_a=0;
            omega_s=0;
        else
            omega_a=rho_a(k-L,q-L)*I_a(k-L,q+1-L)*testloss_a;
            omega_s=rho_s(k-L,q-L)*I_s(k-L,q+1-L)*testloss_s;
        end
        
        if (I_a(k,q)==0)
            I_a(k+1,q+1)=0;
            omega_a=0;
        else
            % If the detection rate had been NaN, the entire population
            % segment gets a positive test result, leaving the class           
            if isnan(omega_a)
                I_a(k+1,q+1)=0;
                omega_a=I_a(k,q)*steploss_a;
            
            % Otherwise, moving a population to the next age, next time
            % step, using an implicit scheme 
            else
                I_a(k+1,q+1) = I_a(k,q)*steploss_a-(omega_a/(p.mu+p.m_a))*(1-steploss_a);
                if (I_a(k+1,q+1)<0) | (isnan(I_a(k+1,q+1)))
                    I_a(k+1,q+1)=0;
                    omega_a=I_a(k,q)*steploss_a;
                end
            end
        end
        
        if (I_s(k,q)==0)
            I_s(k+1,q+1)=0;
            omega_s=0;
        else     
            % If the detection rate had been NaN, the entire population
            % segment gets a positive test result, leaving the class           
            if isnan(omega_s)
                I_s(k+1,q+1)=0;
                omega_s=I_s(k,q)*steploss_s;
            
            % Otherwise, moving a population to the next age, next time
            % step, using an implicit scheme 
            else
               I_s(k+1,q+1) = I_s(k,q)*steploss_s-(omega_s/(p.mu+p.m_s))*(1-steploss_s);
                if I_s(k+1,q+1)<0 | (isnan(I_s(k+1,q+1)))
                   I_s(k+1,q+1)=0;
                   omega_s=I_s(k,q)*steploss_s;
                end
            end
        end
        
        % Updating Detection numbers of capacity exceeds population, or not
        DetInt=DetInt+p.dtau*(omega_a+omega_s);
        
        % New infections caused by this age class added to the overall infection integral 
        InfInt=InfInt+p.dtau*(lambda_a(k+1)*I_a(k+1,q+1)+lambda_s(k+1)*I_s(k+1,q+1));
        
    end
    
    else
        % Simpler update if infection class is empty
        for k=1:p.k
            I_a(k+1)=0;
            I_s(k+1)=0;
        end
        InfInt=0;
    end
    
    % Update the Susceptible class (implicit scheme)
  
    S(q+1)=(S(q)+p.b*p.dt)/(1+p.dt*(p.mu+p.beta*InfInt));
        
    % Increment the cumulative positive test cases
  
    Pos(q+1)= Pos(q)+DetInt*p.dt;
        
    % Update first infection age-class
    
    I_a(1,q+1)=(1-p.f)*p.beta*S(q+1)*InfInt;
    I_s(1,q+1)=p.f*p.beta*S(q+1)*InfInt;
    
    % Update the total number of infections
    
    Cases(q+1)= Cases(q)+(I_a(1,q+1)+I_s(1,q+1))*p.dt;
    
end


end
