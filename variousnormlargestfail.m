% This part of code gives us the average number of failures per 100
% iterations, when we set the company with the largest asset value to fail
% initially, with respect to different values of sigma. Specifically, it is
% used to compare the functions w.r.t. fixed c, theta, and various
% sigma

clearvars

n=100; %number of organizations
iterations=100; %number of iterations for each set of parameter values


% We will consider various parameter values.
%c -- the proportion of an organization that is cross held. 
%theta -- the fraction of an organization's initial value it needs to retain to stay solvent
%d -- the expected degree of an organization (the expected number of other organizations it has cross holdings of and expected number of other organizations that hold it).
%h -- the expected value of each organization's property asset

c_range = [.7:.2:.8];
theta_range = [.9:.03:.92];
d_range = [1:1/2:20];
std_range = [0:0.6:2.4];
AvFailuresStd = zeros(length(d_range), length(std_range));

% We set up arrays to be populated with our simulation output.

AvFailures=zeros(length(d_range),length(theta_range)*length(c_range)); % keeping track of d's, theta's and c's
AvSomeFailures=zeros(length(d_range),length(theta_range)*length(c_range)); %Counter on fraction of iterations with some failures

% These are 3d versions of the above arrays.

AvFailures3d=zeros(length(d_range),length(theta_range),length(c_range));
AvSomeFailures3d=zeros(length(d_range),length(theta_range),length(c_range));

% We will loop over all parameters. First we loop over values of integration c.

std_counter = 0;
for std = std_range
    std_counter = std_counter + 1;
    AvFailures=zeros(length(d_range),length(theta_range)*length(c_range)); % keeping track of d's, theta's and c's
    AvSomeFailures=zeros(length(d_range),length(theta_range)*length(c_range)); %Counter on fraction of iterations with some failu
    AvFailures3d=zeros(length(d_range),length(theta_range),length(c_range));
    AvSomeFailures3d=zeros(length(d_range),length(theta_range),length(c_range));
    c_counter = 0;
    for c=c_range
	    c_counter=c_counter+1;
    
	    % Next we loop over values of theta.
        
	    theta_counter=0; % counter on theta increments
	    
        for theta=theta_range
		    theta_counter=theta_counter+1;
		    
		    % Then we loop over average degrees
		    
            d_counter=0; % counter on average degree increments
		    
            for d=d_range
			    d_counter=d_counter+1;
			    
                % We (re)set this varaibles to zero. For the current parameters, they will record the number of iterations in which at least one organization failed and the total number of failures over all iterations. 
                
                SomeFailures=0; %set the initial counter of cases where some organization fails to 0
			    TotFailures=0; %set the initial counter of total failures to 0
    
			    %The following for loop generates a different random graph for each iteration.
			    
			    for k=[1:iterations]
    
				    %We generate a random graph with n vertices and expected degree d, self-links forbidden; 
				    %First we generate a n by n matrix of iid uniform(0,1) random variables. We then define a threshold t and generate a link from i to j if the associated random number is below the threshold. Finally, we remove self links by setting the diagonal entries to zero.
				    %double(x) converts x from a logical type (which is spit out when we do comparisons, etc.) to numeric (double precision).
				    
				    U=rand(n,n);
				    t=d./(n-1);
				    G=double(U<t).*(ones(n,n)-eye(n));
                    
                    % We randomly generate for each organization a value of
                    % property asset, based on poisson distribution
                    p = normrnd(1,std,n,1);
                    p = max(p,0.01);
    
                    % We generate the asset matrix accordingly                
                    % p = poissrnd(5,n,1);
                    L = zeros(n);
                    for i = 1:n
                        L(:,i) = G(:,i).*p;
                    end
    
                    normL = zeros(n);
                    for i = 1:n
                        normL(:,i) = L(:,i) ./ (sum(L(:,i))+double(sum(L(:,i)) == 0));
                    end
    
				    %Now we can create the crossholdings matrix. To do this we give outside shareholders a share 1-c of each organization and distribute the remaining share of each organizations evenly among is crossholders.
				    %Note that this cross holdings matrix differs from the one defined in the paper as it includes the shares held by outside investors -- i.e. it corresponds to C + hat C from the paper.
                    
                    C = c.*normL;
    
				    %We now find the associated A matrix using the equation derrived in the paper. The difference in the equation is due to the different definition of C. 
				    
                    hatC=(1 - c).*eye(n);
				    A=hatC*(inv(eye(n)-C));
    
				    %Now we find the initial values and set the bankruptcy thresholds relative to these. Initial values are calculated assuming that no organization has failed and are given by A*ones(n,1) because each organization holds a single proprietary asset with value 1. 
                    %Note if an absolute threshold has instead been defined, many banks could be failing even before a shock.
				    
                    underlinev=(A*p).*theta;
    
				    %We now select an underlying asset and drop its value to 0. 
                    %Then we calculate the new organizational values and see if anyone else fails.
                    %We repeat this process until no new organizations fail. 
                    %Failure is modeled as an organization's underlying asset value going to zero. In other words bankrutpcy costs consume all the underlying asset value. Recall that bankrutpcy costs, like the underlying assets, are distributed among the organizations according to the A matrix. Without loss of generality, we always set asset 1's value to 0.
				    
                    Dold=p; %This is the vector of initial underlying asset values.
				    D=Dold;
                    [temp_val, temp_ind] = max(D);
				    D(temp_ind)=0; %We being the cascade by setting organization 1's underlying asset value to zero.
				    
				    %The following loop is our algorithm for calculating the minimum failure set.
				    
                    while ~isequal(D,Dold)%We keep looping until the vectors D and Dold are equal -- such that no new organization has just failed.
					    v=A*D; %These are the current values of organizations after the latest reductions to underlying asset values.
					    Dold=D;
					    D = D.*double(v>underlinev); %We set the underlying asset value to zero for any organizations who's value has fallen below their failure threshold.
                    end
                    
                    %For the current interation we record the number of failures.
                    
                    numberOfFailures = sum(v <= underlinev);
				    
				    %We add these to the total number of failures so far across all iterations at the current parameter values.
				    
                    TotFailures=TotFailures + numberOfFailures;
				    
				    %We also keep track of the number of iterations in which at least one organization has failed.
				    
                    SomeFailures = SomeFailures + double(numberOfFailures > 0); 
    
                end
                
                %To see how the simulations are progressing we print the current set of parameter values once we have done all the iterations for it.
                
                c
			    theta
			    d
                
			    %We now record our results in the arrays we set up earlier.
                
			    AvFailures3d(d_counter,theta_counter,c_counter)=TotFailures./iterations;
			    AvSomeFailures3d(d_counter,theta_counter,c_counter)=SomeFailures./iterations;
    
		    end
    
    
	    end
	    
    end
    AvFailuresStd(:,std_counter) = AvFailures3d(:,:,1);
end
%%
plot(AvFailuresStd,'-o')
legend('sigma = 0', 'sigma = 0.6', 'sigma = 1.2', 'sigma = 1.8', 'sigma = 2.4')
title('Average Failures over 100 Iterations with theta = 0.7, c = 0.7')
