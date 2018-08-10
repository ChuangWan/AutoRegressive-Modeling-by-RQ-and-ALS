function BetaHat = CAViaR_estim(y, THETA)

% *****************************************************************************************
% Set parameters for optimisation.
% *****************************************************************************************
REP			  = 30;          % Number of times the optimization algorithm is repeated.
nInitialVectors = [1000, 3]; % Number of initial vector fed in the uniform random number generator.
nInitialCond = 1;            % Select the number of initial conditions for the optimisation.
MaxFunEvals = 500;           % Parameters for the optimisation algorithm. 
MaxIter     = 500;           % Increase them in case the algorithm does not converge.

%options = optimset('LargeScale', 'off', 'HessUpdate', 'dfp', 'LineSearchType', 'quadcubic','MaxFunEvals', ...
                    %MaxFunEvals, 'display', 'off', 'MaxIter', MaxIter, 'TolFun', 1e-6, 'TolX', 1e-6);
options = optimset('MaxFunEvals',MaxFunEvals, 'display', 'off', 'MaxIter', MaxIter, 'TolFun', 1e-6, 'TolX', 1e-6);
warning('off', 'verbose')
%   
%**************************** Optimization Routine ******************************************  
%
% Compute the empirical THETA-quantile for y.
ysort = sortrows(y(1:100), 1); 
empiricalQuantile = ysort(round(100*THETA));

initialTargetVectors = unifrnd(0, 1, nInitialVectors);
RQfval = zeros(nInitialVectors(1), 1);
for i = 1:nInitialVectors(1)
    RQfval(i) = RQobjectiveFunction(initialTargetVectors(i,:), 1, y, THETA, empiricalQuantile);
end
Results          = [RQfval, initialTargetVectors];
SortedResults    = sortrows(Results,1);

BestInitialCond  = SortedResults(1:nInitialCond,2:4);    
Beta = zeros(size(BestInitialCond)); fval = Beta(:,1); exitflag = Beta(:,1);

for i = 1:size(BestInitialCond,1)
    [Beta(i,:), fval(i,1), exitflag(i,1)] = fminsearch('RQobjectiveFunction', BestInitialCond(i,:), ...
        options, 1, y, THETA, empiricalQuantile);
    for it = 1:REP
        if exitflag(i,1) == 1, break, end
        [Beta(i,:), fval(i,1), exitflag(i,1)] = fminsearch('RQobjectiveFunction', Beta(i,:), ...
            options, 1, y, THETA, empiricalQuantile);
        if exitflag(i,1) == 1, break, end
    end
end
SortedFval  = sortrows([fval, Beta, exitflag, BestInitialCond], 1);
    
BetaHat   = SortedFval(1, 2:4)';
if SortedFval(1,5)~=1, disp('Warning: CAViaR convergence not achieved.'), end
%**************************** End of Optimization Routine ******************************************
   