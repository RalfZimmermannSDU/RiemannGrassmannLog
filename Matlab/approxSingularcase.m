% Approximation of the singular case
clear all;
close all;
disp("The code deals with singular cases.")
warningsYN = input("Do you want to see Matlab warnings? Enter 0/1:")

% save the original warning settings
orig_state = warning;

if warningsYN == 0
    warning('off');
end


n = 1000;
p = 200;
N = 10; % number of runs
r = 1; % number of singular angles
s = RandStream('mt19937ar','Seed',10); % Random stream for reproducability


T = logspace(-20,0,100);
SubspaceErrorLog1 = zeros(N,length(T));
SubspaceErrorLog2 = zeros(N,length(T));
SubspaceErrorLog3 = zeros(N,length(T));

for k = 1:N
    % Create random Stiefel representative U0 with orthogonal completion U0perp
        X = rand(s,n);
        [Q0,~] = qr(X);
        U0 = Q0(:,1:p);
        U0perp = Q0(:,p+1:n);
    % Create a random tangent vector with r largest singular values of pi/2
        B = rand(s,n-p,p);
        Delta = U0perp*B;
        [Q,~,V] = svd(Delta,0);
        S = @(t) diag(sort([pi/2*ones(1,r), pi/2*rand(s,1,p-r)],'descend')).*t;
        Delta = @(t) Q*S(t)*V';

    for i = 1:length(T)
        t = 1-T(i);
        
        % Calculate the subspace associated with the tangent vector
        Deltat = Delta(t);
        U1 = GrassmannExp(U0,Deltat);
        
        % Calculate the new and the standard log
        DeltaLog = GrassmannLogOneSVD(U0,U1);
        DeltaLog_standard = GrassmannLog_standard(U0,U1);
        
        % Project the result of the standard log algorithm onto the
        % horizontal space
        DeltaLog_standardproj = DeltaLog_standard - U0*(U0'*DeltaLog_standard);
        
        % Calculate the subspaces associated with the log results
        U1Log = GrassmannExp(U0,DeltaLog);
        U1Log_standardproj = GrassmannExp(U0,DeltaLog_standardproj);
        U1Log_standard = GrassmannExp(U0,DeltaLog_standard);
        
        % Calculate the subspace errors
        SubspaceErrorLog1(k,i) = subspaceDist(U1,U1Log);
        SubspaceErrorLog2(k,i) = subspaceDist(U1,U1Log_standardproj);
        SubspaceErrorLog3(k,i) = subspaceDist(U1,U1Log_standard);
    end
end


% Plot the results on a log-log plot
axes('XScale', 'log', 'YScale', 'log')
hold on

for k = 1:N
    plot(T,SubspaceErrorLog1(k,:),'*','color',[0, 0.4470, 0.7410]);
    plot(T,SubspaceErrorLog2(k,:),'x','color',[0.8500, 0.3250, 0.0980]);
    plot(T,SubspaceErrorLog3(k,:),'+','color',[0.9290, 0.6940, 0.1250]);
end
xlabel('{\tau}')
ylabel('Error by subspace distance')
legend('New log algorithm', 'Standard log alg. (with horiz. projection)', 'Standard log algorithm')

print -depsc approxSingularcase


% restore original warning state
warning(orig_state);