% comparison between Wirtinger Flow and Kaczmarz method for Gaussian 
clear all
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
% Problem size
n = 256;
m = round(4*n);  

maxiter = 100;
num_tests = 50;
relerrKac_array = zeros(maxiter+1,num_tests);
relerrRKac_array = zeros(maxiter+1,num_tests);

for ll = 1 : num_tests

% make signal and data: Complex
x = randn(n,1) + 1i*randn(n,1);  
A = 1/sqrt(2)*randn(m,n) + 1i/sqrt(2)*randn(m,n);
y = abs(A*x).^2;

% initialization
npower_iter = 50;                           % number of power iterations 
z0 = randn(n,1); z0 = z0/norm(z0,'fro');    % initial guess 
for tt = 1:npower_iter,                     % power iterations
    z0 = A'*(y.* (A*z0)); z0 = z0/norm(z0,'fro');
end

normest = sqrt(sum(y)/numel(y));      
z = normest * z0;                   


% Kaczmarz
%maxiter = 100;
zk = z;
relerrKac = norm(x - exp(-1i*angle(trace(x'*zk))) * zk, 'fro')/norm(x,'fro');

for t = 1 : maxiter
  for r = 1 : m
    nrm2 = norm(A(r,:))^2;
    in = A(r,:)*zk;
    
    zk = zk + (in/abs(in)*sqrt(y(r))-in)*A(r,:)'/nrm2;
  end
  relerrKac = [relerrKac,norm(x - exp(-1i*angle(trace(x'*zk))) * zk, 'fro')/norm(x,'fro')];
end

relerrKac_array(:,ll) = relerrKac;

%relaxed Kaczmarz
%maxiter = 100;
zrk = z;
lambda = 1 + n/m;
relerrRKac = norm(x - exp(-1i*angle(trace(x'*zrk))) * zrk, 'fro')/norm(x,'fro');

for t = 1 : maxiter
  for r = 1 : m
    nrm2 = norm(A(r,:))^2;
    in = A(r,:)*zrk;
    
    zrk = zrk + lambda*(in/abs(in)*sqrt(y(r))-in)*A(r,:)'/nrm2; 
  end
  relerrRKac = [relerrRKac,norm(x - exp(-1i*angle(trace(x'*zrk))) * zrk, 'fro')/norm(x,'fro')];
end

relerrRKac_array(:,ll) = relerrRKac;

%{
fprintf('Kacz: relative error = %g\n', relerrKac(end))
fprintf('Relaxed Kacz: relative error = %g\n', relerrRKac(end))

semilogy(0:maxiter,relerrKac,'r') 
hold on
semilogy(0:maxiter,relerrRKac,'k') 

legend('Kaczmarz', 'Relaxed Kaczmarz')

xlabel('Iteration')
ylabel('Relative error (log10)')
title('Relative error vs. iteration count')
%}
end

save('relerrKac_array.mat','relerrKac_array')
save('relerrRKac_array.mat','relerrRKac_array')
