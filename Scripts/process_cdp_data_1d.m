% process data 1d
function process_cdp_data_1d(fname)

comm_string = [' sed ''s/,//g'' ' fname '> output1.txt'];
system(comm_string);
comm_string = ' sed ''s/Block//g'' output1.txt > output.txt';
system(comm_string);

fid = fopen('output.txt');
% GerSax -> cdp, eig_init, n: 256, m: 512, iter: 2500, t: 0.462963, relres: 0.0212544, relerr: 1.25027
tmp = textscan(fid,'%s %s %s %s %s %d %s %d %s %d %s %f %s %f %s %f');
fclose(fid);
delete('output1.txt');
delete('output.txt');

results = zeros(length(tmp{6}),6);
results(:,1) = tmp{6};  % n
results(:,2) = tmp{8};  % m
results(:,3) = tmp{10}; % iter
results(:,4) = tmp{12}; % time
results(:,5) = tmp{14}; % res
results(:,6) = tmp{16}; % err
clear tmp;


n_m = intersect(results(:,1:2),results(:,1:2),'rows');
len = size(n_m,1);
succ = zeros(len,1);
ave_iter = zeros(len,1);
ave_time = zeros(len,1);
ave_res = zeros(len,1);
ave_err = zeros(len,1);

tol = 1e-5;
num_tests = 10;
for i = 1 : size(n_m,1)
  n = n_m(i,1);
  m = n_m(i,2);
  ind_n = find(results(:,1)==n);
  ind_m = find(results(:,2)==m);
  ind = intersect(ind_n,ind_m);
  
  if length(ind) < num_tests
    fprintf('n: %d, m: %d num of tets: %d\n', n, m, length(ind));
  end
  
  succ(i) = sum(results(ind,6)<=tol);
  ave_iter(i) = round(sum(results(ind,3).*((results(ind,6)<=tol)&~isnan(results(ind,3))))/succ(i));
  ave_time(i) = sum(results(ind,4).*((results(ind,6)<=tol)&~isnan(results(ind,4))))/succ(i);
  ave_res(i) = sum(results(ind,5).*((results(ind,6)<=tol)&~isnan(results(ind,5))))/succ(i);
  ave_err(i) = sum(results(ind,6).*((results(ind,6)<=tol)&~isnan(results(ind,6))))/succ(i);
end

n_m_all = [n_m succ ave_iter ave_time ave_res ave_err];

fname(end-2:end) = 'mat';
save(fname,'results','n_m_all');