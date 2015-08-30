function process_cdp_data_2d(fname)

comm_string = [' sed ''s/,//g'' ' fname '> output1.txt'];
system(comm_string);
comm_string = ' sed ''s/Block//g'' output1.txt > output.txt';
system(comm_string);

fid = fopen('output.txt');
% GerSax -> cdp2d, eig_init, n1: 256, n2: 256, m: 131072, iter: 2500, t: 38.7723, relres: 0.0256928, relerr: 1.53303
tmp = textscan(fid,'%s %s %s %s %s %d %s %d %s %d %s %d %s %f %s %f %s %f');
fclose(fid);
delete('output1.txt');
delete('output.txt');

results = zeros(length(tmp{6}),7);
results(:,1) = tmp{6};  % n1
results(:,2) = tmp{8};  % n2
results(:,3) = tmp{10}; % m
results(:,4) = tmp{12}; % iter
results(:,5) = tmp{14}; % time
results(:,6) = tmp{16}; % res
results(:,7) = tmp{18}; % err
clear tmp;


n1_n2_m = intersect(results(:,1:3),results(:,1:3),'rows');
len = size(n1_n2_m,1);
succ = zeros(len,1);
ave_iter = zeros(len,1);
ave_time = zeros(len,1);
ave_res = zeros(len,1);
ave_err = zeros(len,1);

tol = 1e-5;
num_tests = 10;
for i = 1 : size(n1_n2_m,1)
  n1 = n1_n2_m(i,1);
  n2 = n1_n2_m(i,2);
  m = n1_n2_m(i,3);
  
  ind_n1 = find(results(:,1)==n1);
  ind_n2 = find(results(:,2)==n2);
  ind_m = find(results(:,3)==m);
  ind = intersect(intersect(ind_n1,ind_n2),ind_m);
  
  if length(ind) < num_tests
    fprintf('n1: %d, n2: %d m: %d num of tets: %d\n', n1, n2, m, length(ind));
  end
  
  succ(i) = sum(results(ind,7)<=tol);
  ave_iter(i) = round(sum(results(ind,4).*((results(ind,7)<=tol)&~isnan(results(ind,4))))/succ(i));
  ave_time(i) = sum(results(ind,5).*((results(ind,7)<=tol)&~isnan(results(ind,5))))/succ(i);
  ave_res(i) = sum(results(ind,6).*((results(ind,7)<=tol)&~isnan(results(ind,6))))/succ(i);
  ave_err(i) = sum(results(ind,7).*((results(ind,7)<=tol)&~isnan(results(ind,7))))/succ(i);
end

n1_n2_m_all = [n1_n2_m succ ave_iter ave_time ave_res ave_err];

fname(end-2:end) = 'mat';
save(fname,'results','n1_n2_m_all');