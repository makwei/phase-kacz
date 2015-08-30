% process gaussian data
function process_cdp_gal_data(fname)

comm_string = [' sed ''s/,//g'' ' fname '> output1.txt'];
system(comm_string);
comm_string = ' sed ''s/Block//g'' output1.txt > output.txt';
system(comm_string);

fid = fopen('output.txt');
%GerSax -> cdp2d, eig_init, n1: 1080, n2: 1920, m: 24883200, mode: 1, iter: 254, t: 859.2, relres: 8.07278e-11, relerr: 1.24255e-10
tmp = textscan(fid,'%s %s %s %s %s %d %s %d %s %d %s %d %s %d %s %f %s %f %s %f');
fclose(fid);
delete('output1.txt');
delete('output.txt');

results = zeros(length(tmp{6}),8);
results(:,1) = tmp{6};  % n1
results(:,2) = tmp{8};  % n2
results(:,3) = tmp{10}; % m
results(:,4) = tmp{12}; % mode
results(:,5) = tmp{14}; % iter
results(:,6) = tmp{16}; % time
results(:,7) = tmp{18}; % res
results(:,8) = tmp{20}; % err
clear tmp;

if mod(size(results,1),3)~=0
  disp('data number not right!')
end

num_rows = size(results,1)/3;
results_tmp = zeros(num_rows,8);
for ll = 1 : num_rows
  results_tmp(ll,:) = mean(results(3*(ll-1)+1:3*ll,:));
end

results = results_tmp;
clear results_tmp;

results(:,5) = round(results(:,5));

n1_n2_m = intersect(results(:,1:3),results(:,1:3),'rows');
len = size(n1_n2_m,1);
succ = zeros(len,1);
ave_iter = zeros(len,1);
ave_time = zeros(len,1);
ave_res = zeros(len,1);
ave_err = zeros(len,1);

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
    fprintf('n: %d, m: %d num of tets: %d\n', n, m, length(ind));
  end
  
  
  ind = ind(results(ind,8)<1e-5);
  succ(i) = length(ind);
  ave_iter(i) = round(mean(results(ind,5)));
  ave_time(i) = mean(results(ind,6));
  ave_res(i) = mean(results(ind,7));
  ave_err(i) = mean(results(ind,8));
end

n1_n2_m_all = [n1_n2_m succ ave_iter ave_time ave_res ave_err];

fname(end-2:end) = 'mat';
save(fname,'results','n1_n2_m_all')