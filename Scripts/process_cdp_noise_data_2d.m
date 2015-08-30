% process gaussian data
function process_cdp_noise_data_2d(fname)

comm_string = [' sed ''s/,//g'' ' fname '> output1.txt'];
system(comm_string);
comm_string = ' sed ''s/Block//g'' output1.txt > output.txt';
system(comm_string);

fid = fopen('output.txt');
%GerSax -> cdp2d, eig_init, n1: 256, n2: 256, m: 786432, epsilon: 0.001, iter: 2500, t: 242.236, relres: 0.00131425, relerr: 0.00122092
tmp = textscan(fid,'%s %s %s %s %s %d %s %d %s %d %s %f %s %d %s %f %s %f %s %f');
fclose(fid);
delete('output1.txt');
delete('output.txt');

results = zeros(length(tmp{6}),8);
results(:,1) = tmp{6};  % n1
results(:,2) = tmp{8};  % n2
results(:,3) = tmp{10}; % m
results(:,4) = tmp{12}; % epsilon
results(:,5) = tmp{14}; % iter
results(:,6) = tmp{16}; % time
results(:,7) = tmp{18}; % res
results(:,8) = tmp{20}; % err
clear tmp;


n1_n2_m_eps= intersect(results(:,1:4),results(:,1:4),'rows');
len = size(n1_n2_m_eps,1);
succ = zeros(len,1);
ave_iter = zeros(len,1);
ave_time = zeros(len,1);
ave_res = zeros(len,1);
ave_err = zeros(len,1);

num_tests = 10;
for i = 1 : size(n1_n2_m_eps,1)
  n1 = n1_n2_m_eps(i,1);
  n2 = n1_n2_m_eps(i,2);
  m = n1_n2_m_eps(i,3);
  epsilon = n1_n2_m_eps(i,4);
  
  ind_n1 = find(results(:,1)==n1);
  ind_n2 = find(results(:,2)==n2);
  ind_m = find(results(:,3)==m);
  ind_eps = find(results(:,4)==epsilon);
  ind = intersect(intersect(intersect(ind_n1,ind_n2),ind_m),ind_eps);
  
  if length(ind) < num_tests
    fprintf('n: %d, m: %d num of tets: %d\n', n, m, length(ind));
  end
  
  
  ind = ind(~isnan(results(ind,8)));
  succ(i) = sum(ind);
  ave_iter(i) = round(mean(results(ind,5)));
  ave_time(i) = mean(results(ind,6));
  ave_res(i) = mean(results(ind,7));
  ave_err(i) = mean(results(ind,8));
end

n1_n2_m_eps_all = [n1_n2_m_eps succ ave_iter ave_time ave_res ave_err];



fname(end-2:end) = 'mat';
save(fname,'results','n1_n2_m_eps_all');