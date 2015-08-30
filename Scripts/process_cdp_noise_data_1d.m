% process gaussian data
function process_cdp_noise_data_1d(fname)

comm_string = [' sed ''s/,//g'' ' fname '> output1.txt'];
system(comm_string);
comm_string = ' sed ''s/Block//g'' output1.txt > output.txt';
system(comm_string);

fid = fopen('output.txt');
%GerSax -> cdp, eig_init, n: 256, m: 1536, epsilon: 0.1, iter: 2500, t: 0.932974, relres: 0.102541, relerr: 0.117889
tmp = textscan(fid,'%s %s %s %s %s %d %s %d %s %f %s %d %s %f %s %f %s %f');
fclose(fid);
delete('output1.txt');
delete('output.txt');

results = zeros(length(tmp{6}),7);
results(:,1) = tmp{6};  % n
results(:,2) = tmp{8};  % m
results(:,3) = tmp{10}; % epsilon
results(:,4) = tmp{12}; % iter
results(:,5) = tmp{14}; % time
results(:,6) = tmp{16}; % res
results(:,7) = tmp{18}; % err
clear tmp;


n_m_eps = intersect(results(:,1:3),results(:,1:3),'rows');
len = size(n_m_eps,1);
succ = zeros(len,1);
ave_iter = zeros(len,1);
ave_time = zeros(len,1);
ave_res = zeros(len,1);
ave_err = zeros(len,1);

num_tests = 10;
for i = 1 : size(n_m_eps,1)
  n = n_m_eps(i,1);
  m = n_m_eps(i,2);
  epsilon = n_m_eps(i,3);
  ind_n = find(results(:,1)==n);
  ind_m = find(results(:,2)==m);
  ind_eps = find(results(:,3)==epsilon);
  ind = intersect(intersect(ind_n,ind_m),ind_eps);
  
  if length(ind) < num_tests
    fprintf('n: %d, m: %d num of tets: %d\n', n, m, length(ind));
  end
  
  ind = ind(~isnan(results(ind,7)));
  succ(i) = sum(ind);
  ave_iter(i) = round(mean(results(ind,4)));
  ave_time(i) = mean(results(ind,5));
  ave_res(i) = mean(results(ind,6));
  ave_err(i) = mean(results(ind,7));
end

n_m_eps_all = [n_m_eps succ ave_iter ave_time ave_res ave_err];



fname(end-2:end) = 'mat';
save(fname,'results','n_m_eps_all');