% do robust to noise plots
labsz = 18;
linewd = 1.5;
titsz = 18;
legsz = 18;
marksz = 10;

% ---------------------
% Gaussian real
% ---------------------
clist = 'bkrg';
fnames = {'GerSax_guassian_real_eig.mat', 'wirtinger_guassian_real_eig.mat',...
    'kaczmarz_guassian_real_eig.mat','block_kaczmarz_guassian_real_64_eig.mat'};
folder = '../Noise/';
figure(1); clf; hold on
for i = 1 : 4
  load([folder fnames{i}])
  plot(n_m_eps_all(:,3),n_m_eps_all(:,8),clist(i),'linewidth',linewd)
end
set(gca,'xscale','log')
set(gca,'yscale','log')
h = legend('ER','Wirtinger Flow','Kaczmarz', 'Kaczmarz (64)','location','southeast');
set(h,'fontsize',legsz)
xlabel('Noise level','fontsize',labsz)
ylabel('Relative error','fontsize',labsz)
title('Gaussian, real','fontsize',titsz)
%axis([0, 0.11 0 0.3])
box on

print ../Plots/gaussian_real_er_wirtinger_kacz_noise.pdf -dpdf

% ---------------------
% Gaussian complex
% ---------------------
clist = 'bkrg';
fnames = {'GerSax_guassian_complex_eig.mat', 'wirtinger_guassian_complex_eig.mat',...
    'kaczmarz_guassian_complex_eig.mat','block_kaczmarz_guassian_complex_64_eig.mat'};
folder = '../Noise/';
figure(2); clf; hold on
for i = 1 : 4
  load([folder fnames{i}])
  plot(n_m_eps_all(:,3),n_m_eps_all(:,8),clist(i),'linewidth',linewd)
end
set(gca,'xscale','log')
set(gca,'yscale','log')
h = legend('ER','Wirtinger Flow','Kaczmarz', 'Kaczmarz (64)','location','southeast');
set(h,'fontsize',legsz)
xlabel('Noise level','fontsize',labsz)
ylabel('Relative error','fontsize',labsz)
title('Gaussian, complex','fontsize',titsz)
%axis([0, 0.11 0 0.3])
box on

print ../Plots/gaussian_complex_er_wirtinger_kacz_noise.pdf -dpdf

% ---------------------
% CDP 1D
% ---------------------
clist = 'bkrg';
fnames = {'GerSax_cdp_eig.mat', 'wirtinger_cdp_eig.mat',...
    'kaczmarz_cdp_eig.mat','block_kaczmarz_cdp_256_eig.mat'};
folder = '../Noise/';
figure(3); clf; hold on
for i = 1 : 4
  load([folder fnames{i}])
  plot(n_m_eps_all(:,3),n_m_eps_all(:,8),clist(i),'linewidth',linewd)
end
set(gca,'xscale','log')
set(gca,'yscale','log')
h = legend('ER','Wirtinger Flow','Kaczmarz', 'Kaczmarz (256)','location','southeast');
set(h,'fontsize',legsz)
xlabel('Noise level','fontsize',labsz)
ylabel('Relative error','fontsize',labsz)
title('CDP (1D), complex','fontsize',titsz)
%axis([0, 0.11 0 0.3])
box on

print ../Plots/cdp1d_er_wirtinger_kacz_noise.pdf -dpdf

% ----------------------
% CDP 2D
% ----------------------
clist = 'bkg';
fnames = {'GerSax_cdp2d_eig.mat', 'wirtinger_cdp2d_eig.mat','block_kaczmarz_cdp2d_eig.mat'};
folder = '../Noise/';
figure(4); clf; hold on
for i = 1 : 3
  load([folder fnames{i}])
  plot(n1_n2_m_eps_all(:,4),n1_n2_m_eps_all(:,9),clist(i),'linewidth',linewd)
end
set(gca,'xscale','log')
set(gca,'yscale','log')
h = legend('ER','Wirtinger Flow', 'Kaczmarz (256\times256)','location','southeast');
set(h,'fontsize',legsz)
xlabel('Noise level','fontsize',labsz)
ylabel('Relative error','fontsize',labsz)
title('CDP (2D), complex','fontsize',titsz)
%axis([0, 0.11 0 0.3])
box on

print ../Plots/cdp2d_er_wirtinger_kacz_noise.pdf -dpdf

