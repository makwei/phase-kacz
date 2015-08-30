num_tests = 50;
clist = 'brymcg';
marker = '+o*sdp';
labsz = 18;
linewd = 1.5;
titsz = 18;
legsz = 16;
marksz = 10;

% ---------------------------
% Gaussian: Real
% ---------------------------
fnames1 = {'GerSax_guassian_real_eig.mat','kaczmarz_guassian_real_eig.mat','block_kaczmarz_guassian_real_32_eig.mat',...
    'block_kaczmarz_guassian_real_64_eig.mat','block_kaczmarz_guassian_real_128_eig.mat'};
fnames2 = {'GerSax_guassian_real_rand.mat','kaczmarz_guassian_real_rand.mat','block_kaczmarz_guassian_real_32_rand.mat',...
    'block_kaczmarz_guassian_real_64_rand.mat','block_kaczmarz_guassian_real_128_rand.mat'};
folder = '../Gaussian/';
figure(1), clf, hold on
for i = 1 : 5
  load([folder fnames1{i}])
  if i ~= 3
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,clist(i),'linewidth',linewd)
  else
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'color',[1 0.6 0],'linewidth',linewd)
  end
  
  load([folder fnames2{i}])
  if i ~= 3
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'linestyle','--','color',clist(i),'linewidth',linewd)
  else
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'linestyle','--','color',[1 0.6 0],'linewidth',linewd)
  end
end
h=legend('ER (eig)','ER (rand)','Kaczmarz (eig)', 'Kaczmarz (rand)','Kaczmarz (32, eig)','Kaczmarz (32, rand)',...
    'Kaczmarz (64, eig)', 'Kaczmarz (64, rand)', 'Kaczmarz (128, eig)', 'Kaczmarz (128, rand)','location','southeast');
set(h,'fontsize',legsz)
xlabel('m/n','fontsize',labsz)
ylabel('Successful recovery probability','fontsize',labsz)
title('Gaussian, real','fontsize',titsz)
axis([1.9 6.1 -0.05 1.05])
box on

print ../Plots/gaussian_real_er_kacz_eig_vs_rand.pdf -dpdf

% ---------------------------
% Gaussian: Complex
% ---------------------------
fnames1 = {'GerSax_guassian_complex_eig.mat','kaczmarz_guassian_complex_eig.mat','block_kaczmarz_guassian_complex_32_eig.mat',...
    'block_kaczmarz_guassian_complex_64_eig.mat','block_kaczmarz_guassian_complex_128_eig.mat'};
fnames2 = {'GerSax_guassian_complex_rand.mat','kaczmarz_guassian_complex_rand.mat','block_kaczmarz_guassian_complex_32_rand.mat',...
    'block_kaczmarz_guassian_complex_64_rand.mat','block_kaczmarz_guassian_complex_128_rand.mat'};
folder = '../Gaussian/';
figure(2), clf, hold on
for i = 1 : 5
  load([folder fnames1{i}])
  if i ~= 3
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,clist(i),'linewidth',linewd)
  else
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'color',[1 0.6 0],'linewidth',linewd)
  end
  
  load([folder fnames2{i}])
  if i ~= 3
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'linestyle','--','color',clist(i),'linewidth',linewd)
  else
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'linestyle','--','color',[1 0.6 0],'linewidth',linewd)
  end
end
h=legend('ER (eig)','ER (rand)','Kaczmarz (eig)', 'Kaczmarz (rand)','Kaczmarz (32, eig)','Kaczmarz (32, rand)',...
    'Kaczmarz (64, eig)', 'Kaczmarz (64, rand)', 'Kaczmarz (128, eig)', 'Kaczmarz (128, rand)','location','southeast');
set(h,'fontsize',legsz)
xlabel('m/n','fontsize',labsz)
ylabel('Successful recovery probability','fontsize',labsz)
title('Gaussian, complex','fontsize',titsz)
axis([1.9 6.1 -0.05 1.05])
box on

print ../Plots/gaussian_complex_er_kacz_eig_vs_rand.pdf -dpdf

% ---------------------------
% Unitary: Real
% ---------------------------
fnames1 = {'GerSax_unitary_real_eig.mat','block_kaczmarz_unitary_real_256_eig.mat'};
fnames2 = {'GerSax_unitary_real_rand.mat','block_kaczmarz_unitary_real_256_rand.mat'};
folder = '../Unitary/';
figure(3), clf, hold on
for i = 1 : 2
  load([folder fnames1{i}])
  if i ~= 3
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,clist(i),'linewidth',linewd)
  else
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'color',[1 0.6 0],'linewidth',linewd)
  end
  
  load([folder fnames2{i}])
  if i ~= 3
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'linestyle','--','color',clist(i),'linewidth',linewd)
  else
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'linestyle','--','color',[1 0.6 0],'linewidth',linewd)
  end
end
h=legend('ER (eig)','ER (rand)','Kaczmarz (eig)', 'Kaczmarz (rand)','location','southeast');
set(h,'fontsize',legsz)
xlabel('m/n','fontsize',labsz)
ylabel('Successful recovery probability','fontsize',labsz)
title('Unitary, real','fontsize',titsz)
axis([1.9 6.1 -0.05 1.05])
box on

print ../Plots/unitary_real_er_kacz_eig_vs_rand.pdf -dpdf

% ---------------------------
% Unitary: Complex
% ---------------------------
fnames1 = {'GerSax_unitary_complex_eig.mat','block_kaczmarz_unitary_complex_256_eig.mat'};
fnames2 = {'GerSax_unitary_complex_rand.mat','block_kaczmarz_unitary_complex_256_rand.mat'};
folder = '../Unitary/';
figure(4), clf, hold on
for i = 1 : 2
  load([folder fnames1{i}])
  if i ~= 3
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,clist(i),'linewidth',linewd)
  else
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'color',[1 0.6 0],'linewidth',linewd)
  end
  
  load([folder fnames2{i}])
  if i ~= 3
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'linestyle','--','color',clist(i),'linewidth',linewd)
  else
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'linestyle','--','color',[1 0.6 0],'linewidth',linewd)
  end
end
h=legend('ER (eig)','ER (rand)','Kaczmarz (eig)', 'Kaczmarz (rand)','location','southeast');
set(h,'fontsize',legsz)
xlabel('m/n','fontsize',labsz)
ylabel('Successful recovery probability','fontsize',labsz)
title('Unitary, complex','fontsize',titsz)
axis([1.9 6.1 -0.05 1.05])
box on

print ../Plots/unitary_complex_er_kacz_eig_vs_rand.pdf -dpdf

% ------------------------------
% CDP 1d
% ------------------------------
fnames1 = {'GerSax_cdp_eig.mat','kaczmarz_cdp_eig.mat','block_kaczmarz_cdp_32_eig.mat',...
    'block_kaczmarz_cdp_64_eig.mat','block_kaczmarz_cdp_128_eig.mat','block_kaczmarz_cdp_256_eig.mat'};
fnames2 = {'GerSax_cdp_rand.mat','kaczmarz_cdp_rand.mat','block_kaczmarz_cdp_32_rand.mat',...
    'block_kaczmarz_cdp_64_rand.mat','block_kaczmarz_cdp_128_rand.mat','block_kaczmarz_cdp_256_rand.mat'};
folder = '../CDP1d/';
figure(5), clf, hold on
for i = 1 : 6
  load([folder fnames1{i}])
  if i ~= 3
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,clist(i),'linewidth',linewd,'marker',marker(i),'markersize',marksz)
  else
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'color',[1 0.6 0],'linewidth',linewd,'marker',marker(i),'markersize',marksz)
  end
  
  load([folder fnames2{i}])
  if i ~= 3
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'linestyle','--','color',clist(i),'linewidth',linewd,'marker',marker(i),'markersize',marksz)
  else
    plot(n_m_all(:,2)./n_m_all(:,1),n_m_all(:,3)/num_tests,'linestyle','--','color',[1 0.6 0],'linewidth',linewd,'marker',marker(i),'markersize',marksz)
  end
end
h=legend('ER (eig)','ER (rand)','Kaczmarz (eig)', 'Kaczmarz (rand)','Kaczmarz (32, eig)','Kaczmarz (32, rand)',...
    'Kaczmarz (64, eig)', 'Kaczmarz (64, rand)', 'Kaczmarz (128, eig)', 'Kaczmarz (128, rand)','Kaczmarz (256, eig)', 'Kaczmarz (256, rand)','location','southeast');
set(h,'fontsize',legsz)
xlabel('m/n','fontsize',labsz)
ylabel('Successful recovery probability','fontsize',labsz)
title('CDP (1D), complex','fontsize',titsz)
axis([1.9 12.1 -0.05 1.05])
box on

print ../Plots/cdp1d_er_kacz_eig_vs_rand.pdf -dpdf

% ----------------------------
% CDP 2d
% ----------------------------
fnames1 = {'GerSax_cdp2d_eig.mat', 'block_kaczmarz_cdp2d_eig.mat'};
fnames2 = {'GerSax_cdp2d_rand.mat','block_kaczmarz_cdp2d_rand.mat'};
folder = '../CDP2d/';
figure(6), clf, hold on
for i = 1 : 2
  load([folder fnames1{i}])
  if i ~= 2
    plot(n1_n2_m_all(:,3)./(n1_n2_m_all(:,1).*n1_n2_m_all(:,2)),n1_n2_m_all(:,4)/num_tests,clist(i),'linewidth',linewd,'marker',marker(i),'markersize',marksz)
  else
    plot(n1_n2_m_all(:,3)./(n1_n2_m_all(:,1).*n1_n2_m_all(:,2)),n1_n2_m_all(:,4)/num_tests,clist(i+4),'linewidth',linewd,'marker',marker(i+4),'markersize',marksz)
  end
  
  load([folder fnames2{i}])
  if i ~= 2
    plot(n1_n2_m_all(:,3)./(n1_n2_m_all(:,1).*n1_n2_m_all(:,2)),n1_n2_m_all(:,4)/num_tests,'linestyle','--','color',clist(i),'linewidth',linewd,'marker',marker(i),'markersize',marksz)
  else
    plot(n1_n2_m_all(:,3)./(n1_n2_m_all(:,1).*n1_n2_m_all(:,2)),n1_n2_m_all(:,4)/num_tests,'linestyle','--','color',clist(i+4),'linewidth',linewd,'marker',marker(i+4),'markersize',marksz)
  end
end
h=legend('ER (eig)','ER (rand)', 'Kaczmarz (256\times256, eig)', 'Kaczmarz (256\times256, rand)','location','southeast');
set(h,'fontsize',legsz)
xlabel('m/n_1n_2','fontsize',labsz)
ylabel('Successful recovery probability','fontsize',labsz)
title('CDP (2D), complex','fontsize',titsz)
axis([1.9 12.1 -0.05 1.05])
box on

print ../Plots/cdp2d_er_kacz_eig_vs_rand.pdf -dpdf