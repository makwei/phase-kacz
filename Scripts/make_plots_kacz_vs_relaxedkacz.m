% Kaczmarz vs Relaxed Kaczmarz
labsz = 18;
linewd = 1.5;
titsz = 18;
legsz = 18;

load('relerrKac_array.mat')
load('relerrRKac_array.mat')

ind1 = find(relerrKac_array(end,:)<1e-12);
ind2 = find(relerrRKac_array(end,:)<1e-12);
ind = intersect(ind1,ind2);

relerrKac_array = relerrKac_array(:,ind);
relerrRKac_array = relerrRKac_array(:,ind);

figure(1)
clf;
semilogy(mean(relerrKac_array,2),'r','linewidth',linewd)
hold on 
semilogy(mean(relerrRKac_array,2),'k','linewidth',linewd)

h=legend('simple Kaczmarz','simple Kaczmarz with relaxation');
set(h,'fontsize',legsz)
%title('Convergence rate comparison','fontsize',titsz)
xlabel('Cycle','fontsize',labsz)
ylabel('Relative error','fontsize',labsz)
box on

print ../Plots/kacz_vs_relaxedkacz.pdf -dpdf