clear all
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));

% problem size
n = 256;
m = 6*n;
noise_level = [1e-1,0.5e-1,1e-2,0.5e-2,1e-3,0.5e-3,1e-4,0.5e-4,1e-5];
type = {'real','complex'};

num_tests = 10;

for tt = 1 : length(type)
  for ll = 1 : length(noise_level)
  epsilon = noise_level(ll);
  for num = 1 : num_tests
  %----------------------------
    % generate data
    if strcmp(type{tt},'real') == 1 % real
      x = randn(n,1);
      A = randn(m,n);   
      y = abs(A*x).^2;
      
      % generate additive noise
      e = randn(m,1);
      e = epsilon*norm(y)*e/norm(e);
      
      y = max(y+e,0);
    else                            % complex
      x = randn(n,1) + 1i*randn(n,1);
      A = 1/sqrt(2)*randn(m,n)+1i/sqrt(2)*randn(m,n);
      y = abs(A*x).^2 ;
      
      % generate additive noise
      e = randn(m,1);
      e = epsilon*norm(y)*e/norm(e);
      y = max(y+e,0);
    end
    
    % initilization 
    z0 = randn(n,1);
    z0 = z0/norm(z0,'fro');
    for ii = 1 : 50
      z0 = A'*(y.* (A*z0));
      z0 =  z0/norm(z0,'fro');
    end
    normest = sqrt(sum(y)/numel(y));
    z = normest * z0;
    
    %--- Wirtinger flow ---%
    maxiter = 2500;
    tol = epsilon;
    tau0 = 330;
    mu = @(t) min(1-exp(-t/tau0), 0.2);
    zw = z;
    
    tic;
    for iter = 1 : maxiter
      Azw = A*zw;
      
      diff = abs(Azw).^2-y;
      if norm(diff)/norm(y) <= tol 
        break;
      end
      
      grad  = 1/m* A'*( diff .* Azw );
      zw = zw - mu(iter)/normest^2 * grad; 
    end
    t = toc;
    
    % output results of wirtinger flow
    fname = ['wirtinger_gaussian_' type{tt} '_eig.txt'];
    fid = fopen(fname,'a');
    fprintf(fid, ['Wirtinger -> gaussian, ' type{tt} ', eig_init, n: %d, m: %d, epsilon: %g, iter: %d, t: %g, relres: %g, relerr: %g\n'], n, m, epsilon, iter, t, norm(abs(Azw).^2-y)/norm(y), ... 
        norm(x - exp(-1i*angle(trace(x'*zw))) * zw, 'fro')/norm(x,'fro'));
    fclose(fid);
    
    %--- Gerchberg-Saxton ---%
    maxiter = 2500;
    tol = epsilon;
    ATA = A'*A;
    zg = z;
    
    tic;
    for iter = 1 : maxiter
      Azg= A*zg;
      
      if norm(abs(Azg).^2-y)/norm(y) <= tol
        break;
      end
      
      yzg = sqrt(y).*(Azg./abs(Azg));
      %Azg = ATA\(A'*yzg);
      [zg,~] = pcg(ATA,A'*yzg,1e-6,n,[],[],zg);
    end
    t = toc;

    % output results of Gerchberg-Saxton
    fname = ['GerSax_gaussian_' type{tt} '_eig.txt'];
    fid = fopen(fname,'a');
    fprintf(fid, ['GerSax -> gaussian, ' type{tt} ', eig_init, n: %d, m: %d, epsilon: %g, iter: %d, t: %g, relres: %g, relerr: %g\n'], n, m, epsilon, iter, t, norm(abs(Azg).^2-y)/norm(y), ... 
        norm(x - exp(-1i*angle(trace(x'*zg))) * zg, 'fro')/norm(x,'fro'));
    fclose(fid);
    

    %--- Kaczmarz method ---%
    maxiter = 500;
    tol = epsilon;
    zk = z;
    
    tic;
    for iter = 1 : maxiter
      maxrelres = 0;
      for r = 1 : m
        nrm2 = norm(A(r,:))^2;
        Arzk = A(r,:)*zk;
        
        maxrelres = max(maxrelres,abs(y(r)-abs(Arzk)^2)/y(r)); 
        
        zk = zk + (Arzk/abs(Arzk)*sqrt(y(r))-Arzk)*A(r,:)'/nrm2;
      end
      
      if maxrelres < tol
        break;
      end
    end
    t = toc;
    
    % output results of Kaczmarz method
    fname = ['kaczmarz_guassian_' type{tt} '_eig.txt'];
    fid = fopen(fname,'a');
    fprintf(fid, ['Kaczmarz -> gaussian, ' type{tt} ', eig_init, n: %d, m: %d, epsilon: %g, iter: %d, t: %g, relres: %g, relerr: %g\n'], n, m, epsilon, iter, t, norm(y-abs(A*zk).^2)/norm(y), ... 
        norm(x - exp(-1i*angle(trace(x'*zk))) * zk, 'fro')/norm(x,'fro'));
    fclose(fid);
   
    
    %--- Block Kaczmarz ---%
    maxiter = 500;
    tol = epsilon;
    blockSS = round(n/4);
    
    for bls = 1 : length(blockSS)
      blockSize = blockSS(bls);
      nblocks = ceil(m/blockSize);
      zbk = z;
      
      tic; 
      for iter = 1 : maxiter
        maxrelres = 0;
        for r = 1 : nblocks
          Asub = A((r-1)*blockSize+1:min(r*blockSize,m),:);
          ysub = y((r-1)*blockSize+1:min(r*blockSize,m));
          Asubzbk = Asub*zbk;
          
          maxrelres = max(maxrelres,norm(ysub-abs(Asubzbk).^2)/norm(ysub));
          
          zbk = zbk +pinv(Asub)*((Asubzbk./abs(Asubzbk)).*sqrt(ysub)-Asubzbk);
        end
        
        if maxrelres < tol
           break;
        end
      end
      t = toc;
      
      % output results for block kaczmarz
      fname = ['block_kaczmarz_guassian_' type{tt} '_' num2str(blockSize) '_eig.txt'];
      fid = fopen(fname,'a');
      alg_name = ['Block Kaczmarz-' num2str(blockSize)];
      fprintf(fid, [alg_name ' -> guassian, ' type{tt} ', eig_init, n: %d, m: %d, epsilon: %g, iter: %d, t: %g, relres: %g, relerr: %g\n'], n, m, epsilon,iter, t, norm(y-abs(A*zbk).^2)/norm(y), ... 
        norm(x - exp(-1i*angle(trace(x'*zbk))) * zbk, 'fro')/norm(x,'fro'));
      fclose(fid);
    end
  %--------------------------
  end
  end
end


