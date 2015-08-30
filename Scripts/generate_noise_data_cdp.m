clear all
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));

% problem size
n = 256;
L = 6; 
noise_level = [1e-1,0.5e-1,1e-2,0.5e-2,1e-3,0.5e-3,1e-4,0.5e-4,1e-5];

num_tests = 10;

dict = [1i -1i 1 -1];

for ll = 1 : length(noise_level)
  epsilon = noise_level(ll);
  for num = 1 : num_tests
    % generate the mask
    mask = zeros(n,L);
    for i = 1 : n
      for j = 1 : L
        mask(i,j) = dict(randi(4));
      end
    end
    
    temp = rand(size(mask));
    mask = mask .* ( (temp <= 0.2)*sqrt(3) + (temp > 0.2)/sqrt(2) );
    
    % define the operators
    A = @(I)  fft(conj(mask) .* repmat(I,[1 L]));  % Input is n x 1 signal, output is n x L array
    At = @(Y) mean(mask .* ifft(Y), 2);            % Input is n x L array, output is n x 1 signal          
    Ainv = @(Y) sum(mask.*ifft(Y),2)./sum(abs(mask).^2,2);
    
    % generate data
    x = randn(n,1) + 1i*randn(n,1);
    Y = abs(A(x)).^2; 
    
    % generate noise 
    e = randn(size(Y));
    e = epsilon*norm(Y(:))*e/norm(e(:));
    
    % noisy data
    Y = max(Y+e,0);
    
    % for the simple Kaczmarz 
    Amat = [];
    for i = 1 : L
      Amat = [Amat; fft(diag(conj(mask(:,i))))];
    end
    yv = Y(:);
    
    % initilization
    npower_iter = 50;                           
    z0 = randn(n,1); z0 = z0/norm(z0,'fro'); 
    for tt = 1:npower_iter, 
      z0 = At(Y.*A(z0)); z0 = z0/norm(z0,'fro');
    end
    normest = sqrt(sum(Y(:))/numel(Y)); % Estimate norm to scale eigenvector  
    z = normest * z0;   
   
    % --- Wirtinger flow --- %
    maxiter = 2500;          
    tol = epsilon;
    tau0 = 330;                         
    mu = @(t) min(1-exp(-t/tau0), 0.2);    
    zw = z;
    
    tic;
    for iter = 1 : maxiter
      Azw = A(zw);
      
      diff = abs(Azw).^2-Y;
      if norm(diff,'fro')/norm(Y,'fro') <= tol
        break;
      end
      
      C  = diff .* Azw;
      grad = At(C);
      zw = zw - mu(iter)/normest^2 * grad;
    end
    t = toc;
    
    % output results of wirtinger flow
    fname = 'wirtinger_cdp_eig.txt';
    fid = fopen(fname,'a');
    fprintf(fid, 'Wirtinger -> cdp, eig_init, n: %d, m: %d, epsilon: %g, iter: %d, t: %g, relres: %g, relerr: %g\n', n, L*n, epsilon, iter, t, norm(abs(Azw).^2-Y,'fro')/norm(Y,'fro'), ... 
        norm(x - exp(-1i*angle(trace(x'*zw))) * zw, 'fro')/norm(x,'fro'));
    fclose(fid);
    
    %--- Gerchberg-Saxton ---%
    maxiter = 2500;
    tol = epsilon;
    zg = z;
    
    tic;
    for iter = 1 : maxiter
      Azg = A(zg);
      
      if norm(abs(Azg).^2-Y,'fro')/norm(Y,'fro') <= tol
        break;
      end
      
      Yzg = sqrt(Y).*(Azg./abs(Azg));
      zg = Ainv(Yzg);
    end
    t = toc;
    
    % output results of Gerchberg-Saxton
    fname = 'GerSax_cdp_eig.txt';
    fid = fopen(fname,'a');
    fprintf(fid, 'GerSax -> cdp, eig_init, n: %d, m: %d, epsilon: %g, iter: %d, t: %g, relres: %g, relerr: %g\n', n, L*n, epsilon,iter, t, norm(abs(Azg).^2-Y,'fro')/norm(Y,'fro'), ... 
        norm(x - exp(-1i*angle(trace(x'*zg))) * zg, 'fro')/norm(x,'fro'));
    fclose(fid);
   
    % --- Kaczmarz --- %
    maxiter = 500;
    tol = epsilon;
    zk = z;
    
    tic;
    for iter = 1 : maxiter
      maxrelres = 0;
      for r = 1 : n*L
        nrm2 = norm(Amat(r,:))^2;
        Arzk = Amat(r,:)*zk;
        
        maxrelres = max(maxrelres,abs(yv(r)-abs(Arzk)^2)/yv(r));
        
        zk = zk + (Arzk/abs(Arzk)*sqrt(yv(r))-Arzk)*Amat(r,:)'/nrm2;
      end
      
      if maxrelres < tol
        break;
      end
    end
    t = toc;
    
    % output results of Kaczmarz method
    fname = 'kaczmarz_cdp_eig.txt';
    fid = fopen(fname,'a');
    fprintf(fid, 'Kaczmarz -> cdp, rand_init, n: %d, m: %d, epsilon: %g, iter: %d, t: %g, relres: %g, relerr: %g\n', n, L*n, epsilon,iter, t, norm(yv-abs(Amat*zk).^2)/norm(yv), ... 
        norm(x - exp(-1i*angle(trace(x'*zk))) * zk, 'fro')/norm(x,'fro'));
    fclose(fid);
    
    
    %--- Block Kaczmarz ---%
    maxiter = 500;
    tol = epsilon;
    blockSS = round(n);
    
    for bls = 1 : length(blockSS)
      blockSize = blockSS(bls);
      nblocks = ceil(n/blockSize);
      zbk = z;  
      
      tic;
      for iter = 1 : maxiter
        maxrelres = 0;
        for kk = 1 : L
          d = conj(mask(:,kk));
          for r = 1 : nblocks
            Gamma = (r-1)*blockSize+1:min(r*blockSize,n);
            Azbk = Fd(zbk,d,Gamma);
          
            maxrelres = max(maxrelres,norm(Y(Gamma,kk)-abs(Azbk).^2,'fro')/norm(Y(Gamma,kk),'fro'));
          
            zbk = zbk + invFd(((Azbk./abs(Azbk)).*sqrt(Y(Gamma,kk))-Azbk),d,Gamma,n);
          end
        end
        
        if maxrelres < tol
           break;
        end
      end
      t = toc;
      
      % output results for block kaczmarz
      fname = ['block_kaczmarz_cdp_' num2str(blockSize) '_eig.txt'];
      fid = fopen(fname,'a');
      alg_name = ['Block Kaczmarz-' num2str(blockSize)];
      fprintf(fid, [alg_name ' -> cdp, eig_init, n: %d, m: %d, epsilon: %g, iter: %d, t: %g, relres: %g, relerr: %g\n'], n, L*n, epsilon,iter, t, norm(abs(A(zbk)).^2-Y,'fro')/norm(Y,'fro'), ... 
        norm(x - exp(-1i*angle(trace(x'*zbk))) * zbk, 'fro')/norm(x,'fro'));
      fclose(fid);
    end
  %----------------------
  end
end
