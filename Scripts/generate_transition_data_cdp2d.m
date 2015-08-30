clear all
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));

% problem size
n1 = 256;
n2 = 256;

D = [2 3 4 5 6 7 8 9 10 11 12];  

num_tests = 50;

dict = [1i -1i 1 -1];

for ll = 1 : length(D)
  L = D(ll);
  for num = 1 : num_tests
    % generate the mask
    mask = zeros(n1,n2,L);
    for p = 1 : L
      for i = 1 : n1
        for j = 1 : n2
          mask(i,j,p) = dict(randi(4));
        end
      end
    end
    
    temp = rand(size(mask));
    mask = mask .* ( (temp <= 0.2)*sqrt(3) + (temp > 0.2)/sqrt(2) );
    
    % define the operators
    A = @(I)  fft2(conj(mask) .* reshape(repmat(I,[1 L]), size(I,1), size(I,2), L));  
    At = @(Y) mean(mask .* ifft2(Y), 3);                      
    denom = sum(abs(mask).^2,3); 
    Ainv = @(Y) sum(mask.*ifft2(Y),3)./denom;
    
    % generate data
    x = randn(n1,n2) + 1i*randn(n1,n2);
    Y = abs(A(x)).^2; 
    
    % initialization
    npower_iter = 50;                           
    z0 = randn(n1,n2); z0 = z0/norm(z0,'fro'); 
    for tt = 1:npower_iter, 
      z0 = At(Y.*A(z0)); z0 = z0/norm(z0,'fro');
    end
    normest = sqrt(sum(Y(:))/numel(Y));
    z = normest * z0;
    
    % random initialization
    %z = randn(n1,n2)+1i*randn(n1,n2);
   
    % --- Wirtinger flow --- %
    maxiter = 2500;          
    tol = 1e-8;
    tau0 = 330;                         
    mu = @(t) min(1-exp(-t/tau0), 0.4);    
    zw = z;
    
    tic;
    for iter = 1 : maxiter
      Azw = A(zw);
      
      diff = abs(Azw).^2-Y;
      if norm(diff(:),'fro')/norm(Y(:),'fro') <= tol
        break;
      end
      
      C  = (abs(Azw).^2-Y) .* Azw;
      grad = At(C);
      zw = zw - mu(iter)/normest^2 * grad;
    end
    t = toc;
    
    % output results of wirtinger flow
    fname = 'wirtinger_cdp2d_eig.txt';
    fid = fopen(fname,'a');
    diff = abs(Azw).^2-Y;
    fprintf(fid, 'Wirtinger -> cdp2d, eig_init, n1: %d, n2: %d, m: %d, iter: %d, t: %g, relres: %g, relerr: %g\n', n1, n2, L*n1*n2, iter, t, norm(diff(:),'fro')/norm(Y(:),'fro'), ... 
        norm(x - exp(-1i*angle(trace(x'*zw))) * zw, 'fro')/norm(x,'fro'));
    fclose(fid);
    
    
    %--- Gerchberg-Saxton ---%
    maxiter = 2500;
    tol = 1e-8;
    zg = z;
    
    tic;
    for iter = 1 : maxiter
      Azg = A(zg);
      
      diff = abs(Azg).^2-Y;
      if norm(diff(:),'fro')/norm(Y(:),'fro') <= tol
        break;
      end
      
      Yzg = sqrt(Y).*(Azg./abs(Azg));
      zg = Ainv(Yzg);
    end
    t = toc;
    
    % output results of Gerchberg-Saxton
    fname = 'GerSax_cdp2d_eig.txt';
    fid = fopen(fname,'a');
    diff = abs(Azg).^2-Y;
    fprintf(fid, 'GerSax -> cdp2d, eig, n1: %d, n2: %d, m: %d, iter: %d, t: %g, relres: %g, relerr: %g\n', n1, n2, L*n1*n2, iter, t, norm(diff(:),'fro')/norm(Y(:),'fro'), ... 
        norm(x - exp(-1i*angle(trace(x'*zg))) * zg, 'fro')/norm(x,'fro'));
    fclose(fid);
   
    %--- Block Kaczmarz ---%
    maxiter = 500;
    tol = 1e-7;
    zbk = z;
    
    tic;
    for iter = 1 : maxiter
      maxrelres = 0;
      for kk = 1 : L
        d = conj(mask(:,:,kk));
        Azbk = fft2(d.*zbk);
        
        maxrelres = max(maxrelres,norm(Y(:,:,kk)-abs(Azbk).^2,'fro')/norm(Y(:,:,kk),'fro'));
        
        zbk = (1./d).* ifft2(((Azbk./abs(Azbk)).*sqrt(Y(:,:,kk))));
      end
      
      if maxrelres < tol
        break;
      end
    end    
    t = toc;
      
    % output results for block kaczmarz
    fname = 'block_kaczmarz_cdp2d_eig.txt';
    fid = fopen(fname,'a');
    diff = abs(A(zbk)).^2-Y;
    fprintf(fid, 'Block Kaczmarz -> cdp2d, eig_init, n1: %d, n2: %d, m: %d, iter: %d, t: %g, relres: %g, relerr: %g\n', n1, n2, L*n1*n2, iter, t, norm(diff(:),'fro')/norm(Y(:),'fro'), ... 
      norm(x - exp(-1i*angle(trace(x'*zbk))) * zbk, 'fro')/norm(x,'fro'));
    fclose(fid);
  %----------------------
  end
end