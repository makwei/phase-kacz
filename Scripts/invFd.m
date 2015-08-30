% partail inverse Fourier transform
function x = invFd(f,d,Omega,n)
ff = zeros(n,1);
ff(Omega) = f;
x = (1./d).*ifft(ff);