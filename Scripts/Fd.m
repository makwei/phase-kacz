% partial Fourier transform
function f = Fd(x,d,Omega)
f = fft(d.*x);
if nargin == 3
  f = f(Omega);
end