% process all the data

% process Gaussian transition data
folder = '../Gaussian/';
files = dir(folder);
for ll = 3 : length(files)
  process_gaussian_data([folder files(ll).name]);
end



% process CDP transition data 1D
folder = '../CDP1d/';
files = dir(folder);
for ll = 3 : length(files)
  process_cdp_data_1d([folder files(ll).name]);
end



% process CDP transition data 2D
folder = '../CDP2d/';
files = dir(folder);
for ll = 3 : length(files)
  process_cdp_data_2d([folder files(ll).name]);
end



% process Unitary transition data
folder = '../Unitary/';
files = dir(folder);
for ll = 3 : length(files)
  process_gaussian_data([folder files(ll).name]);
end



% process timing data
folder = '../Timing/';
files = dir([folder '*_guassian_*']);
for ll = 1 : length(files)
  process_gaussian_data([folder files(ll).name]);
end

files = dir([folder '*_cdp_*']);
for ll = 1 : length(files)
  process_cdp_data_1d([folder files(ll).name]);
end

files = dir([folder '*_cdp2d_*']);
for ll = 1 : length(files)
  process_cdp_data_2d([folder files(ll).name]);
end



% process noise data
folder = '../Noise/';
files = dir([folder '*_guassian_*']);
for ll = 1 : length(files)
  process_gaussian_noise_data([folder files(ll).name]);
end

files = dir([folder '*_cdp_*']);
for ll = 1 : length(files)
  process_cdp_noise_data_1d([folder files(ll).name]);
end

files = dir([folder '*_cdp2d_*']);
for ll = 1 : length(files)
  process_cdp_noise_data_2d([folder files(ll).name]);
end



% process image data
folder = '../Image/';
files = dir([folder '*_gal*']);
for ll = 1 : length(files)
  process_cdp_gal_data([folder files(ll).name]);
end

files = dir([folder '*_mol*']);
for ll = 1 : length(files)
  process_cdp_data_2d([folder files(ll).name]);
end
