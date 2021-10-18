function [mat] = load_alphafile(file)
%Load alpha file and interpolate to desired grid
%
data = load(file);

alpha = data(:,1);
alpha_std = data(:,2);

mat.alpha = alpha;
mat.alpha_std = alpha_std;

end

