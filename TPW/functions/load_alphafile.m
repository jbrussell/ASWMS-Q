function [mat] = load_alphafile(file,alpha_ref)
%Load alpha file and interpolate to desired grid
%
data = load(file);

alpha = data(:,1);
alpha_std = data(:,2);

% Add reference alpha back in
alpha = alpha + alpha_ref;

mat.alpha = alpha;
mat.alpha_std = alpha_std;

end

