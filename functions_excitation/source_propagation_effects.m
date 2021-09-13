% Multiply source excitation amplitude by factors that account for propagation 
% effects associated with attenuation and geometrical spreading
%
% jbrussell - 09/21
%
% A = raw source excitation amplitude
% w = angular frequency (rad/s)
% X = epicentral distance (km)
% phv = phase velocity (km/s)
% grv = group velocity (km/s)
% Q = quality factor

function A_corr = source_propagation_effects(A,w,X,phv,grv,Q)
    
    % wavenumber
    k = w./phv
    
    % Calculate attenuation factor
    fac_atten = exp(-w.*X./(2.*grv.*Q));
                    
    % Calculate geometrical spreading factor
    fac_gs = (8*pi*k.*abs(sind(km2deg(X))))^(-0.5) ./ (grv.*phv);
    
    A_corr = A .* fac_atten .* fac_gs;
end
