function [ ones_v ] = check_eignorm( eig,scale,TYPE )
% [ ones_v ] = check_eignorm( eig,scale,TYPE )
% Check that the eigenfunctions are normalized according to equation ##
% from Dahlen & Tromp - Theoretical Global Seismology

rho = eig(1).rho;
r = eig(1).r;

ones_v = ones(length(eig),1);
for iper = 1:length(eig)
    if TYPE == 'S'
        % Spheroidal Eigenfunctions
        omega = eig(iper).w;
        u = eig(iper).u;
        v = eig(iper).v;
        ones_v(iper) = trapz(r , rho.*((scale*u).^2 + (scale*v).^2).*r.^2) * omega^2;
    elseif TYPE == 'T'
        % Toroidal Eigenfunctions
        omega = eig(iper).w;
        wl = eig(iper).wl;
        ones_v(iper) = trapz(r , rho.*((scale*wl).^2).*r.^2) * omega^2;
    end
end


end

