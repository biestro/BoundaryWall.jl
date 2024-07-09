function wave = Uo(z_,w_,sigma_,a_)
% PARABOLIC BEAM (even)
wave = 2 / (sqrt(2)* pi) * abs2(gammac(0.75 + 0.5j*a_))*Po(sigma_*z_,a_) .* Po(sigma_*w_,-a_);
end