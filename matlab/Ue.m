function wave = Ue(z_,w_,sigma_,a_)
% PARABOLIC BEAM (even)
wave = 1 / sqrt(2) / pi * abs2(gammac(0.25 + 0.5j*a_))*Pe(sigma_*z_,a_) .* Pe(sigma_*w_,-a_);

end