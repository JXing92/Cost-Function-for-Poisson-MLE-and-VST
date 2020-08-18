function [f,grad] = L1_VST(x,PSF_a,z,lambda)
% L1_VST calculates data fidelity term using Anscombe's VST
% More details in A proximal Iteration for deconvolving poisson noisy
% images using sparse representations

% z acquired image, after Anscombe's VST

% PSF_a: template
% x: solution
% lambda: L1 regularization term, set to zero for non-regularized
% processing
sim_im=PSF_a*x+3/8;
f=sum((z-2*sqrt(sim_im)).^2)./2500+lambda*norm(x(1:end-1),1);
grad=sum(-2.*z./sqrt(sim_im).*PSF_a+4.*PSF_a,1)./2500+lambda.*[sign(x(1:end-1));0].';
end
