function [f,grad] = MLEL1( sol,PSF_a,im3D,lambda )
% Full-on MLE with L1 regularization, set lambda=0 for non-regularized
% processing

sim_im=PSF_a*sol;
f=gather(sum(sim_im-im3D.*log(sim_im))+lambda*norm(sol(1:end-1),1));
grad=gather(sum(bsxfun(@times,PSF_a,1-im3D./sim_im),1)+lambda.*[sign(sol(1:end-1));0].');


end
