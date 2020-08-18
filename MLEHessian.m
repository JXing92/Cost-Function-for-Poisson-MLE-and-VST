function hessian = MLEHessian(x,lambda,PSF_a,im3D)
% MLEHessian calculates the Hessian matrix of the maximum liklihood problem
sim_im=PSF_a*x;
Temp=bsxfun(@times,PSF_a,im3D./sim_im);
hessian=gather(Temp.'*Temp);

end