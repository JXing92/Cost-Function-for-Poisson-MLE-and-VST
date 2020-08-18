function hessian = L1_VST_Hessian(x,lambda,C,d)
%L1_VST_Hessian calculates Hessian for L1_VST.m
sim_im=((C*x+3/8).^(-3/4));
Temp=bsxfun(@times,C,d.*sim_im);
hessian=gather(Temp.'*Temp)./2500;

end

