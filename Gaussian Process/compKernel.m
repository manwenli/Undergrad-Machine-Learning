function [kernel_val] = compKernel(x,xp,tao_sq)
    temp = -(x-xp).^2/(2*tao_sq);
    kernel_val = exp(temp);

end

