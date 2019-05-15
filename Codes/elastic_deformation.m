function reg_img = deformation(i1, i2, sigma, lambda, n_iter, step_size)

%i1 = 255.*i1;
%i2 = 255.*i2;

lookup_field_y = repmat([1:size(i1,2)] ,size(i1,1),1);
lookup_field_x = repmat([1:size(i1,1)]',1,size(i1,2));

def_field_x = zeros(size(i1));
def_field_y = zeros(size(i1));


[r,c] = size(i1);
for u = 0:r-1
    for v = 0:c-1
        L(u+1,v+1) = -2*(cos(2*pi*u/r)+cos(2*pi*v/c))+4;
    end
end

%L(1,1) = 0.1;

[i2_x, i2_y] = grad(i2);

for x = 1:n_iter
    %i2_u = apply_field(i2, lookup_field_x, def_field_x, lookup_field_y, def_field_y);
    
    hx = lookup_field_x + def_field_x;
    hy = lookup_field_y + def_field_y;
    
    i2_u = interp2(i2,hy,hx,'linear',0);
    diff = i2_u - i1;
    
    if (rem(x,10) == 0)
        figure(1);
        imagesc(i2_u);
        drawnow;

        figure(2);
        imagesc(diff);
        drawnow;
        ssd = sum(diff(:).^2)

        figure(3); 
        plot(hy(1:5:end,1:5:end),hx(1:5:end,1:5:end),'b');
        hold on;
        plot(hy(1:5:end,1:5:end)',hx(1:5:end,1:5:end)','b');
        hold off;
     
    end
    %[grad_i2_x, grad_i2_y] = grad(i2_u);
    grad_i2_x = interp2(i2_x,hy,hx,'linear',0);
    grad_i2_y = interp2(i2_y,hy,hx,'linear',0);
    
    der_y = diff .* grad_i2_x;
    der_x = diff .* grad_i2_y;
    
    % Now change to Fourier Domain.
   
    der_x_f = fftn(der_x)./(L+0.0001);
    der_y_f = fftn(der_y)./(L+0.0001);
    
    def_field_x = (def_field_x + (1/(sigma^2 * lambda))*step_size.*real(ifftn(der_x_f)));
    def_field_y = (def_field_y + (1/(sigma^2 * lambda))*step_size.*real(ifftn(der_y_f)));
    
   
end

%reg_img = apply_field(i2, lookup_field_x, def_field_x, lookup_field_y, def_field_y);
reg_img = i2_u;