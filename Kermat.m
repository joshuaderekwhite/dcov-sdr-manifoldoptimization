function K=Kermat(x)
    
    % Input
    % x: an n by p matrix
    % 
    % Output
    % K: an n by n matrix with B_{ij}=||x_i-x_j||_2
    
    
    sx=sum(x.^2,2);
    K=real((bsxfun(@minus,sx',bsxfun(@minus,2*(x*x'),sx))).^(1/2)); % K_ij=||x_i-x_j||_2
    
end