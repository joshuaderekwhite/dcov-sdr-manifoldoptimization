function [Syx, beta] = mmdcov(X,Y)
    % Precompute B, Z from X and Y step 2
    constants = struct();
    constants.X = X;
    constants.Y = Y;
    [constants.n,constants.p] = size(X);
    constants.d = size(constants.Y)(2);
    B=Kermat(constants.Y);
    constants.B=B-mean(B,2)-mean(B,1)+mean(mean(B));
    constants.Z=constants.X/sqrtm(cov(constants.X));   % Standardized X
    constants.Index=(Kermat(constants.Z)~=0);
    constants.epsilon=1e-10; % Perturbed constant

    [gamma, gamma_cost, info, options] = StiefelOpt(constants);
    gamma_est = real(gamma);

    % Sort largest to smalles and remove smallest values.
    [~,ind] = sort(-mean(abs(real(gamma)')));
    rank = sum(abs(eig(X'*X)) > 1E-7);
    Syx = ind(1:rank);
    beta = real(gamma_est(Syx,:)'*inv(sqrtm(cov(constants.X(:,Syx)))))';
end
