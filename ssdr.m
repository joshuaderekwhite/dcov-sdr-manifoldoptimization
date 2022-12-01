function [X, beta, labels] = ssdr(X,Y,labels)
    centerX = mean(X);
    centerY = mean(Y);
    X = X-centerX;
    Y = Y-centerY;

    count = 0;
    reduced = true;
    Syx = [];
    while (reduced)
        [n,p] = size(X);
        Syxi = [];
        for section = 1:ceil(p/(n - ceil(p/n)))
            pnew = ((section - 1)*(n-1)+1):min(section*(n-1),p);
            X_sec = X(:,pnew);
            [Syx_sec, beta_sec] = mmdcov(X_sec,Y);
            Syx_sec = Syx_sec + (min(pnew) - 1);
            % if (size(Syx_sec)(2) == 0)
            %     %disp('No further reductions after %i iterations. Subspace dimension is (%i X %i)', )
            %     reduced = false;
            % end
            Syxi = union(Syxi, Syx_sec);
        end
        X = X(:,Syxi);
        centerX = centerX(Syxi);
        labels = labels(Syxi);
        Syx = union(Syx, Syxi);
        rank = sum(abs(eig(X'*X)) > 1E-7);
        if (rank(X) == size(X)(2))
            reduced = false;
        end
    endwhile
    X = X + centerX;
    beta = beta_sec;
end