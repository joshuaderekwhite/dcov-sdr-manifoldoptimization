function [gamma, gamma_cost, info, options]=StiefelOpt(constants)

    % Create the problem structure
    St = stiefelcomplexfactory(constants.p, constants.d);
    problem.M = St;

    function [Q,L] = qlMatrices(gamma, constants)
        % Calculate Q, L iteratively during linesearch: step 4-5
        % Coefficients of Quadratic Term in Surrogate Function
        size(gamma);
        xgamma=Kermat(constants.Z*gamma);
        xgamma(constants.Index)=1./(xgamma(constants.Index)+constants.epsilon);
        xgamma=xgamma.*constants.Index;
        C=(constants.B.*(constants.B<0)).*xgamma;
        Q=constants.Z'*(2*(diag(sum(C,2))-C)/(constants.n^2))*constants.Z;
        
        % Coefficients of Linear Term in Surrogate Function
        D=(constants.B.*(constants.B>0)).*xgamma;
        L=constants.Z'*(2*(diag(sum(D,2))-D)/(constants.n^2))*constants.Z*gamma;
    end

    % % Cost Function and Gradient
    % problem.cost = @(gamma) -real((0.5)*trace(gamma'*Q*gamma) + trace(gamma'*L));
    % problem.egrad = @(gamma) -(Q*gamma + L);
    % %checkgradient(problem);
    % problem.ehess = @(gamma,xi) -(Q*xi);
    % %checkhessian(problem);

    # Advanced method from: https://github.com/NicolasBoumal/manopt/blob/master/examples/elliptope_SDP.m
    function store = prepare(gamma, store)
        if ~isfield(store, 'Q')
            [store.Q, store.L] = qlMatrices(gamma, constants);
        end
    end

    function store = update(gamma, store)
        [store.Q,store.L] = qlMatrices(gamma, constants);
    end

    % Define the cost function to be /maximized/.
    problem.cost = @cost;
    function [f, store] = cost(gamma, store)
        store = prepare(gamma, store);
        f = -real((0.5)*trace(gamma'*store.Q*gamma) + trace(gamma'*store.L));
    end

    % Define the Riemannian gradient.
    problem.egrad = @grad;
    function [G, store] = grad(gamma, store)
        store = prepare(gamma, store);
        G = -(store.Q*gamma + store.L);
    end

    % Define the Hessian
    problem.ehess = @hess;
    function [H, store] = hess(gamma, xi, store)
        store = prepare(gamma, store);
        H = -(store.Q*xi);
    end

    % Solve
    options.solver = @conjugategradient;
    options.maxiter = 10;

    % linesearch using as reference:
    % https://www.manopt.org/reference/examples/low_rank_tensor_completion.html
    problem.linesearch = @linesearch_helper;
    function [gamma_t, store] = linesearch_helper(gamma, xi, store)
        # Calculate gamma from xi iteratively: steps 7-12
        s = 1; alpha = 1e-20;
        while PDistCov(constants.Z*qr(gamma + s*xi,0)(1),constants.Y,constants.epsilon) < ...
            PDistCov(constants.Z*gamma,constants.Y,constants.epsilon) + ...
            alpha*s*norm(xi,'fro')^2
            s = 0.5*s;
        end
        gamma_t = qr(gamma + s*xi,0)(1);
        store = update(gamma_t, store);
    end
    % options.statsfun = statsfunhelper('QLupdate', @(gamma, store) update(gamma, store));
    [gamma, gamma_cost, info, options] = manoptsolve(problem, [], options);

end

function dcovXY = PDistCov(x,y,ep)

    % Input
    % x: an n by p matrix
    % y: an n by q matrix
    
    % Output
    % dcovXY: the perturbed distance covariance between the x vector and the
    % y vector;
    
    
    %------------------------------------------------initialize the calculation
    sy=sum(y.^2,2);
    B=sqrt(real((bsxfun(@minus,sy',bsxfun(@minus,2*(y*y'),sy))))); % b_kl=||y_k-y_l||
    B=B-mean(B,2)-mean(B,1)+mean(mean(B)); % normalization
    
    sx=sum(x.^2,2);
    A=sqrt(real((bsxfun(@minus,sx',bsxfun(@minus,2*(x*x'),sx))))); % a_kl=||x_k-x_l||
    A=A-ep.*log(1+A./ep); % Perturbed Objective Function
    %------------------------------------------------
    
    dcovXY =mean(mean(A.*B));
    
    
end