function [X_proj_new,Efinal] = frankWolfeSolver(dist1,dist2,params)
% tic
params.null = [];
numIters = getoptions(params,'numIters',100);
threshold = getoptions(params,'threshold',10^-6);
Xinit = getoptions(params,'Xinit',eye(size(dist1,1)));
doLineSearch = getoptions(params,'doLineSearch',0);
docheckTOptimality = getoptions(params,'docheckTOptimality',0);
scalingFactor = getoptions(params,'scalingFactor',10^3);


cX = Xinit;
diff = inf;
iter = 1;
E = [];
while (iter <=  numIters) && (diff > threshold)
    E(iter) = -2*trace(dist1*cX*dist2'*cX');
    
    % calc gradient d/dX (-2tr(AXBX^T))=-4AXB
    C = 4*dist1*cX*dist2;
    % solve linear assignmnet
    if isempty(scalingFactor)
        max_payoff_value = 10^7;
        scalingFactor = max_payoff_value/max(C(:));
    end
    Cscaled = sparse(C*scalingFactor+eps);
    [~, nX] = ...
        sparseAssignmentProblemAuctionAlgorithm(Cscaled);
    % do line search if necessary
    if doLineSearch && iter>=2
        % solve min 0<=t<=1 [(t*Cx+(1-t)*nextX)'*W(t*Cx+(1-t)*nextX))]
        f = @(t) (-2*trace(dist1*(t*cX+(1-t)*nX)*dist2*(t*cX+(1-t)*nX)'));
        t_opt = (-4*trace(dist1*nX*dist2*nX')+4*trace(dist1*nX*dist2*cX'))/...
            (-4*trace(dist1*cX*dist2*cX')-4*trace(dist1*nX*dist2*nX')+8*trace(dist1*nX*dist2*cX'));
        f_0 = f(0);
        f_1 = f(1);        
        a = (-4*trace(dist1*cX*dist2*cX')-4*trace(dist1*nX*dist2*nX')+8*trace(dist1*nX*dist2*cX'));       
        %  check that t_opt is indeed optimal
        if docheckTOptimality
            ts = linspace(-3,3,300);
            fs = [];
            for jj = 1:numel(ts)
                fs(jj) = f(ts(jj));
            end
            figure
            plot(ts,fs)
            hold on
            scatter(t_opt,f(t_opt))
        end
        
        if (0<=t_opt)&&(t_opt<=1)
            f_opt = f(t_opt);
        else
            f_opt = inf;
        end
        % select minimal value
        t = 0;
        if f_0<=f_1 && f_0<=f_opt
            t= 0;
        elseif f_1<=f_0 && f_1<=f_opt
            t=1;
        elseif f_opt<=f_0 && f_opt<=f_1
            t=t_opt;
        end
        if a>0
            nonConcaveIterations = nonConcaveIterations+1;
        end
        
        nX = t*cX+(1-t)*nX;
    end
    
    % prepare for next iteration
    diff = max(abs(nX(:)-cX(:)));
    %     fprintf('norm diff from last iter = %.4f...\n',norm(nX(:)-cX(:)))
    cX = nX;
    iter = iter+1;
end
X_proj_new = cX;
Efinal = E(end);
% toc
end
