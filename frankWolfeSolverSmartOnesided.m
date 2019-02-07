function [X_proj_new,Efinal] = frankWolfeSolverSmartOnesided(dist1,dist2,params)
params.null = [];
numIters = getoptions(params,'numIters',100);
Xinit = getoptions(params,'Xinit',eye(size(dist1,1)));
k = size(Xinit,1);
n = size(Xinit,2);
row_max = getoptions(params,'row_max',[]);
X_trans = getoptions(params,'X_trans',[]);


% calc translation vector
if isempty(row_max)
    X_trans = -2*sum(dist1,2)*sum(dist2,1)+...
    sum(dist1.^2,2)*n*ones(1,n)+...
    ones(k,1)*k*sum(dist2.^2,1);
    row_max = max(X_trans,[],2);
    X_trans = repmat(row_max,[1,n]);
end

cX = Xinit;
iter = 1;
E = [];
contFlag = 1;

while (iter <=  numIters) && contFlag
    E(iter) = trace(dist1.^2*cX*ones(n)*cX')+trace(ones(k)*cX*dist2.^2*cX')-2*trace(dist1*cX*dist2*cX');
    
    orig_grad = -4*dist1*cX*dist2...
        +2*dist1.^2*cX*ones(n)...
        +2*ones(k)*cX*dist2.^2;
    [~,cXind] = max(cX,[],2);
    differences = orig_grad(sub2ind([k n], (1:k)', cXind)) - min(orig_grad,[],2);
    [translations,minidx] = sort(differences);  
    translations = [0;sort(translations)./(2*row_max(minidx))-10^-8];
    
    for jj = 1:numel(translations)
        % calc gradient
        % d/dX (tr(A.^2*X*11^T*X^T)+tr(11^TXB.^2X^T)-2tr(AXBX^T))=
        % -4*A*X*B+2*A.^2*X*1*1'+2*1*1'*X*B.^2
        grad = orig_grad-2*translations(jj)*X_trans.*cX;
        % solve linear assignmnet (1-sided MRF version)
        [~, inds] = min(grad,[],2);
        inds = sub2ind([k n], (1:k)', inds); %inds of ones in solution
        nX = zeros(size(cX));
        nX(inds)=1;
        nE = trace((dist1.^2*nX*ones(n)+ones(k)*nX*dist2.^2-2*(dist1*nX*dist2))*nX');
        doProgress = nE < E(iter);
        if doProgress
            cX = nX;
            break
        end
        if  jj == numel(translations)
            contFlag = 0;            
        end
    end
    iter = iter+1;
end


X_proj_new = cX;
Efinal = E(end);
end
