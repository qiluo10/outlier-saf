%% Implementation of the median truncated Wirtinger Flow (median-TWF) algorithm which is adapted from 
%  TWF by Y. Chen and E. J. Candès.
%% The author is Huishuai Zhang.

function [Relerrs] = medianRWF(y, x, Params, A, At)    
%% Initialization
    Relerrs=zeros(Params.T+1,1);
    npower_iter = Params.npower_iter;           % Number of power iterations 
    z0 = randn(Params.n1,Params.n2); z0 = z0/norm(z0,'fro');    % Initial guess 
    normest = sqrt(median(y(:))/0.455);    % Estimate norm to scale eigenvector  
    ymag    = sqrt(y);
    mu = 0.8;
    alpha_h = 3;
    
    for tt = 1:npower_iter                     % Truncated power iterations
        ytr = y.* (abs(y) <= Params.alpha_y^2 * normest^2 );
        z0 = At( ytr.* A(z0) ); z0 = z0/norm(z0,'fro');
    end
    
    z = normest * z0;                   % Apply scaling 
    Relerrs(1) = norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro'); % Initial rel. error
    
    %% Loop

    m=Params.m;
    for t = 1: Params.T
        Az=A(z);
        absyz=abs(Az);

       
        
        
        delta_y = absyz - ymag;
        Mt = median(abs(delta_y))*alpha_h;
        Eh = abs(delta_y) <= Mt;
        grad  = 1/m* At(  ( Az-ymag.*sign(Az) ).*Eh  );
        z = z - mu * grad;             % Gradient update
        
        currenterrs=norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro');
        
        if currenterrs< Params.trun_err
            break;
        end
        Relerrs(t+1) = currenterrs;  
    end
    for t=t:Params.T,
        Relerrs(t+1)= currenterrs;  
    end
