%% Implementation of the Smooth Amplitude Flow algorithm proposed in the paper
%  Phase Retrieval via Smooth Amplitude Flow
%  by  Q. Luo, H. Wang, and S. Lin from NUDT, China
%  The code below is adapted from implementation of the Wirtinger Flow
%  algorithm implemented by E. Candes, X. Li, and M. Soltanolkotabi.

function [outs, z] = med_saf1d(y, x, Params, A, At)



%% Initialization
if Params.cplx_flag
    z0 = randn(Params.n1, Params.n2) + 1i*randn(Params.n1, Params.n2);
else
    z0 = randn(Params.n1, Params.n2);
end

z0      = z0 / norm(z0, 'fro');    % Initial guess
% AmatT  = Amatrix';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normest = sqrt(median(y(:))/0.455);    % Estimate norm to scale eigenvector      % Estimate norm to scale eigenvector
m       = Params.m;
ymag    = sqrt(y);
order   = Params.order;

epsilon = Params.epsilon;
ymag_order = epsilon^order * ymag.^order;
tt = (1+epsilon^order)^(1/order)*ymag;
trun_err= Params.trun_err;

if isfield(Params,'s')
    p_prior = ceil(  0.9*100*(1-Params.s));
else
    p_prior = 50;
end


%% The weighted maximal correlation initialization can be computed using power iterations
%% or the Lanczos algorithm, and the latter performs well when m/n is small

if ~Params.no_ini
    for i = 1:Params.npower_iter                     % Truncated power iterations
        ytr = y.* (abs(y) <= Params.alpha_y^2 * normest^2 );
        z0 = At( ytr.* A(z0) ); z0 = z0/norm(z0,'fro');
    end
end

z = normest * z0;


Relerrs = norm(x - exp(-1i * angle(trace(x' * z))) * z, 'fro') / norm(x, 'fro'); % Initial rel. error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









beta = 0.2;
% alpha_h = Params.alpha_h;
alpha_h = 3;

for t = 1: Params.T
    tau = Params.mu;
    Az    = A(z);
    delta_y = abs(Az) - ymag;
    Mt = median(abs(delta_y))*alpha_h;
%     Mt = prctile(abs(delta_y), p_prior )*alpha_h;
%     Mt = median(abs(delta_y) )*alpha_h;
    Eh = abs(delta_y) <= Mt;
    
%     [f,gradf] = SAF_Objective(Eh); 
    
    
%     f_now = f(Az);
    
    f_now = 0.5 * norm( ((abs(Az).^order+ymag_order).^(1/order)-tt).* Eh   ).^2/m;
    
%     grad_now = At(gradf(Az));
    grad_now = At(((abs(Az).^order+ymag_order).^(1/order)-tt).*(abs(Az).^order+...
            ymag_order).^(1/order-1).* abs(Az).^(order-1).*sign(Az).* Eh/m);
    
    count = 0;
    ngrad2 = norm(grad_now)^2;
    
    while count<=1
        znew  = z - tau  * grad_now;
        Aznew = A(znew);
        
%         f_znew = f(Aznew);
        f_znew = 0.5 * norm( ((abs(Aznew).^order+ymag_order).^(1/order)-tt).* Eh   ).^2/m;
        
        if f_znew-f_now <= - tau*beta*ngrad2
            
            %             disp(count)
            break;
        end
        
        
        
        tau = tau*0.2;
        count = count+1;
    end
    z = znew;
    %     disp(tau)
    %     if count==0
    %         tau = tau*1.4;
    %     end
    
    %     if tau<2e-15
    %         break
    %     end
    
    errcnt = norm(x - exp(-1i * angle(trace(x' * z))) * z, 'fro') / norm(x, 'fro');
    
    Relerrs = [Relerrs; errcnt]; %#ok<AGROW>
    taus(t)    = tau;
    
    
    if errcnt<= trun_err
        break
    end
end


outs.Relerrs = Relerrs;
outs.taus    = taus; % output the stepth at each iteration



    function [f, gradf] = SAF_Objective(Eh)

        
        f = @(z) 0.5 * norm( ((abs(z).^order+ymag_order).^(1/order)-tt).* Eh   ).^2/m;
        gradf = @(z) ((abs(z).^order+ymag_order).^(1/order)-tt).*(abs(z).^order+...
            ymag_order).^(1/order-1).* abs(z).^(order-1).*sign(z).* Eh/m;
    end






end




