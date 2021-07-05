% This script file conduct PR from Gaussian model with sparse outliers

clear
close all

Params.n2          = 1;
if isfield(Params, 'n1')          == 0,  Params.n1          = 200; end             % signal dimension
if isfield(Params, 'm')           == 0,  Params.m           = round(6* Params.n1);  end     % number of measurements
if isfield(Params, 'cplx_flag')   == 0,  Params.cplx_flag   = 1;    end             % real: cplx_flag = 0;  complex: cplx_flag = 1;


if isfield(Params, 'alpha_lb')    == 0,  Params.alpha_lb    = 0.3;  end
if isfield(Params, 'alpha_ub')    == 0,  Params.alpha_ub    = 5;    end
if isfield(Params, 'alpha_h')     == 0,  Params.alpha_h     = 5;    end
if isfield(Params, 'alpha_y')     == 0,  Params.alpha_y     = 3;    end
if isfield(Params, 'T')           == 0,  Params.T           = 900;  end  % number of iterations
if isfield(Params, 'mu')          == 0,  Params.mu          = 7;  end		% step size / learning parameter
if isfield(Params, 'npower_iter') == 0,  Params.npower_iter = 100;   end	% number of power iterations
Params.order = 4; Params.epsilon = 1; Params.trun_err=1e-5;
Params.no_ini = 0;


n           = Params.n1;
m           = Params.m;
cplx_flag	= Params.cplx_flag;  % real-valued: cplx_flag = 0;  complex-valued: cplx_flag = 1;
% display(Params)

%% Make signal and data (noiseless)
x = randn(n,1)  + cplx_flag * 1i * randn(n,1);
Amatrix = (randn(m,n) + cplx_flag * 1i * randn(m,n)) / (sqrt(2)^cplx_flag);
A  = @(I) Amatrix  * I;
At = @(Y) Amatrix' * Y;
b = abs(A(x));
Params.Amat = Amatrix;
ErrorNorm= 2*norm(x);


s = 0.1;

[~,out_idx]  =sort(randn(m,1));
outliers  = zeros(m,1); outliers(out_idx(1:ceil(m*s))) = 1;

NumOutliers=norm(outliers,1);

eta=ErrorNorm*rand(m,1);
outliers=outliers.*eta;
b_clean   = b;
b_outlier = b_clean+outliers;

snr= inf;
NoiseNorm= 10^(-snr/20) * norm(b);
noise=randn(m,1);
noise = noise / norm(noise) * NoiseNorm;
b_noise = b_outlier+noise;

y = b_noise .^2;   % intensity measurements

%% Check results and Report Success/Failure
tic
        [outs] = med_saf1d(y,  x, Params, A, At); Relerrs=outs.Relerrs;
%         [outs] = saf1d(y, x, Params, A, At);Relerrs=outs.Relerrs;
        [Relerrs] = medianTWF(y, x, Params, A, At);
        [Relerrs] = medianRWF(y, x, Params, A, At);
dt = toc;


if Relerrs(end)<=Params.trun_err
    state = 'SUCCESS';
else
    state = 'FAIL';
end

fprintf('Final relative err %3.2e after %d iters within %2.2f sec \n Params: outlier fraction s:%2.2f and eta:%3.2f, %s\n', Relerrs(end),length(Relerrs),s,dt, ErrorNorm/norm(x)^2,state)

semilogy(Relerrs)
drawnow
xlabel('Iteration'), ylabel('Relative error (log10)'), ...
    title('Relative error vs. iteration count')
hold on
