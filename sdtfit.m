function [w,nll,grad,H] = sdtfit(A,b,r1,r2,plotflag)
% [w,nll,grad,H] = sdtfit(A,b,r1,r2)
%
%   Maximum likelihood fit of a sdt-like categorization model.
%   
%   A is the designmatrix and b the corresponding vector, r1 is the number
%   of category 1 responses, r2 is the number of category 2 responses. nll
%   is the negative log likelihood, grad is the gradient at the solution
%   and H is the Hessian.
%
%   Example: Say we have 4 stim and the first two (with means mu1, mu2) are 
%   category 1 and the second two (mu3, mu4) are category two. The subject
%   has one criterion to decide whether it's category 1 or category 2 where
%   everything on the left of the criterion is a category 1 response. Hence 
%   we would like to find a vector w=[mu1 mu2 mu3 mu4 c]'. r1 is a vector
%   with category 1 responses for each stimulus r1 = [r11 r12 r13 r14]' and
%   the same for r2:
%
%     w  = [mu1 mu2 mu3 mu4 c]'
%     r1 = [r11 r12 r13 r14]'
%     r2 = [r21 r22 r23 r24]'
%
%   The design matrix A for this case is:
%
%     A  = [-1  0  0  0  1; ...
%            0 -1  0  0  1; ...
%            0  0 -1  0  1; ...
%            0  0  0 -1  1]
%
%   and b contains only zeros,
%   so that we have the following model to predict the response
%   probabilities p1 for r1 and p2 for r2:
%
%     p1 = normcdf(A*w+b)
%     p2 = 1 - p1;
%
%   Note that A does not have full rank. The origin is arbitrary, hence you
%   might either set one stimulus (say mu1) to be zero or choose c to
%   be zero by dropping the respective column from the designmatrix.
%
%   The optimization problem is convex but the unique maximum might lie
%   at infinity. By adding a small regularizer as a penality to the
%   likelihood function this problem is avoided. So strictly speaking, we
%   are doing penalized likelihood maximization here.
%
%   I've designed the code for the case where we have sparse (mostly
%   binary) design matrices A---as it happens in Thurstonian-scaling-like
%   situations as described in the example; and especially when you go to
%   trial-by-trial modeling of the data.
%
%   Version 1.1 --- March 2021, Frank J??kel, fjaekel@uos.de
%                               Christina Koß, christina.koss@tu-darmstadt.de
%

if nargin == 5
    RESULTPLOT = plotflag;
else
    RESULTPLOT = 0;
end

% Regularization parameter (this could also be a vector if you wanted
% different values for each dimension):
reg = 10^-5;

% Note that the designmatrix has to have "full" rank for the Newton
% method to work, hence you have to explicitly deal with the
% arbitrary origin in SDT/Thurstonian models by fixing one of the
% parameters to 0 (i.e. dropping this column in the design matrix).
if not(rank(A)==size(A,2))
    error('design matrix does not have independent columns')
end

% if you use a non-zero lapse the optimization problem becomes non-convex
lapse = 0;
[w,nll] = mlfit_newton(A,b,r1,r2,lapse,reg);

if RESULTPLOT
    figure
    % plot data
    z = A*w+b;
    p = r1./(r1+r2);
    plot(z,p,'b.');
    % plot fit
    hold on
    z = linspace(-3.5,3.5,1000);
    p = F(z,lapse);
    % plot confidence regions around fit assuming that all data points have
    % roughly the same number of trials N
    N = round(mean(r1+r2));
    lower = binoinv(0.025,N,p)/N;
    upper = binoinv(1-0.025,N,p)/N;
    plot(z,p,'b-');
    plot(z,lower,'b:');
    plot(z,upper,'b:');
    xlabel('predicted z')
    ylabel('category 1 response prob.')
    title(sprintf('nll = %6.8g',nll))
    xlim([min(z), max(z)])
end

if nargout > 2
    [~, grad, H] = negloglikelihood(w,A,b,r1,r2,lapse,reg);
end

% -------------------------------------------------------------------------

function p = F(z,lapse)
% The psychometric function
p = (1-lapse) * normcdf(z) + lapse/2;

function z = Finv(p,lapse)
% The inverse of the psychometric function
z = norminv(p-lapse/2)/(1-lapse);

function p = f(z,lapse)
% The first derivative of the psychometric function
p = (1-lapse)*normpdf(z);

function p = df(z,lapse)
% The second derivative of the psychometric function
% (in an earlier version I forgot the (-1), so be careful when you go back
% to an earlier version of this code!)
p = (1-lapse)*(-1)*z.*normpdf(z);

% -------------------------------------------------------------------------

function [w, sse] = lsqfit(A,b,r1,r2,lapse)
% least squares fit
p = r1./(r1+r2);
% scale the probabilites so that they have twice the lapse rate and we are
% sure stay in a nice regime where we don't get infinities anywhere
if lapse < 0.01;
    lapse = 0.01;
end
p = p*(1-2*lapse)+lapse;  
z = Finv(p,lapse);
if not(rank(A)==size(A,2))
    w = pinv(A)*(z-b);
else
    w = A\(z-b);
end   
% return the summed square error
sse = sum((z-F(A*w+b,lapse)).^2);

% -------------------------------------------------------------------------

function [nll, grad, H] = negloglikelihood(w,A,b,r1,r2,lapse,reg)
% regularized negative log of the likelihood function, also returns the
% gradient and the Hessian for parameters w if requested
if length(reg)==1
    reg = reg * ones(size(w));
end
z = A*w+b;
p = F(z,lapse);
q = 1-p;
nll = - sum(r1 .* log(p) + r2 .* log(q)) + 0.5 * (reg.*w)'*w;
if nargout > 1
    pdf = f(z,lapse);
    grad = - A' * (r1.*pdf./p - r2.*pdf./q) + reg.*w;
end
if nargout > 2
    dpdf = df(z,lapse);
    F1 = (dpdf.*p-pdf.^2) ./ p.^2;
    F2 = (-dpdf.*q-pdf.^2) ./ q.^2;
    B = sparse(1:size(A,1),1:size(A,1),r1.*F1 + r2.*F2);  
    H = full((-1) * A'*B*A) + diag(reg);
end

% -------------------------------------------------------------------------
function [w, nll] = mlfit_newton(A,b,r1,r2,lapse,reg,quick)
% maximum likelihood fit using Newton's method, but do a quick and dirty
% prefit for initializing the weights; if we have a lapse>0 parameter, it's 
% better to initialize with the ml-fit with a lapse of 0  
if lapse == 0
    w = lsqfit(A,b,r1,r2,lapse);
else
    w = mlfit_newton(A,b,r1,r2,0,reg,1);
end
% using sparse matrix multiplication speeds up the algorithm if A is
% sparse; in our applications where A is the design matrix with only a few
% 1s in each row this speed-up is substantial
A = sparse(A);
% calculate gradient and Hessian
[nll, grad, H] = negloglikelihood(w,A,b,r1,r2,lapse,reg);
% set a few optimization parameters
if nargin<7
    quick = 0;
end
if quick
    TolFun = 10^-1;
    TolDir = 10^-1;
else
    TolFun = eps;
    TolDir = 10^-5;
end
oldnll = Inf;
dir = Inf;
c = 10^-4;
i = 0;
PLOT = 0;
ILLCOUNT = 0; 
ILLLIMIT = 100;
while (oldnll-nll) > TolFun || norm(dir) > TolDir || isnan(nll)
    i = i+1;
    % plotting is helpful for debugging
    if PLOT 
        allnll(i) = nll;
        plot(allnll)
        drawnow
    end
    % calculate direction to change w. If H is ill-conditioned do gradient
    % descent instead of a Newton step
    if rcond(H) > eps
        dir = (-1) * H\grad;
        ILLCOUNT = 0;
    else
        dir = -grad;
        if ILLCOUNT==0
          fprintf('nll=%g, Hessian was ill-conditioned in step %d\n',nll,i)
        end
        ILLCOUNT = ILLCOUNT+1; 
    end
    % if we had too many ill-conditioned matrices it's probably hopeless!
    % so we give up and just return NaN
    if ILLCOUNT > ILLLIMIT
      fprintf('nll=%g, Hessian was ill-conditioned for %d steps -> stop optimization\n',nll,ILLLIMIT)
      break
    end
    % check that we are going downhill
    d = grad'*dir;
    if d>0
        % If the problem is convex, i.e. for lapse=0, this should not
        % happen. But it can happen if lapse>0. Just resort to slow
        % gradient descent then
        dir = -grad;
        d = grad'*dir;
        fprintf('Hessian is negative in step %d\n',i)
    end
    % choose stepsize by backtracking search so that it satisfies the first
    % Wolfe condition (e.g., see page 41f in Nocedal and Wright, 1999)
    stepsize = 1;
    while(negloglikelihood(w+stepsize*dir,A,b,r1,r2,lapse,reg) > ...
                nll+c*stepsize*d)
            stepsize = stepsize * 0.5;
    end
    % update w
    w = w + stepsize * dir;
    oldnll = nll;
    [nll, grad, H] = negloglikelihood(w,A,b,r1,r2,lapse,reg);
end