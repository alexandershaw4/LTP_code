function [x,i] = solvefixedpoint(P,M,input)


% inital states
%------------------------------------------------------------------
ns   = size(P.A{1},1);                           % number of sources

if isfield(M,'x') && ~isempty(M.x)
    np = size(M.x,2);
else
    np   = size(P.H,1);                              % number of populations
end

%np = 7;

% find fp in presence of an input (DC), or not
if nargin < 3
    M.u = sparse(ns,1);%3
else
    M.u = input(1);
end
    
% create (initialise voltage at -50mV)
%--------------------------------------------------------------------------
if isfield(M,'x')
    nk = size(M.x,3);
else
    nk = 7;
end

x        = zeros(ns,np,nk) ;
x(:,:,1) = -70;
        
M.g   = {};
M.x   = x;
M.pE  = P;
M.n   = length(spm_vec(x));

% solve for steady state
%--------------------------------------------------------------------------
warning off
[x,i] = solvefixed(P,M);
f     = M.f;
warning on

end

function [x,i] = solvefixed(P,M)

% solve for fixed point
%------------------------------------------------------------------
ns    = size(P.A{1},1);     % number of sources (endogenous inputs)
a     = 2;                  % regulariser
dnx   = 0;
for i = 1:128
    
    % solve under locally linear assumptions
    %--------------------------------------------------------------
    [f,dfdx] = feval(M.f,M.x,M.u,P,M);
    dx       = - dfdx\f;
    
    % regularise
    %--------------------------------------------------------------
    ndx   = norm(dx,Inf);
    if ndx < dnx
        a = a/2;
    end
    dnx    = ndx;
    
    % update and convergence
    %--------------------------------------------------------------
    M.x    = spm_unvec(spm_vec(M.x) + exp(-a)*dx,M.x);
    if dnx < 1e-12, break, end
    
end

x = M.x;
end