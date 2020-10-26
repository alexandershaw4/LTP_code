function [f,J,Q,D,DV] = tc_dev_dense(x,u,P,M,m)
% State equations for an extended canonical thalamo-cortical neural-mass model.
%
% This model implements a conductance-based canonical thalamo-cortical circuit,
% with cytoarchitecture inspired by Gilbert & Wiesel (1983), Douglas & 
% Martin (2004) and Traub (2004) models.
%
% The equations of motion are Moris Lecar-esque equations, similar to Moran
% (2011), but with conductances for AMPA, NMDA, GABA-A, & GABA-B channels. 
% These 'channels' feature their own reversal poentials and rate constants:
%
% K  = -70           (Leak)
% Na =  60  / 4 ms   (AMPA)
% Cl = -90  / 16 ms  (GABA-A)
% Ca =  10  / 100 ms (NMDA)   + voltage mag switch
% B  = -100 / 200 ms (GABA-B)
% f  = -40
%
% FORMAT [f,J,Q,D] = atcm.tcm_nomh(x,u,P,M)
%
% x - states and covariances
%
% x(i,j,k)        - k-th state of j-th population of i-th source
%                   i.e., running over sources, pop. and states
%
%   population: 1  - Spint stellates (L4)
%               2  - Superficial pyramids (L2/3)
%               3  - Inhibitory interneurons (L2/3)     
%               4  - Deep pyramidal cells (L5)
%               5  - Deep interneurons (L5)
%               6  - Thalamic projection neurons (pyramid) (L6)
%               7  - Reticular cells (Thal)
%               8  - Thalamo-cortical relay cells (Thal)
%
%
%        state: 1 V  - voltage
%               2 gE - conductance: AMPA   (excitatory)
%               3 gI - conductance: GABA-A (inhibitory)
%               4 gN - conductance: NMDA   (excitatory)
%               5 gB - conductance: GABA-B (inhibitory)
%
%      outputs: f = model states as a vector - hint: spm_unvec(f,M.x) 
%               J = system Jacobian - dfdx
%               Q = delay operator  - Q = inv(1 - D.*dfdx)*f(x(t))
%               D = states delay matrix
%
%
%
% Alexander Shaw 2019: ShawA10@cardiff.ac.uk
%
% Notes, changes, updates:
%
%
%
%--------------------------------------------------------------------------



% Flag: include M- & H- channels on L6 TP & Thalamic Relay cells, or not
%--------------------------------------------------------------------------
IncludeMH = 1;
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
ns   = size(M.x,1);                      % number of sources
np   = size(M.x,2);                      % number of populations per source
nk   = size(M.x,3);                      % number of states per population
x    = reshape(x,ns,np,nk);              % hidden states 


% extrinsic connection strengths
%==========================================================================
 
% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
for i = 1:length( P.A )
    A{i}  = exp(P.A{i});
    AN{i} = exp(P.AN{i});
%     if i > 2
%         A{i}  = A{i}  *0;
%         AN{i} = AN{i} *0;
%     end
end

C     = exp(P.C); 
 

% detect and reduce the strength of reciprocal (lateral) connections
%--------------------------------------------------------------------------
for i = 1:length(A)
    L    = (A{i} > exp(-8)) & (A{i}' > exp(-8));
    A{i} = A{i}./(1 + 8*L);
end

            
% intrinsic connection strengths
%==========================================================================
G    = full(P.H);
G    = exp(G);

% rep G 
n_pop = 10;

G = kron(G,ones(n_pop));

% dispersion of within population effects (10 cells per pop)
for i   = 1:8
    ind = (i*n_pop) - (n_pop-1);
    dg  = G( ind:ind+(n_pop-1),ind:ind+(n_pop-1) )  ;
    %dg  = diag(diag(dg)) + triu(dg,-1) - triu(dg);
    dg = triu(dg) - triu(dg,2);
    G( ind:ind+(n_pop-1),ind:ind+(n_pop-1) ) = dg*(2*n_pop);
end



% if ~all(size(P.G)==np)
%     for i = 1:size(G,3)
%         % trial specific intrinsic effects !
%         Gtrial   = diag( (P.G));
%         G(:,:,i) = G(:,:,i) + Gtrial; 
%     end
% elseif all(size(P.G)==np)
%     % a full 8*8 connectivity for this trial
%     for i = 1:size(G,3)
%         G(:,:,i) = G(:,:,i) + P.G;
%     end
% end
% 



% connectivity switches
%==========================================================================
% 1 - excitatory spiny stellate cells (granular input cells)
% 2 - superficial pyramidal cells     (forward  output cells)
% 3 - inhibitory interneurons         (intrisic interneuons)
% 4 - deep pyramidal cells            (backward output cells)
% 5 - deep interneurons               
% 6 - thalamic projection pyramidal cells (with m- and h- currents)
% 7 - thalamic reticular cells
% 8 - thalamic relay cells (with m- and h- currents)
%
% Thalamic cells attached to different cortical regions (models) are laterally connected


% % extrinsic connections (F B) - from superficial and deep pyramidal cells
% %--------------------------------------------------------------------------
%       SP  DP  tp  rt  rc
SA   = [1   0   0   0   0;   %  SS
        0   1   0   0   0;   %  SP
        0   1   0   0   0;   %  SI
        1   0   0   0   0;   %  DP
        0   0   0   0   0;   %  DI
        0   0   0   0   0;   %  TP
        0   0   0   0   0;   %  rt
        0   0   0   0   0]/8;%  rc
    
% % extrinsic NMDA-mediated connections (F B) - from superficial and deep pyramidal cells
% %--------------------------------------------------------------------------    
SNMDA = [1   0   0   0   0;   %  SS
         0   1   0   0   0;   %  SP
         0   1   0   0   0;   %  SI
         1   0   0   0   0;   %  DP
         0   0   0   0   0;   %  DI
         0   0   0   0   0;   %  TP
         0   0   0   0   0;   %  rt
         0   0   0   0   0]/8;%  rc

% intrinsic connectivity switches
%--------------------------------------------------------------------------    
%   population: 1  - Spint stellates (L4)                : e
%               2  - Superficial pyramids (L2/3)         : e
%               3  - Inhibitory interneurons (L2/3)      : i
%               4  - Deep pyramidal cells (L5)           : e
%               5  - Deep interneurons (L5)              : i
%               6  - Thalamic projection neurons -L6     : e
%               7  - Reticular cells (Thal)              : i
%               8  - Thalamo-cortical relay cells (Thal) : e

GEa = zeros(8,8);
GIa = zeros(8,8);

% Excitatory (np x np): AMPA & NMDA
%--------------------------------------------------------------------------
% This is a simplified, predictive-coding friendly excitatory architecture
%           ss  sp  si  dp  di  tp  rt  rl   
GEa(1,:) = [0   0   0   0   0   2   0   2]/1;
GEa(2,:) = [4   0   0   0   0   0   0   0]/1;
GEa(3,:) = [4   4   0   0   0   0   0   0]/1; 
GEa(4,:) = [0   4   0   0   0   0   0   0]/1;
GEa(5,:) = [0   0   0   4   0   0   0   0]/1;
GEa(6,:) = [0   0   0   2   0   0   0   1/4]/1; % added RL->TP [Ghodrati 2017]
GEa(7,:) = [0   0   0   0   0   0   0   2]/1; 
GEa(8,:) = [0   0   0   0   0   2   0   0]/1;


%GEa = GEa .* ~eye(np);
GEa = GEa + eye(np/n_pop);
GEn = GEa;

% Inhibitory connections (np x np): GABA-A & GABA-B
%--------------------------------------------------------------------------
%           ss  sp  si  dp  di  tp  rt  rl
GIa(1,:) = [8   0   2   0   0   0   0   0 ];
GIa(2,:) = [0   16  16  0   0   0   0   0 ];
GIa(3,:) = [0   0   32  0   0   0   0   0 ];
GIa(4,:) = [0   0   0   8   8   0   0   0 ];
GIa(5,:) = [0   0   0   0   16  0   0   0 ];
GIa(6,:) = [0   0   0   0   8   8   0   0 ];
GIa(7,:) = [0   0   0   0   0   0   32  0 ];
GIa(8,:) = [0   0   0   0   0   0   8   32];

%GIa = GIa*2;
GIb      = GIa;


% rep these up to n-cells per pop
%--------------------------------------------------------------------------

GEa = kron(GEa,ones(n_pop));
GIa = kron(GIa,ones(n_pop));

% for i   = 1:8
%     ind = (i*n_pop) - (n_pop-1);
%     GIa( ind:ind+(n_pop-1),ind:ind+(n_pop-1) ) = ...
%         GIa( ind:ind+(n_pop-1),ind:ind+(n_pop-1) ) + ( .001*rand(10) );%.* randi([0 1],10);    
% end

GEa = GEa * 2;
GIa = GIa * 4;

GEn      = GEa;
GIb      = GIa;






if IncludeMH
    
    % M- & H- channel conductances (np x np) {L6 & Thal Relay cells only}
    %----------------------------------------------------------------------
    VM   = -70;                            % reversal potential m-channels          
    VH   = -30;                            % reversal potential h-channels 

    %GIm  = eye(8)*4/10;                    % local TP & RL expression only
    %GIm  = sparse([6 8],[6 8],1/10,8,8);
    
    GIm  = kron(eye(8),ones(n_pop)) * 4/10;
    Mh   = kron(diag(exp(P.Mh)),ones(n_pop));

    %GIh      = full(sparse([6 8],[6 8],1/10,8,8));
    GIh = kron(eye(8),ones(n_pop)) * 4/10;
    
    %Hh       = exp(P.Hh);
    Hh = kron(diag(exp(P.Hh)),ones(n_pop));
    
    KM    = (exp(-P.m)*1000/160) ;               % m-current opening + CV
    KH    = (exp(-P.h)*1000/100) ;               % h-current opening + CV
    h     = 1 - spm_Ncdf_jdw(x(:,:,1),-100,300); % mean firing for h-currents
end

% Channel rate constants [decay times]
%--------------------------------------------------------------------------
KE  = exp(-P.T(:,1))*1000/4;            % excitatory rate constants (AMPA)
KI  = exp(-P.T(:,2))*1000/16;           % inhibitory rate constants (GABAa)
KN  = exp(-P.T(:,3))*1000/100;          % excitatory rate constants (NMDA)
KB  = exp(-P.T(:,4))*1000/200;          % excitatory rate constants (NMDA)


% Trial effects on time constants: AMPA & NMDA only
if isfield(P,'T1')
    KE = KE + P.T1(1);
    KN = KN + P.T1(2);
end


% Voltages [reversal potentials] (mV)
%--------------------------------------------------------------------------
VL   = -70;                               % reversal  potential leak (K)
VE   =  60;                               % reversal  potential excite (Na)
VI   = -90;                               % reversal  potential inhib (Cl)
VR   = -40;                               % threshold potential (firing)
VN   =  10;                               % reversal Ca(NMDA)   
VB   = -100;                              % reversal of GABA-B

% membrane capacitances {ss  sp  ii  dp  di  tp   rt  rl}
%--------------------------------------------------------------------------
CV   = exp(P.CV).*      [128 32  32  128 64  128  256 64*8]/1000;  


CV = spm_vec(repmat(CV,[n_pop 1]))';


% leak conductance - fixed
%--------------------------------------------------------------------------
GL   = 1;          

% mean-field effects:
%==========================================================================

% neural-mass approximation to covariance of states: trial specific
%----------------------------------------------------------------------
Vx   = exp(P.S)*32;
Vx   = spm_vec(repmat(Vx,[n_pop 1]))';
if nargin < 5
    % compute only if not passed by integrator
    m    =     spm_Ncdf_jdw(x(:,:,1),VR,Vx);
end

%m = m*n_pop.^2;
%m = (m/n_pop) + rand(size(m));

% extrinsic effects
%--------------------------------------------------------------------------
a       = zeros(ns,5);
an      = zeros(ns,5); 
a(:,1)  = A{1}*m(:,2);                      % forward afference  AMPA
a(:,2)  = A{2}*m(:,4);                      % backward afference AMPA
a(:,3)  = A{3}*m(:,6);                      % thalamic projection pyramids
a(:,4)  = A{4}*m(:,7);                      % reticular AMPA
a(:,5)  = A{5}*m(:,8);                      % relay AMPA
an(:,1) = AN{1}*m(:,2);                     % forward afference  NMDA
an(:,2) = AN{2}*m(:,4);                     % backward afference NMDA
an(:,3) = AN{3}*m(:,6);                     % thalamic projection pyramids
an(:,4) = AN{4}*m(:,7);                     % reticular NMDA
an(:,5) = AN{5}*m(:,8);                     % relay NMDA

% Averge background activity and exogenous input
%==========================================================================
BE     = exp(P.E)*0.8;

% input(s)
%--------------------------------------------------------------------------
if isfield(M,'u')
      U =   u;%(:); % endogenous input
else; U = C*u;%(:); % exogenous input
end

% flow over every (ns x np) subpopulation
%==========================================================================
f     = x;

% Thalamo-cortical flow [eq. motion] over modes, populations, states...
%--------------------------------------------------------------------------
for i = 1:ns
   
        % allow switching of thalamus on/off
        %HasThal = M.HasThal(i);
        
    
        % input scaling: Main input = RL->SS, but weak ?--> RL, SS, SP & DP
        %------------------------------------------------------------------
        if any(full(U(:))) && size(U,1) >= i
            dU = u(1)*( C(i,:).*[1 1/64 1/128 1/128] );
        else
            dU = [0 0 0 0];
        end
        
        % intrinsic coupling - parameterised
        %------------------------------------------------------------------
        
%         % sep excitatory to inhibitory self-mods
%         G0(:,:,i) = G(:,:,i) .* ~eye(8);
%         G0(:,:,i) = G0(:,:,i) + diag(P.GE);
%         
%         E      = ( G0(:,:,i).*GEa)*m(i,:)';
%         ENMDA  = ( G0(:,:,i).*GEn)*m(i,:)';
        
        E      = ( G(:,:,i).*GEa)*m(i,:)'; % AMPA currents
        ENMDA  = ( G(:,:,i).*GEn)*m(i,:)'; % NMDA currents
        I      = ( G(:,:,i).*GIa)*m(i,:)'; % GABA-A currents
        IB     = ( G(:,:,i).*GIb)*m(i,:)'; % GABA-B currents
        
        if IncludeMH
            
            % intrinsic coupling - non-parameterised: intrinsic dynamics
            %--------------------------------------------------------------
            Im     = (Mh(:,:).*GIm)*m(i,:)'; % M currents
            Ih     =           GIh *h(i,:)'; % H currents
        end
        
        % extrinsic coupling (excitatory only) and background activity
        %------------------------------------------------------------------
        SAa = spm_vec( repmat( (SA*a(i,:)')', [n_pop 1]) )';
        E     = (E     +  BE  + SAa' )*2;
        
        SNMDAa = spm_vec( repmat( (SNMDA*an(i,:)')', [n_pop 1]) )';
        ENMDA = (ENMDA +  BE  + SNMDAa')*2;
      
        % and exogenous input(U): 
        %------------------------------------------------------------------
        %input_cell        = [8 1 2 4];
        
        cell_id = spm_vec(repmat(1:8,[n_pop 1]))';
        
        input_cell = [ find(cell_id==8) ...
                       find(cell_id==1) ...
                       find(cell_id==2) ...
                       find(cell_id==4) ];
        
        E(input_cell)     = E(input_cell)     +spm_vec( repmat(dU'/n_pop,[n_pop,1]) );%    +dU';
        ENMDA(input_cell) = ENMDA(input_cell) +spm_vec( repmat(dU'/n_pop,[n_pop,1]) );   % +dU';
        
%         try   E(input_cell)          = E(input_cell)         +U';%  + U(i,1);
%               ENMDA(input_cell)      = ENMDA(input_cell)     +U';%  + U(i,1);
%         catch E(input_cell)          = E(input_cell)         +U';%  + U(1,1);
%               ENMDA(input_cell)      = ENMDA(input_cell)     +U';%  + U(1,1);
%         end
        
        % Voltage equation
        %==================================================================
        if ~IncludeMH
            
          f(i,:,1) =         (GL*(VL - x(i,:,1))+...
                       x(i,:,2).*(VE - x(i,:,1))+...
                       x(i,:,3).*(VI - x(i,:,1))+...
                       x(i,:,5).*(VB - x(i,:,1))+...
                       x(i,:,4).*(VN - x(i,:,1)).*mg_switch(x(i,:,1)))./CV;
            
        elseif IncludeMH
            
          f(i,:,1) =         (GL*(VL - x(i,:,1))+...
                       x(i,:,2).*(VE - x(i,:,1))+...
                       x(i,:,3).*(VI - x(i,:,1))+...
                       x(i,:,5).*(VB - x(i,:,1))+...
                       x(i,:,6).*(VM - x(i,:,1))+...
                       x(i,:,7).*(VH - x(i,:,1))+...
                       x(i,:,4).*(VN - x(i,:,1)).*mg_switch(x(i,:,1)))./CV;
        end
                   
        % Conductance equations
        %==================================================================        
        f(i,:,2) = (E'     - x(i,:,2)).*KE(i,:);
        f(i,:,3) = (I'     - x(i,:,3)).*KI(i,:);
        f(i,:,5) = (IB'    - x(i,:,5)).*KB(i,:);
        f(i,:,4) = (ENMDA' - x(i,:,4)).*KN(i,:);
        
        if IncludeMH
            f(i,:,6) = (Im'    - x(i,:,6)).*(KM(i,:) );
            f(i,:,7) = (Ih'    - x(i,:,7)).*(KH(i,:) );
        end

        % Restrict the states flow rate by the conduction rate of the population
        % c.f. synaptic delays + conduction delays
        %------------------------------------------------------------------
        DV       = 1./[1 1 1 2.2 1 2 8 8]; 
        DV       = 1./[2 1 1 2.2 1 2 1 2]; 
        DV       = DV.*exp(P.TV);
        
        DV = spm_vec( repmat(DV,[n_pop,1]) )';
        
        f(i,:,2) = f(i,:,2) .* DV;  % AMPA
        f(i,:,3) = f(i,:,3) .* DV;  % GABA-A
        f(i,:,4) = f(i,:,4) .* DV;  % NMDA
        f(i,:,5) = f(i,:,5) .* DV;  % GABA-B
        
        if IncludeMH
            f(i,:,6) = f(i,:,6) .* DV;  % M
            f(i,:,7) = f(i,:,7) .* DV;  % H
        end 
                
end


% vectorise equations of motion
%==========================================================================
f = spm_vec(f);


 
if nargout < 2, return, end

% Only compute Jacobian (gradients) if requested
%==========================================================================
J = spm_cat(spm_diff(M.f,x,u,P,M,1));

%J = sparse(jaco_dcm(M.f,x,u,P,M,1));


if nargout < 3, return, end

% Only compute Delays if requested
%==========================================================================
% Delay differential equations can be integrated efficiently (but 
% approximately) by absorbing the delay operator into the Jacobian
%
%    dx(t)/dt     = f(x(t - d))
%                 = Q(d)f(x(t))
%
%    J(d)         = Q(d)df/dx
%--------------------------------------------------------------------------
% [specified] fixed parameters
%--------------------------------------------------------------------------

D  = [.6 16];
d  = -D.*full(exp(P.D(1:2)))/1000;
Sp = kron(ones(nk,nk),kron( eye(np,np),eye(ns,ns)));  % states: same pop.
Ss = kron(ones(nk,nk),kron(ones(np,np),eye(ns,ns)));  % states: same source

% Thalamo cortical interactions: ~80ms round trip: 20 ms T->C, 60 ms C->T
%--------------------------------------------------------------------------
Tc              = zeros(np,np);
Tc([7 8],[1:6]) = 60  * exp(P.D0(1)); % L6->thal
Tc([1:6],[7 8]) = 20  * exp(P.D0(2)); % thal->ss

Tc = -Tc / 1000;
Tc = Tc .* ~~(GEa | GIa);

Tc = kron(ones(nk,nk),kron(Tc,eye(ns,ns)));


% intrisc delays
%ID = [.6 .2 .1 .6 .2 .6 .2 .6];
%ID = -ID.*exp(P.ID)/1000;
%ID = kron(ones(nk,nk),kron(diag(ID),eye(ns,ns)));

%Tc = Tc + ID;

% if ~isfield(M,'HasThal')
%     M.HasThal = ones(ns,1);
% end
% 
% if all(M.HasThal)
%     Tc = kron(ones(nk,nk),kron(Tc,eye(ns,ns)));
% else
%     Tc = kron(ones(nk,nk),kron(Tc,diag(M.HasThal)));
% end

% Mean intra-population delays, inc. axonal etc. Seem to help oscillation
%--------------------------------------------------------------------------
Dp = ~Ss;                            % states: different sources
Ds = ~Sp & Ss;                       % states: same source different pop.
%Ds = Ds.*(~(Ds & Tc));              % remove t-c and c-t from intrinsic
D  = d(2)*Dp + d(1)*Ds + Tc  ;       %+ Dself;% Complete delay matrix

% Implement: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = spm_inv(speye(length(J)) - D.*J);

