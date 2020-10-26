function DCM = parameters(DCM,Ns)
% collect initial (prior) parameters
%
% Need to rewrite this function by writing out the priors saved in matfile

%load('april23_priors.mat','pE','pC');

% load('april25_priors.mat','pE','pC');
% 
% pC.G = pC.G*0;
% 
% load('SomePrior','Pp');

if Ns == 1
    %load('May7.mat');

    %load NewDelayPriors.mat
    pth = fileparts(mfilename('fullpath'));
    load([pth '/AugSpectralPriors']);
    
    DCM.M.pE = PX.pE;
    DCM.M.pC = PX.pC;
    
else
    
    % NEEDS UPDATING FOR MULTINODE/REGION
    
    
    
    % Extrinsics: restructure adjacency matrices
    %----------------------------------------------------------------------
    A     = DCM.A;
    A{1}  = A{1} | A{3};                              % forward
    A{2}  = A{2} | A{3};                              % backward

    % [F] SP -> SS & DP
    pE.A{1} = A{1}*32 - 32;
    pC.A{1} = A{1}/8;

    % [B] DP -> SP & SI
    pE.A{2} = A{2}*32 - 32; 
    pC.A{2} = A{2}/8;

    % [Extrinsic T->C] (none)
    pE.A{3} = A{1}*32 - 32; 
    pC.A{3} = A{1}/8;

    % [Extrinsic T->T] rt -> rt & rl
    AA      = A{1} | A{2};
    
    pE.A{4} = AA*32 - 32; 
    pC.A{4} = AA/8;
    
    % [Extrinsic T->T] rl -> rt & rl
    pE.A{5} = AA*32 - 32; 
    pC.A{5} = AA/8;

    % NMDA = AMPA
    for i = 1:length(pE.A)
        pE.AN{i} = pE.A{i};
        pC.AN{i} = pC.A{i};
    end

    % Beta's: input-dependent scaling
    %----------------------------------------------------------------------
    B     = DCM.B;
    for i = 1:length(B)
        B{i} = ~~B{i};
        pE.B{i} = B{i} - B{i};
        pC.B{i} = B{i}/8;
    end

    for i = 1:length(B)
        B{i} = ~~B{i};
        pE.BN{i} = B{i} - B{i};
        pC.BN{i} = B{i}/8;
    end
    % Intrinsic (local) parameters, per region
    %----------------------------------------------------------------------
    ns = length(A{1}); % number of regions / nodes
    np = 8;            % number of populations per region
    nk = 7;            % number of states per population
    % exogenous inputs
    %----------------------------------------------------------------------
    C     = DCM.C;
    C     = ~~C;
    pE.C  = C*32 - 32;
    pC.C  = C/8;
    
    % new multi-input model
    %pE.C = repmat(pE.C,[1 4]);
    pC.C = repmat(pC.C,[1 4]);
    pE.C = repmat([-0.1164 -0.0025 -0.0017 -0.0011],[ns 1]);

    
    
    % Average Firing
    pE.S = zeros(1,8);
    pC.S = ones(1,8)*0.0625;
    
    % Average (Self) Excitation Per Population
    pE.G = zeros(ns,np);
    pC.G = zeros(ns,np); 
    pC.G(1:2) = 1;
    
    % Average Membrane Capacitance
    pE.CV =  [-0.0382646336383898
               0.00217382054723993
               0.0200769192441413
               0.0128455748594877
               0.0186176599910718
              -0.0369774537374279
               0.000783992915601112
               0.000912603133131219]';
    pC.CV = ones(1,8) * 0.0625;
    
    % Average Background Activity
    pE.E = 0;
    pC.E = 0;
    
    % Average Intrinsic Connectivity
    pE.H = repmat([...
    0.2448         0    0.0100         0         0         0         0    0.2698
   -0.0002    0.0330    0.1282         0         0         0         0         0
         0   -9.4951   -0.0575         0         0         0         0         0
         0         0         0   -8.7794         0         0         0         0
         0         0         0    0.0735   -0.0604         0         0         0
         0         0         0   -0.0298         0    0.0875         0    0.0006
         0         0         0         0         0         0    0.0186         0
         0         0         0         0         0    0.0955         0   -0.1479], [1 1 ns]);
    %pC.H = ~~pE.H / 8;
    pC.H =repmat([...
    0.0625         0    0.0625         0         0    0.0625         0    0.0625
    0.0625    0.0625    0.0625         0         0         0         0         0
    0.0625    0.0625    0.0625         0         0         0         0         0
         0    0.0625         0    0.0625    0.0625         0         0         0
         0         0         0    0.0625    0.0625         0         0         0
         0         0         0    0.0625    0.0625    0.0625         0    0.0625
         0         0         0         0         0         0    0.0625    0.0625
         0         0         0         0         0    0.0625    0.0625    0.0625], [1 1 ns]);
    
    % Receptor Time Constants (Channel Open Times)
    pE.T = repmat([0.3431   -0.0682    0.5655    3.1907],[ns,1]);
    pC.T =  ( ones(ns,4) / 8);
    
    % Parameters on input bump: delay, scale and width
    pE.R = [0 0 0];
    pC.R = [0 0 0];
    
    % Delays: states - intrinsic & extrinsic
    pE.D = [0 0];
    pC.D = [0 0];
    
    pE.D0 = [0 0];
    pC.D0 = [1 1]/8;
    
    % Dipole position (not in use)
    pE.Lpos = sparse(3,0);
    pC.Lpos = sparse(3,0);
    
    % Electrode gain
    pE.L = repmat(0.0534          ,[ns,1]);
    pC.L = repmat(64              ,[ns,1]);
    
    % Contributing states
    J           = zeros(1,np,nk) - 1000;
    J(:,[1 2 4 6],1)    = log([.2 .8 .2 .2]);
    J(isinf(J)) = -1000;
    pE.J        = spm_vec(J)';
    pC.J        = spm_vec(J)' * 0;
    
    % Noise Components - a, b & c
    pE.a = repmat([-1.7311;0],[1 ns]);
    pC.a = repmat([1/16;0]   ,[1 ns]);
    
    pE.b = [-0.0664;0];
    pC.b = [1/16   ;0];
    
    pE.c = repmat([-24.4382;-0.3475],[1 ns]);
    pC.c = repmat([0;0],                      [1 ns]);
    
    % Neuronal innovations
    pE.d = [ -1.1933
   -0.4832
    0.6237
   -0.2421
   -0.0032
   -0.0055
   -0.0004
    0.0047];
          
    pC.d = ones(length(pE.d),1) / 8;
       
    pE.h = repmat(0,  [ns 1]);
    pC.h = repmat(0,  [ns 1]);
    
    pE.m = repmat(0,  [ns 1]);
    pC.m = repmat(0,  [ns,1]);
    
    pE.TV = [-0.0173   -0.3861    0.3930   -0.0739    0.0054    0.2167   -0.1259    0.0094];
    pC.TV = ones(1,8)/16;
    
    pE.Mh = zeros(8,1);
    pC.Mh = zeros(8,1);
    
    pE.Hh = [0 0];
    pC.Hh = [0 0];
    
    % Pack pE & pC into M
    DCM.M.pE = pE;
    DCM.M.pC = pC;
    
end
    
    
end
