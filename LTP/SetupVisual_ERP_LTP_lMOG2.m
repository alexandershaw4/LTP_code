function SetupVisual_ERP_LTP_lMOG2(s)

addpath('~/code/'); amodel.start;

% Data & Design
%--------------------------------------------------------------------------
Data.Datasets     = 'Datas_LTP.txt';  % textfile list of LFP SPM datasets (.txt)
Data.Design.X     = [-1 0 1]';              % std/dev
Data.Design.name  = {'Undefined'};         % condition names

% NON-TETANISED:
Data.Design.tCode = [4 5 6];             % condition codes in SPM

% TETANISED
%Data.Design.tCode = [1 2 3];             % condition codes in SPM

Data.Design.Ic    = [1];             % channel indices
Data.Design.Sname = {'V1'};         % channel (node) names
Data.Prefix       = 'FixNew_-101_';      % outputted DCM prefix
Data.Datasets     = ReadverifyDatasets(Data.Datasets);

% Model space - T = ns x ns, where 1 = Fwd, 2 = Bkw
%--------------------------------------------------------------------------

% region order: L MOG, R MOG, L ITG, R ITG, L SPL, R SPL

%    LM  
T = [0 ]; % L Mog
F = (T==1);
B = (T==2);
C = [1]';      % inputs
L = sparse(1,1); 

nt = size(Data.Design.X,1); % num trials modelling
nb = size(Data.Design.X,2); % num contrasts/betas

% NEW: APPEAND THE MEAN DATASET AND FIT THAT FIRST
%--------------------------------------------------------------------------
Data.Datasets = ['LFPs_PlasticityMean.mat'  Data.Datasets];


% Set up, over subjects
%--------------------------------------------------------------------------
%for s = 1:length(Data.Datasets)
    
    % Data Naming & Design Matrix
    %----------------------------------------------------------------------
    DCM          = [];
    [fp fn fe]   = fileparts(Data.Datasets{s});
    DCM.name     = [Data.Prefix fn fe];
    
    DCM.xY.Dfile = Data.Datasets{s};  % original spm datafile
    Ns           = length(F);         % number of regions / modes
    DCM.xU.X     = Data.Design.X;     % design matrix
    DCM.xU.name  = Data.Design.name;  % condition names
    tCode        = Data.Design.tCode; % condition index (in SPM)
    DCM.xY.Ic    = Data.Design.Ic;    % channel indices
    DCM.Sname    = Data.Design.Sname; % channel names
        
    % Extrinsic Connectivity - Model Space
    %----------------------------------------------------------------------
    DCM.A{1} = F;
    DCM.A{2} = B;
    DCM.A{3} = L;
    DCM.B{1} = DCM.A{1} | DCM.A{2};
    DCM.B(2:length(DCM.xU.X)) = DCM.B;
    DCM.C    = C;
    
    % Function Handles
    %----------------------------------------------------------------------
    DCM.M.f  = @atcm.tc_dev;
    DCM.M.IS = 'spm_gen_erp';
    %DCM.options.SpecFun = @atcm.fun.Afft;
    
    % Print Progress
    %----------------------------------------------------------------------
    fprintf('Running Dataset %d / %d\n',s,length(Data.Datasets));
    
    % Prepare Data
    %----------------------------------------------------------------------
    fprintf('Preparing ERP data\n');

    DCM.M.U            = sparse(diag(ones(Ns,1)));  %... ignore [modes]
    DCM.options.trials = tCode;                     %... trial code [GroupDataLocs]
    DCM.options.Tdcm   = [0 350];                   %... peristimulus time
    DCM.options.Fdcm   = [1 50];                    %... frequency window
    DCM.options.D      = 1;                         %... downsample
    DCM.options.han    = 1;    %1                     %... apply hanning window
    DCM.options.h      = 0;%4;                         %... number of confounds (DCT)
    DCM.options.DoData = 1;                         %... leave on [custom]
    %DCM.options.Bdcm   = [-200 0];                  %... baseline times [new!]
    %DCM.options.Fltdcm = [1 30];                    %... bp filter [new!]

    DCM.options.analysis      = 'ERP';              %... analyse type
    DCM.xY.modality           = 'LFP';              %... ECD or LFP data? [LFP]
    DCM.options.spatial       = 'LFP';              %... spatial model [LFP]
    DCM.options.model         = 'tcm';              %... neural model
    DCM.options.Nmodes        = length(DCM.M.U);    %... number of modes

    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.spm_dcm_erp_data(DCM);
    
    % polarity flip the ERPs for Rachael's data
    if s > 1
        % only flip subject data, not the mean: alrwady flipped
        for ip = 1:length(DCM.xY.y)
            DCM.xY.y{ip} = DCM.xY.y{ip} * - 1;
        end
    end

    DCM.options.DATA = 1 ;      
    
    % Subfunctions
    %----------------------------------------------------------------------
    fprintf('Setting prior parameters\n');

    DCM = atcm.parameters(DCM,Ns);       % gets latet priors for tc nmm 
    
    % new extras: thalamo-cortical and cortico-thalamic delays
    DCM.M.pE.D0   = [0 0];
    DCM.M.pC.D0   = [0 0]+1/8;
    DCM.M.HasThal = [1];       % switch for including thalamic 
    DCM.M.pC.R = [1 1 1]/8;
    
    % new - trial specific intrinsic self mods
    DCM.M.pE.G = zeros(8,8);
    DCM.M.pC.G = zeros(8,8);
    
    Mods = [1 0 1 0 0 0 0 1;
            1 1 1 0 0 0 0 0;
            1 1 1 0 0 0 0 0;
            0 1 0 0 0 0 0 0;
            0 0 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0;
            0 0 0 0 0 0 0 0;
            0 0 0 0 0 1 0 0]/8;
    DCM.M.pC.G = Mods;
    
    % rep over 3 trials
    %DCM.M.pE.G(:,:,2) = DCM.M.pE.G(:,:,1);
    %DCM.M.pE.G(:,:,3) = DCM.M.pE.G(:,:,1);
    %DCM.M.pC.G(:,:,2) = DCM.M.pC.G(:,:,1); 
    %DCM.M.pC.G(:,:,3) = DCM.M.pC.G(:,:,1); 

    % new - trial specific AMPA and NMDA TCs
    DCM.M.pE.T1 = zeros(nb,2);
    DCM.M.pC.T1 = ones (nb,2)/8;
    
    for i = 1:nb
        DCM.M.pE.B{i}  = 1;
        DCM.M.pC.B{i}  = 1/8;
        DCM.M.pE.BN{i} = 1;
        DCM.M.pC.BN{i} = 1/8;
    end
    
    
    % intrinsics (base)
    
    H = zeros(8,8);
    H = [1 0 1 0 0 0 0 1;% ss
         1 1 1 0 0 0 0 0;% sp
         1 1 1 0 0 0 0 0;% si
         0 1 0 1 1 0 0 0;% dp
         0 0 0 1 1 0 0 0;% di
         0 0 0 0 1 1 0 1;% tp
         0 0 0 0 0 0 1 1;% rt
         0 0 0 0 0 1 0 1]/8;% rl
    DCM.M.pC.H=H;
    
    DCM.M.gC.J([1 2 4 6])=1/8;
    DCM.M.pC.D(1) = 1/16;
        
    fprintf('Completing model specification\n');
    DCM = atcm.complete_erp(DCM);         % complete network specification
            
    if s > 1
        % load the mean fit and use its posteriors as priors
        fprintf('Loading mean model to use as prior on individual fits\n');
        meanmodel = load([Data.Prefix 'LFPs_PlasticityMean']);
        DCM.M.pE = meanmodel.DCM.Ep;
        DCM.M.gE = meanmodel.DCM.Eg;
    end
    
    
    % load the prior variances file
    %fprintf('Reading prior variances file\n');
    %XV = load('PriorVarsNew.mat','V');
    %DCM.M.pC = XV.V;
    
    % load some new priors i saved
    %fprintf('Reading prior parameters file\n');
    %PX = load('NewAugPriors.mat');
    %DCM.M.pE = PX.pE;
    
    % more
    %load('SomePriors_Aug5.mat');
    %DCM.M.pE = P0;
    %DCM.M.pC = PC;
    %DCM.M.gE = G0;
    %DCM.M.gC = GC;
    
    % show inversion plots (0=yes, 1=no)
    %---------------------------------------------------------------------
    DCM.M.nograph = 0;
    
    % only integrate over the same time period as the real data of interest
    %---------------------------------------------------------------------
    fprintf('Setting simulation time period & sample rate\n');
    DCM.M.sim.pst = DCM.xY.pst;
    DCM.M.sim.dt  = DCM.xY.dt;
    
    
    % Or, Fit the model using DCM
    %---------------------------------------------------------------------
    fprintf('Inverting\n');
    DCM = atcm.optim.dcm_erp_invert(DCM);
    close; drawnow;
    
    AllModels(s) = DCM;
    
    
end