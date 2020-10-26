
cd('/Users/Alex/Dropbox/tcm/LTP');
files = dir('ReInv*'); files = {files.name}';

% Loop over and unpack all of the models
%--------------------------------------------------------------------------
for i = 1:length(files)
    clear DCM; load(files{i});
    
    % unpack parameters
    p(i)  = atcm.fun.unpack_parameters(DCM,DCM.Ep);
    
    % unpack model predictions and error
    mod.pre(i,:)    = spm_vec(DCM.H{1});
    mod.early(i,:)  = spm_vec(DCM.H{2});
    mod.late(i,:)   = spm_vec(DCM.H{3});
    
    data.pre(i,:)   = spm_vec(DCM.H{1}) + spm_vec(DCM.R{1});
    data.early(i,:) = spm_vec(DCM.H{2}) + spm_vec(DCM.R{2});
    data.late(i,:)  = spm_vec(DCM.H{3}) + spm_vec(DCM.R{3});
    
end

% Get sampletimes
pst = DCM.xY.pst;

% condition the matrices - remove NaNs for averaging
cmod = spm_vec(mod);   cmod(isnan(cmod))=0;   mod = spm_unvec( full(cmod), mod);
cdata = spm_vec(data); cdata(isnan(cdata))=0; data = spm_unvec( full(cdata), data);


% note: p(1) is the parameter structure for subject 1, etc.



% Make ERP plots for each condition
%--------------------------------------------------------------------------
figure;
subplot(1,3,1);
kPlotShadedAverage(pst,data.pre,[0 0 1],[0 0 1],.6); hold on;
kPlotShadedAverage(pst,mod.pre,[1 0 0],[1 0 0],.6);  hold off;
title('Pre');ylim([-1.5 1.1]);

subplot(1,3,2);
kPlotShadedAverage(pst,data.early,[0 0 1],[0 0 1],.6); hold on;
kPlotShadedAverage(pst,mod.early,[1 0 0],[1 0 0],.6);  hold off;
title('Early');ylim([-1.5 1.1]);

subplot(1,3,3);
kPlotShadedAverage(pst,data.late,[0 0 1],[0 0 1],.6); hold on;
kPlotShadedAverage(pst,mod.late,[1 0 0],[1 0 0],.6);  hold off;
title('Late');ylim([-1.5 1.1]);

set(findall(gcf,'-property','FontSize'),'FontSize',18)


% Reconstruct the trial-specific intrinsic connectivity parameters for each
% of the three conditions:
%--------------------------------------------------------------------------
for s = 1:20
    Intr_Pre(s,:,:)  = p(s).Intrinsic + ( DCM.xU.X(1) * p(s).G(:,:,1) );
    Intr_Early(s,:,:)= p(s).Intrinsic + ( DCM.xU.X(2) * p(s).G(:,:,2) );
    Intr_Late(s,:,:) = p(s).Intrinsic + ( DCM.xU.X(3) * p(s).G(:,:,3) );
end

% These are 8x8, with cells in this order:
%       L4 stellate, L2/3 pyra, L2/3 inter, L5 pyra, L5 inter L6 pyra, thal ret, thal relay
pops = {'ss'         'sp'       'si'        'dp'     'di'     'tp'     'rt'      'rl'};


% so squeeze(nanmean(Intr_Pre,1)) is the group mean connectivity for pre

% to plot: 
% imagesc( squeeze(nanmean(Intr_Pre,1)) )
% set(gca,'ytick',1:8,'yticklabels',pops);
% set(gca,'xtick',1:8,'xticklabels',pops);


% Extract and quantify the AMPA and NMDA current contributions to the ERP
%-------------------------------------------------------------------------
for i = 1:length(files)
    clear DCM; load(files{i});
    
    states_pre(i,:,:)   = DCM.x{1};
    states_early(i,:,:) = DCM.x{2};
    states_late(i,:,:)  = DCM.x{3};
end

states_pre(isnan(states_pre)) = 0;
states_early(isnan(states_early))=0;
states_late(isnan(states_late))=0;

% these states matrices are: subject x time x model state
% where the states are membrane voltages and receptor currents for each
% population (8 pops, 7 states each = 56 states)
% states are oredered 8xvoltage, 8xAMPA, 8*GABAA, 8*NMDA, 8*GABAB,
% 8*m_current, 8*h_current, where 8 reflect the populations, in the order
% of 'pops' above

% Extract the AMPA & NMDA currents
Pre_AMPAs   = squeeze(states_pre(:,:,9:16));
Early_AMPAs = squeeze(states_early(:,:,9:16));
Late_AMPAs  = squeeze(states_late(:,:,9:16));

Pre_NMDAs   = squeeze(states_pre(:,:,17:24));
Early_NMDAs = squeeze(states_early(:,:,17:24));
Late_NMDAs  = squeeze(states_late(:,:,17:24));


% plot the ERPs with the currents:
%---------------------------------
% normalise so they're on the same scale [nope]
f0 = @(x) ( x - min(x) )./ ( max(x) - min(x) );

afigure, % AMPA
subplot(2,3,1); plot(pst,mean(data.pre)) ; title('Pre ERP');
subplot(2,3,4); plot(pst,squeeze(mean(Pre_AMPAs,1)) )
title('Pre: AMPA currents');

subplot(2,3,2); plot(pst,mean(data.early)) ; title('Early ERP');
subplot(2,3,5); plot(pst,squeeze(mean(Early_AMPAs,1)) )
title('Early: AMPA currents');

subplot(2,3,3); plot(pst,mean(data.late)) ; title('Late ERP');
subplot(2,3,6); plot(pst,squeeze(mean(Late_AMPAs,1)) )
title('Late: AMPA currents');legend(pops)


afigure, % NMDA
subplot(2,3,1); plot(pst,mean(data.pre)) ; title('Pre ERP');
subplot(2,3,4); plot(pst,squeeze(mean(Pre_NMDAs,1)) )
title('Pre: NMDA currents');

subplot(2,3,2); plot(pst,mean(data.early)) ; title('Early ERP');
subplot(2,3,5); plot(pst,squeeze(mean(Early_NMDAs,1)) )
title('Early: NMDA currents');

subplot(2,3,3); plot(pst,mean(data.late)) ; title('Late ERP');
subplot(2,3,6); plot(pst,squeeze(mean(Late_NMDAs,1)) )
title('Late: NMDA currents');legend(pops)


% make matrices for nxn correlations
for i = 1:88 % time points
    
    % AMPA--------------------------
    for j = 1:8 % pops
        
        x = data.pre(:,i);
        y = Pre_AMPAs(:,i,j);
        
        [r,p] = corr(x,y);
        
        pre_ampa_r0(i,j) = r;
        pre_ampa_p0(i,j) = p;
    end
    for j = 1:8 % pops
        
        x = data.early(:,i);
        y = Early_AMPAs(:,i,j);
        
        [r,p] = corr(x,y);
        
        early_ampa_r0(i,j) = r;
        early_ampa_p0(i,j) = p;
    end    
    for j = 1:8 % pops
        
        x = data.late(:,i);
        y = Late_AMPAs(:,i,j);
        
        [r,p] = corr(x,y);
        
        late_ampa_r0(i,j) = r;
        late_ampa_p0(i,j) = p;
    end    

    % NMDA--------------------------
    for j = 1:8 % pops
        
        x = data.pre(:,i);
        y = Pre_NMDAs(:,i,j);
        
        [r,p] = corr(x,y);
        
        pre_nmda_r0(i,j) = r;
        pre_nmda_p0(i,j) = p;
    end
    for j = 1:8 % pops
        
        x = data.early(:,i);
        y = Early_NMDAs(:,i,j);
        
        [r,p] = corr(x,y);
        
        early_nmda_r0(i,j) = r;
        early_nmda_p0(i,j) = p;
    end    
    for j = 1:8 % pops
        
        x = data.late(:,i);
        y = Late_NMDAs(:,i,j);
        
        [r,p] = corr(x,y);
        
        late_nmda_r0(i,j) = r;
        late_nmda_p0(i,j) = p;
    end    
    
    
end

% imagesc
afigure;
subplot(2,3,1);
imagesc(pre_ampa_r0'.^2*100);
set(gca,'ytick',1:8,'yticklabels',pops);
set(gca,'xtick',(pst(4:4:end))/pst(end) * length(pst),'xticklabels',pst(4:4:end));
rotateXLabels(gca,45);
xlabel('time (ms)');
title('% of ERP due to AMPA current (Pre)');
caxis([0 60]);%colorbar

% place another axis on top, for the ERP overlay
ax1_pos = get(gca,'position'); 
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
axes(ax2);
plot(pst,mean(data.pre),'w','linewidth',3);
set(ax2,'visible',0);


subplot(2,3,2);
imagesc(early_ampa_r0'.^2*100);
set(gca,'ytick',1:8,'yticklabels',pops);
set(gca,'xtick',(pst(4:4:end))/pst(end) * length(pst),'xticklabels',pst(4:4:end));
rotateXLabels(gca,45);
xlabel('time (ms)');
title('% of ERP due to AMPA current (Early)');
caxis([0 60]);%colorbar

% place another axis on top, for the ERP overlay
ax1_pos = get(gca,'position'); 
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
axes(ax2);
plot(pst,mean(data.early),'w','linewidth',3);
set(ax2,'visible',0);


subplot(2,3,3);
imagesc(late_ampa_r0'.^2*100);
set(gca,'ytick',1:8,'yticklabels',pops);
set(gca,'xtick',(pst(4:4:end))/pst(end) * length(pst),'xticklabels',pst(4:4:end));
rotateXLabels(gca,45);
xlabel('time (ms)');
title('% of ERP due to AMPA current (Late)');
caxis([0 60]);%colorbar

% place another axis on top, for the ERP overlay
ax1_pos = get(gca,'position'); 
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
axes(ax2);
plot(pst,mean(data.late),'w','linewidth',3);
set(ax2,'visible',0);


% NMDA
subplot(2,3,4);
imagesc(pre_nmda_r0'.^2*100);
set(gca,'ytick',1:8,'yticklabels',pops);
set(gca,'xtick',(pst(4:4:end))/pst(end) * length(pst),'xticklabels',pst(4:4:end));
rotateXLabels(gca,45);
xlabel('time (ms)');
title('% of ERP due to NMDA current (Pre)');
caxis([0 60]);%colorbar

% place another axis on top, for the ERP overlay
ax1_pos = get(gca,'position'); 
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
axes(ax2);
plot(pst,mean(data.pre),'w','linewidth',3);
set(ax2,'visible',0);

subplot(2,3,5);
imagesc(early_nmda_r0'.^2*100);
set(gca,'ytick',1:8,'yticklabels',pops);
set(gca,'xtick',(pst(4:4:end))/pst(end) * length(pst),'xticklabels',pst(4:4:end));
rotateXLabels(gca,45);
xlabel('time (ms)');
title('% of ERP due to NMDA current (Early)');
caxis([0 60]);%colorbar

% place another axis on top, for the ERP overlay
ax1_pos = get(gca,'position'); 
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
axes(ax2);
plot(pst,mean(data.early),'w','linewidth',3);
set(ax2,'visible',0);

subplot(2,3,6);
imagesc(late_nmda_r0'.^2*100);
set(gca,'ytick',1:8,'yticklabels',pops);
set(gca,'xtick',(pst(4:4:end))/pst(end) * length(pst),'xticklabels',pst(4:4:end));
rotateXLabels(gca,45);
xlabel('time (ms)');
title('% of ERP due to NMDA current (Late)');
caxis([0 60]);%colorbar

% place another axis on top, for the ERP overlay
ax1_pos = get(gca,'position'); 
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
axes(ax2);
plot(pst,mean(data.late),'w','linewidth',3);
set(ax2,'visible',0);


set(findall(gcf,'-property','FontSize'),'FontSize',18)
colormap(jet);

%pre_ampa_erp = [ full(mean(data.pre));
%                squeeze(mean(Pre_AMPAs,1))'] ;
%corr_mat(pre_ampa_erp)







