
cd '/cubric/scratch/sapas10/tcm/LTP'
files = dir('NewInt_LFPs*.mat'); files = {files.name}'; files = files(2:end);

for i = 1:length(files)
    clear DCM; load(files{i});
    
    % Data features
    [ymod,Yreal,pst] = atcm.fun.evalERP(DCM,DCM.Ep,DCM.Eg);
    
    ERP(i,:,:)  = full(spm_cat(Yreal));
    PRED(i,:,:) = full(spm_cat(ymod));
    
    % Get parameters
    p(i) = atcm.fun.unpack_parameters(DCM,DCM.Ep);
    
end

% plots
afigure;
for i = 1:20
    subplot(5,4,i);
    plot(pst,squeeze(ERP(i,:,:)),':',pst,squeeze(PRED(i,:,:)),'linewidth',2);

    if i == 1
        legend({'Pre [data]' 'Early [data]' 'Late [data]' ....
            'Pre [model]' 'Early [model]' 'Late [model]'});
    end
end

afigure;
for i = 1:20
    subplot(5,4,i);
    
    plot( (ERP(i,:)) ,':'); hold on;
    plot( (PRED(i,:)),'b','linewidth',2);
    
    [r,p] = corr( ERP(i,:)' , PRED(i,:)' );

    if i == 1
        legend({'Data' 'Model'});
    end
    title(sprintf('Var Expl %d%%',round(100*(r^2))));
end
set(findall(gcf,'-property','FontSize'),'FontSize',14)