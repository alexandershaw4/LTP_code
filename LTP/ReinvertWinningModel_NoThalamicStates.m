% Re-run the winning (combination) LTP model, but WITHOUT the thalamic
% cells.
% Rather than writing a new model function, I have added an option within
% the integration function that removes the states corresponding to the
% thalamic populations.

cd '/cubric/scratch/sapas10/tcm/LTP/NewSept2019';

mods = dir('FixNew_-101-110*'); mods = {mods.name}';

for i = 1:length(mods)
    
    load(mods{i});
    
    DCM = rmfield(DCM,'Ep');
    DCM = rmfield(DCM,'Eg');
    
    DCM.M.HasThal = 0;
    
    DCM.name = ['NoThal_' DCM.name];
    
    fprintf('Inverting dataset %d/%d\n',i,length(mods));
    DCM = atcm.optim.dcm_erp_invert(DCM);
    close; drawnow;
        
    clear DCM
    
end