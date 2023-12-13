% Representational similarity analysis for interpreting the superiority of contrastive features

clear
loadClinicalData

allMeasure = allMeasure(target>=7,:); % keep patient data
age_ASD = target_age(target>=7);

load(['CPCA300_allPatients_' EEGband '.mat']) % load your CPCA results
load(['compiledData_Russ31_CPCA_' EEGband '.mat']) % load your EEG data
CoV_patient = cov(Fe_ASD);
CoV_HC = cov(Fe_HC);
CoV_contrastive = CoV_patient - CoV_HC;

% alpha = 1; % the contrastive parameter that you are investigating
iAlpha = 1; % align this with the index of the current alpha in your list
Fe_contrastive = Fe_testArchive(:,:,iAlpha); % retrieve the corresponding contrastive features, edit as needed

nPC = 300;
[coeff,score,latent] = pca(Fe_ASD,'NumComponents',nPC);
[r_patient, p_patient] = corr(score,allMeasure);
[r_specific, p_specific] = corr(squeeze(Fe_contrastive,allMeasure);

score = score(:,1:nPC);
score_cpca = squeeze(Fe_contrastive);
score_cpca = score_cpca(:,1:nPC);

[coeff_shared,~,~] = pca(Fe_HC,'NumComponents',nPC);
coeff_shared = coeff_shared(:,1:min(nPC,size(coeff_shared,2)));
score_shared = Fe_ASD * coeff_shared;

nf = 10;
nTrial = 100;
nFe = 25; % Number of measures you want to examine. Here 25 = 24 (clinical measures) + 1 (age)
mdlFit_pca = zeros(nFe,nTrial);
mdlFit_cpca = zeros(nFe,nTrial);
mdlFit_shared = zeros(nFe,nTrial);

% RSA
nSub = length(target_ASD);
kSub = round(0.9*nSub); % take 90% of the subjects to run the analysis, to evaluate the stability of RSA results
score_archive = score;
score_cpca_archive = score_cpca;
score_shared_archive = score_shared;
for iTrial = 1:nTrial
    disp(iTrial)
    idx_train = randperm(nSub,kSub);
    RDM_standard = zeros(kSub,kSub);
    RDM_cpca = zeros(kSub,kSub);
    RDM_shared = zeros(kSub,kSub);
    RDM_scale = zeros(kSub,kSub,25);
    for iSub = 1:kSub
        for jSub = 1:kSub
            RDM_standard(iSub,jSub) = norm(score(idx_train(iSub),:)-score(idx_train(jSub),:));
            RDM_shared(iSub,jSub) = norm(score_shared(idx_train(iSub),:)-score_shared(idx_train(jSub),:));
            RDM_cpca(iSub,jSub) = norm(score_cpca(idx_train(iSub),:)-score_cpca(idx_train(jSub),:));
            for iFe = 1:nFe
                if iFe~=nFe
                    y = allMeasure(:,iFe);
                else
                    y = age_ASD;
                end
                RDM_scale(iSub,jSub,iFe) = norm(y(idx_train(iSub))-y(idx_train(jSub)));
            end
        end
    end
    RDM_standard = reshape(RDM_standard,[],1);
    RDM_shared = reshape(RDM_shared,[],1);
    RDM_cpca = reshape(RDM_cpca,[],1);
    RDM_scale = reshape(RDM_scale,length(RDM_cpca),25);
    for iFe = 1:nFe
        mdlFit_pca(iFe,iTrial) = corr(RDM_standard,RDM_scale(:,iFe),'type','spearman');
        mdlFit_cpca(iFe,iTrial) = corr(RDM_cpca,RDM_scale(:,iFe),'type','spearman');
        mdlFit_shared(iFe,iTrial) = corr(RDM_shared,RDM_scale(:,iFe),'type','spearman');
    end
end

% save your results
% save('interpretContrastiveFeatures_ALPHA.mat','mdlFit_pca','mdlFit_cpca','mdlFit_shared')

% run statistical test
p_pca = zeros(nFe,1);
p_shared = zeros(nFe,1);
for i = 1:nFe
    [p_pca(i),~] = signrank(mdlFit_cpca(i,:),mdlFit_pca(i,:));
    [p_shared(i),~] = signrank(mdlFit_cpca(i,:),mdlFit_shared(i,:));
end

% FDR correction
[q]=mafdr([p_pca;p_shared],'BHFDR',true);
q_pca = q(1:nFe);
q_shared = q(nFe+1:end);
q = [q_pca,q_shared];
