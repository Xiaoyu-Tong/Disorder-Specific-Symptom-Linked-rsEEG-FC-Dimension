clear

loadClinicalData

%% test consistency using CV
[subName,IAA,IBB]=unique(subName,'stable');
templabel=[]; tempscore=[]; IDX_test=[];
nRound = 5;
nf=10;
ncomp = 24;
nPC = 300;

allMeasure = allMeasure(target>=7,:);% target>7 means ASD patients here

% standardization of FC and clinical scores
target_all = allMeasure - repmat(mean(allMeasure),size(allMeasure,1),1);
target_all= target_all./repmat(std(target_all),size(target_all,1),1);
% target_all = allMeasure;

max_iter = 10000;
load(['CPCA300_allPatients_' EEGband '.mat']) % load the .mat file you saved with CPCA_allPatients.py

% alpha = 1; % the contrastive parameter that you are investigating
iAlpha = 1; % align this with the index of the current alpha in your list
Fe_all = Fe_testArchive(:,:,iAlpha); % retrieve the corresponding contrastive features, edit as needed

% ***
% DO NOT further normalize your contrastive connectivity features -- the
% script below automatically deals with the variance and de-mean is strongly
% discouraged because it violates the overall linearity of the
% transformation, thus breaking the linear interpretability of this system.

rng('default')

% fine-tune the sparsity parameters for your data. The sparsity parameters
% for contrastive features and clinical measures may be different.
log_lambda1 = -0.4; % fine-tune as needed
lambda1 = 10^log_lambda1;
log_lambda2 = -0.2; % fine-tune as needed
lambda2 = 10^log_lambda2;

r_test = zeros(nRound,ncomp);
p_test = zeros(nRound,ncomp);
WxArchive = zeros(nPC,ncomp,nf);
WyArchive = zeros(ncomp,ncomp,nf);

% Construct reference dimensions
[Wx_ref,Wy_ref,R_ref]=SparseCCA(Fe_all,target_all,ncomp,max_iter,lambda1,lambda2);
% Add sorting here
[val,idx]=sort(R_ref,'descend');
Wy_ref = Wy_ref(:,idx);
Wx_ref = Wx_ref(:,idx);
R_ref = val;
X = Fe_all;
Y = target_all;
X=X./repmat(std(X),size(X,1),1);
Y=Y-repmat(mean(Y),size(Y,1),1);
Y=Y./repmat(std(Y),size(Y,1),1);
V = Y*Wy_ref;
U = X*Wx_ref;

U_test = zeros(length(target_ASD),ncomp,nRound);
V_test = zeros(length(target_ASD),ncomp,nRound);
for iRound = 1:nRound
    cvidx=crossvalind('Kfold',1:size(FC,1),nf);
    for cv=1:nf
        idx_test=find(cvidx==cv); Fe_test=Fe_all(idx_test,:);
        idx_train=find(cvidx~=cv);
        Fe_train=Fe_all(idx_train,:); target_train=target_all(idx_train,:);
        target_test=allMeasure(idx_test,:);
   
        [Wx,Wy,R]=SparseCCA_PMD_noDevariance(Fe_train,target_train,ncomp,max_iter,lambda1,lambda2);
    % intermediate
        r_temp = zeros(ncomp,ncomp);
        rSign_temp = zeros(ncomp,ncomp);
        for i = 1:ncomp
            for j = 1:ncomp
%                             r_temp(i,j) = Wx_ref(:,i)' * Wx(:,j);
                r_temp(i,j) = Wy_ref(:,i)' * Wy(:,j);
                rSign_temp(i,j) = Wy_ref(:,i)' * Wy(:,j);
            end
        end
        idx_temp = zeros(ncomp,1);
        cTemp = zeros(1,ncomp);
        for i = 1:ncomp 
            [~,temp] = max(abs(r_temp(i,:)));
            idx_temp(i) = temp(1);
            cTemp(i) = sign(rSign_temp(i,idx_temp(i)));
%                         cTemp(i) = 1;
        end
        Wx = Wx(:,idx_temp) .* repmat(cTemp,nPC,1);
        Wy = Wy(:,idx_temp) .* repmat(cTemp,ncomp,1);
        R = R(idx_temp);
        % intermediate end
        WxArchive(:,:,(iRound-1)*nf+cv) = Wx;
        WyArchive(:,:,(iRound-1)*nf+cv) = Wy;
        % CV evaluate
        X = Fe_test;
        Y = target_test;
        X=X./repmat(std(Fe_train),size(X,1),1);
        Y=Y-repmat(mean(target_train),size(Y,1),1);
        Y=Y./repmat(std(target_train),size(Y,1),1);
        U_testTemp = X*Wx;
        V_testTemp = Y*Wy;
        rFold(:,(iRound-1)*nf+cv) = diag(corr(U_testTemp,V_testTemp));
        U_test(idx_test,:,iRound) = U_testTemp;
        V_test(idx_test,:,iRound) = V_testTemp;
    end
    for i = 1:ncomp
        [r_test(iRound,i),p_test(iRound,i)] = corr(U_test(:,i,iRound),V_test(:,i,iRound));%,'type','spearman');
    end
end



