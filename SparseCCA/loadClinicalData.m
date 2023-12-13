
EEGband = 'ALPHA'; % specify the frequency band
% Load EEG data
load('ABC_Austim_EEG_powenv_Russ31ROIs.mat') % load your own EEG data
% the goal is to align clinical data with neuroimaging data and exclude
% subjects with missing data.

ROIConn=ROIConn_ALPHA;%band
nSub=size(ROIConn,3);
nROI=31;
Fe=zeros(nROI*(nROI-1)/2,nSub);
for iSub=1:nSub
    FeTemp=ROIConn(:,:,iSub);
    FeTemp=tril(FeTemp,-1);
    Fe(:,iSub)=FeTemp(FeTemp~=0);
end

% Process subject names
subName=cellfun(@(x) x(1:end-9),subjectID,'uniformoutput',false);
[subName,IAA,IBB]=unique(subName,'stable');
dateEEG=NaT(length(subjectID),1);
for i=1:length(subjectID)
    dateEEG(i)=datetime(subjectID{i}(end-7:end),'InputFormat','yyyyMMdd');
end

% % % session number
target_session = zeros(length(IBB),1);
for iSub = 1:length(subName)
    sessionIdx = find(IBB == iSub);
    target_session(sessionIdx) = 1:length(sessionIdx);
end

% load ASD severity data
data3=readtable('PhenotypicData/ados3_201201.txt');
targetID3=data3.subjectkey;
targetID=targetID3;
% targetID=[targetID1;targetID2;targetID3];
targetID_raw=targetID;
[targetID,IAA_target,IBB_target]=unique(targetID,'stable');
date3=data3.interview_date;
dateList=date3;
ados3=data3.scoresumm_overalltotal;
ados=ados3;
age=data3.interview_age;
sex = data3.sex;
sa = data3.scoresumm2_abtotal;
rrb = data3.scoresumm_adtotal;

[subName,IAA_ados,IBB_ados]=unique(subName,'stable');

% AIM
data_aim=readtable('PhenotypicData/aim01.txt');
targetID_aim = data_aim.subjectkey;
freq_rb = data_aim.freq_restricted_behaviors;
freq_comm = data_aim.frequency_com_language;
freq_sr = data_aim.freq_social_em_reciprocity;
freq_ab = data_aim.freq_odd_atypical_behavior;
impact_rb = data_aim.impact_restricted_behaviors;
impact_comm = data_aim.impact_com_language;
impact_sr = data_aim.impact_social_em_reciprocity;
impact_ab = data_aim.impact_odd_atypical_behavior;
freq_total = data_aim.frequency_domain_total;
impact_total = data_aim.impact_domain_total;
dateList_aim=data_aim.interview_date;
targetID_raw_aim=targetID_aim;
[targetID_aim,IAA_target_aim,IBB_target_aim]=unique(targetID_aim,'stable');

% target_aim = zeros(length(target),8);
target_aim_all = [freq_rb,freq_comm,freq_sr,freq_ab,impact_rb,impact_comm,impact_sr,impact_ab];

% SRS
data_srs=readtable('PhenotypicData/srs02.txt');
maleSRS = table2array(data_srs(:,20:24));
femaleSRS = table2array(data_srs(:,32:36));
target_srs = maleSRS;
target_srs(isnan(maleSRS(:,1)),:) = femaleSRS(isnan(maleSRS(:,1)),:);
targetID_srs = data_srs.subjectkey;

dateList_srs=data_srs.interview_date;
targetID_raw_srs=targetID_srs;
[targetID_srs,IAA_target_srs,IBB_target_srs]=unique(targetID_srs,'stable');

% vineland
data_vl=readtable('PhenotypicData/vinland301.txt');
target_vl = data_vl(:,542:550);
targetID_vl = data_vl.subjectkey;

dateList_vl=data_vl.interview_date;
targetID_raw_vl=targetID_vl;
[targetID_vl,IAA_target_vl,IBB_target_vl]=unique(targetID_vl,'stable');


% Eliminate redundant data for ADOS
IDidx=ones(length(IBB_target),1);
for sub=1:length(targetID)
    IBB_targetTemp=find(IBB_target==sub);
    uniDate=unique(dateList(IBB_targetTemp));
    for iDate=1:length(uniDate)
        idxTemp=find(dateList(IBB_targetTemp)==uniDate(iDate));
        if length(idxTemp)>1
            IDidx(IBB_targetTemp(idxTemp(2:end)))=NaN;
        end
    end
end
targetID_raw(isnan(IDidx))=[];
dateList(isnan(IDidx))=[];
ados(isnan(IDidx))=[];
age(isnan(IDidx))=[];
sex(isnan(IDidx))=[];
sa(isnan(IDidx))=[];
rrb(isnan(IDidx))=[];
% Eliminate redundant data for AIM
IDidx=ones(length(IBB_target_aim),1);
for sub=1:length(targetID)
    IBB_targetTemp=find(IBB_target_aim==sub);
    uniDate=unique(dateList_aim(IBB_targetTemp));
    for iDate=1:length(uniDate)
        idxTemp=find(dateList_aim(IBB_targetTemp)==uniDate(iDate));
        if length(idxTemp)>1
            IDidx(IBB_targetTemp(idxTemp(2:end)))=NaN;
        end
    end
end
targetID_raw_aim(isnan(IDidx))=[];
dateList_aim(isnan(IDidx))=[];
target_aim_all(isnan(IDidx),:)=[];

% Eliminate redundant data for SRS
IDidx=ones(length(IBB_target_srs),1);
for sub=1:length(targetID)
    IBB_targetTemp=find(IBB_target_srs==sub);
    uniDate=unique(dateList_srs(IBB_targetTemp));
    for iDate=1:length(uniDate)
        idxTemp=find(dateList_srs(IBB_targetTemp)==uniDate(iDate));
        if length(idxTemp)>1
            IDidx(IBB_targetTemp(idxTemp(2:end)))=NaN;
        end
    end
end
targetID_raw_srs(isnan(IDidx))=[];
dateList_srs(isnan(IDidx))=[];
target_srs(isnan(IDidx),:)=[];

% Eliminate redundant data for vinelind
IDidx=ones(length(IBB_target_vl),1);
for sub=1:length(targetID)
    IBB_targetTemp=find(IBB_target_vl==sub);
    uniDate=unique(dateList_vl(IBB_targetTemp));
    for iDate=1:length(uniDate)
        idxTemp=find(dateList_vl(IBB_targetTemp)==uniDate(iDate));
        if length(idxTemp)>1
            IDidx(IBB_targetTemp(idxTemp(2:end)))=NaN;
        end
    end
end
targetID_raw_vl(isnan(IDidx))=[];
dateList_vl(isnan(IDidx))=[];
target_vl(isnan(IDidx),:)=[];

% remove very early clinic visit, which is strange 
earliestEEGdate=min(dateEEG);
earliestMeasureDate=earliestEEGdate-13;
invalidDateIdx=dateList<earliestMeasureDate;
dateList(invalidDateIdx)=[];
ados(invalidDateIdx)=[];
sa(invalidDateIdx)=[];
age(invalidDateIdx)=[];
sex(invalidDateIdx)=[];
rrb(invalidDateIdx)=[];
targetID_raw(invalidDateIdx)=[];
[targetID,IAA_target,IBB_target]=unique(targetID_raw,'stable');
%AIM
invalidDateIdx=dateList_aim<earliestMeasureDate;
dateList_aim(invalidDateIdx)=[];
target_aim_all(invalidDateIdx,:)=[];
targetID_raw_aim(invalidDateIdx)=[];
[targetID_aim,IAA_target_aim,IBB_target_aim]=unique(targetID_raw_aim,'stable');

%SRS
invalidDateIdx=dateList_srs<earliestMeasureDate;
dateList_srs(invalidDateIdx)=[];
target_srs(invalidDateIdx,:)=[];
targetID_raw_srs(invalidDateIdx)=[];
[targetID_srs,IAA_target_srs,IBB_target_srs]=unique(targetID_raw_srs,'stable');
srs = target_srs;

%vinland
invalidDateIdx=dateList_vl<earliestMeasureDate;
dateList_vl(invalidDateIdx)=[];
target_vl(invalidDateIdx,:)=[];
targetID_raw_vl(invalidDateIdx)=[];
[targetID_vl,IAA_target_vl,IBB_target_vl]=unique(targetID_raw_vl,'stable');
vl = table2array(target_vl);

% subFirst=[]; % ADOS
EEG=[];
dateEEG_new=[];
sessionTemp = target_session;
target_session = [];
target=[];
target_diag = [];
target_age = [];
target_sex = [];
% target_aim = [];
target_sa=[];
target_rrb=[];
% target_srs = [];
% target_vl = [];
targetID_inUse=cell(0);
timeThres=14; % collection date of EEG and ados
for sub=1:length(targetID)
%     if isempty(IBB_target==sub)
%         continue;
%     end
    dateTemp=dateList(IBB_target==sub);
%     if length(dateTemp)==3
%         dateTemp=dateTemp(1:2);
%     end
%     if length(dateTemp)==4
%         error('Too many dates!! Invalid data?')
%     end
    idd=find(contains(subjectID,targetID{sub}));
    dateEEGtemp=dateEEG(idd);
    for iTemp=1:length(idd)
        if ~isempty(find(dateTemp-dateEEGtemp(iTemp)<timeThres&dateTemp-dateEEGtemp(iTemp)>-timeThres))
%             idxDateTemp=find(dateTemp==dateEEGtemp(iTemp));
            idxDateTemp=find(dateTemp-dateEEGtemp(iTemp)<timeThres&dateTemp-dateEEGtemp(iTemp)>-timeThres);
            temp=find(IBB_target==sub);
            EEG=[EEG,Fe(:,idd(iTemp))];
            target_session = [target_session;sessionTemp(idd(iTemp))];
            target=[target;ados(temp(idxDateTemp))];
            target_sa=[target_sa;sa(temp(idxDateTemp))];
            target_rrb=[target_rrb;rrb(temp(idxDateTemp))];
            target_age=[target_age;age(temp(idxDateTemp))];
            target_sex=[target_sex;sex(temp(idxDateTemp))];
            dateEEG_new = [dateEEG_new;dateEEG(idd(iTemp))];
%             target_diag=[target_diag;Diag(temp(idxDateTemp))];
            targetID_inUse=[targetID_inUse;subjectID(idd(iTemp))];
        end
    end
end
% subName=cellfun(@(x) x(1:end-9),targetID_inUse,'uniformoutput',false);



% subFirst=[]; % AIM
target_aim = [];
% target_srs = [];
% target_vl = [];
targetID_inUse_aim = cell(0);
timeThres=14; % collection date of EEG and ados
for sub=1:length(targetID_aim)
%     if isempty(IBB_target==sub)
%         continue;
%     end
    dateTemp=dateList_aim(IBB_target_aim==sub);
%     if length(dateTemp)==3
%         dateTemp=dateTemp(1:2);
%     end
%     if length(dateTemp)==4
%         error('Too many dates!! Invalid data?')
%     end
    idd=find(contains(subjectID,targetID_aim{sub}));
    dateEEGtemp=dateEEG(idd);
    for iTemp=1:length(idd)
        if ~isempty(find(dateTemp-dateEEGtemp(iTemp)<timeThres&dateTemp-dateEEGtemp(iTemp)>-timeThres))
%             idxDateTemp=find(dateTemp==dateEEGtemp(iTemp));
            idxDateTemp=find(dateTemp-dateEEGtemp(iTemp)<timeThres&dateTemp-dateEEGtemp(iTemp)>-timeThres);
            temp=find(IBB_target_aim==sub);
            target_aim=[target_aim;target_aim_all(temp(idxDateTemp(end)),:)];
%             target_diag=[target_diag;Diag(temp(idxDateTemp))];
            targetID_inUse_aim=[targetID_inUse_aim;subjectID(idd(iTemp))];
        end
    end
end
% subName=cellfun(@(x) x(1:end-9),targetID_inUse_aim,'uniformoutput',false);

% subFirst=[]; % SRS
target_srs = [];
% target_srs = [];
% target_vl = [];
targetID_inUse_srs = cell(0);
timeThres=14; % collection date of EEG and ados
for sub=1:length(targetID_srs)
%     if isempty(IBB_target==sub)
%         continue;
%     end
    dateTemp=dateList_srs(IBB_target_srs==sub);
%     if length(dateTemp)==3
%         dateTemp=dateTemp(1:2);
%     end
%     if length(dateTemp)==4
%         error('Too many dates!! Invalid data?')
%     end
    idd=find(contains(subjectID,targetID_srs{sub}));
    dateEEGtemp=dateEEG(idd);
    for iTemp=1:length(idd)
        if ~isempty(find(dateTemp-dateEEGtemp(iTemp)<timeThres&dateTemp-dateEEGtemp(iTemp)>-timeThres))
%             idxDateTemp=find(dateTemp==dateEEGtemp(iTemp));
            idxDateTemp=find(dateTemp-dateEEGtemp(iTemp)<timeThres&dateTemp-dateEEGtemp(iTemp)>-timeThres);
            temp=find(IBB_target_srs==sub);
            target_srs=[target_srs;srs(temp(idxDateTemp),:)];
%             target_diag=[target_diag;Diag(temp(idxDateTemp))];
            targetID_inUse_srs=[targetID_inUse_srs;subjectID(idd(iTemp))];
        end
    end
end
% subName=cellfun(@(x) x(1:end-9),targetID_inUse_srs,'uniformoutput',false);

% subFirst=[]; % Vineland
target_vl = [];
% target_vl = [];
% target_vl = [];
targetID_inUse_vl = cell(0);
timeThres=14; % collection date of EEG and ados
for sub=1:length(targetID_vl)
%     if isempty(IBB_target==sub)
%         continue;
%     end
    dateTemp=dateList_vl(IBB_target_vl==sub);
%     if length(dateTemp)==3
%         dateTemp=dateTemp(1:2);
%     end
%     if length(dateTemp)==4
%         error('Too many dates!! Invalid data?')
%     end
    idd=find(contains(subjectID,targetID_vl{sub}));
    dateEEGtemp=dateEEG(idd);
    for iTemp=1:length(idd)
        if ~isempty(find(dateTemp-dateEEGtemp(iTemp)<timeThres&dateTemp-dateEEGtemp(iTemp)>-timeThres))
%             idxDateTemp=find(dateTemp==dateEEGtemp(iTemp));
            idxDateTemp=find(dateTemp-dateEEGtemp(iTemp)<timeThres&dateTemp-dateEEGtemp(iTemp)>-timeThres);
            temp=find(IBB_target_vl==sub);
            target_vl=[target_vl;vl(temp(idxDateTemp(end)),:)];
%             target_diag=[target_diag;Diag(temp(idxDateTemp))];
            targetID_inUse_vl=[targetID_inUse_vl;subjectID(idd(iTemp))];
        end
    end
end
% subName=cellfun(@(x) x(1:end-9),targetID_inUse_vl,'uniformoutput',false);

% Get data points in common
targetID_common = intersect(intersect(targetID_inUse,targetID_inUse_aim),intersect(targetID_inUse_srs,targetID_inUse_vl));
idx_ados = find(contains(targetID_inUse,targetID_common));
id_ados = targetID_inUse(contains(targetID_inUse,targetID_common));
target_ados = [target_sa,target_rrb];
[~,idx] = ismember(targetID_common,id_ados);
id_ados = id_ados(idx);
EEG = EEG(:,idx_ados);
EEG = EEG(:,idx);
target_session = target_session(idx_ados);
target_session = target_session(idx);
target = target(idx_ados);
target = target(idx);
target_age = target_age(idx_ados);
target_age = target_age(idx);
target_sex = target_sex(idx_ados);
target_sex = target_sex(idx);
target_ados = target_ados(idx_ados,:);
target_ados = target_ados(idx,:);
idx_aim = find(contains(targetID_inUse_aim,targetID_common));
id_aim = targetID_inUse_aim(contains(targetID_inUse_aim,targetID_common));
[~,idx] = ismember(targetID_common,id_aim);
target_aim = target_aim(idx_aim,:);
target_aim = target_aim(idx,:);
idx_srs = find(contains(targetID_inUse_srs,targetID_common));
id_srs = targetID_inUse_srs(contains(targetID_inUse_srs,targetID_common));
[~,idx] = ismember(targetID_common,id_srs);
target_srs = target_srs(idx_srs,:);
target_srs = target_srs(idx,:);
idx_vl = find(contains(targetID_inUse_vl,targetID_common));
id_vl = targetID_inUse_vl(contains(targetID_inUse_vl,targetID_common));
[~,idx] = ismember(targetID_common,id_vl);
target_vl = target_vl(idx_vl,:);
target_vl = target_vl(idx,:);

% allMeasure = [target_sa,target_rrb,target_aim];
allMeasure= [target_ados,target_aim,target_srs,target_vl];

