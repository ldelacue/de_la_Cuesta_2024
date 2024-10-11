%% Fit criterion by using learning models
 
% We first track the session-by-session criteria by using the One-criterion-per-session model
% and then we fit these data with the criteria resulting from using the fitted parameters 
% obtained through the different SDT learning-models.

% Main possible models included in the manuscript are: 

% Integrate Rewards                                                             Model = 'model 3a'
% Integrate Reward Omissions                                                    Model = 'model 1a'
% Integrate Rewards & Rew. Omissions                                            Model = 'model 2'
% Integrate Rewards w/ Stim. Learning rates                                     Model = 'model 3b'
% Integrate Rewards w/ Stim. Learning rates (reduced)                           Model = 'model 3b (red)'
% Integrate Rewards w/ Relative Rew.Rates                                       Model = 'model 3a NUNP'
% Integrate Rewards w/ Stim. Learning rates & Relative Rew.Rates                Model = 'model 3b NUNP'
% Integrate Rewards w/ Stim. Learning rates (reduced) & Relative Rew.Rates      Model = 'model 3b (red) NUNP'


%% Reconstruct criterion in a session by session fashion (using Multiple Linear Regression)

% Choose Dataset and Model 
clear all
Model    = 'model 3a'; % see possible models above
Dataframe = 2;          % 0 For Rat Long schedule, 2 for Pigeons 
ProvideStimMeans = 0;  % You provide Stim Means from KDB (1) or you fit them de novo (0)

if Dataframe == 0
    load ('Long_struct.mat')
    load ('AllFittedParams_BICs_4Rats.mat')
    Dataset = Long_schedule;
    Subjects = 'Rat';
elseif Dataframe ==2
    load ('Pigeon_struct.mat')
    load ('AllFittedParams_BICs_4Pigeons.mat')
    load ('AllFittedParams_BICs_3b_red4Pigeons.mat') % this set of parameters carries the reduced model versions that include the SLR (with only 2LR instead of one per stimulus)
    Dataset = Pigeon_schedule;
    Subjects = 'Pigeon';
end

for iDataset = 1:length(Dataset)

allBDataLAK = Dataset(iDataset).allBDataKDB_LAK;
sequence    = Dataset(iDataset).sequence;
response    = allBDataLAK(:,3);                     % was the response 1 or 2
responded   = not(isnan(response)|response==0);
allBDataLAK = allBDataLAK(responded,:);
Datamatrix  = allBDataLAK;

% Use the KDB fitted parameters to reconstruct criterion

Actual_reward = zeros(length(allBDataLAK(:,1)),1);
Actual_no_reward = zeros(length(allBDataLAK(:,1)),1);
Actual_no_reward( allBDataLAK(:,3)==1 & allBDataLAK(:,2)==2 )  = +1;
Actual_no_reward( allBDataLAK(:,3)==1 & allBDataLAK(:,4)==0 )  = +1;
Actual_no_reward( allBDataLAK(:,3)==2 & allBDataLAK(:,2)==1 )  = -1;
Actual_no_reward( allBDataLAK(:,3)==2 & allBDataLAK(:,4)==0 )  = -1;
Actual_reward( allBDataLAK(:,3)==1 & allBDataLAK(:,2)==1 &  allBDataLAK(:,4) ==1) = -1;
Actual_reward( allBDataLAK(:,3)==2 & allBDataLAK(:,2)==2 &  allBDataLAK(:,4) ==1) = +1;
Actual_stimulus  = allBDataLAK(:,1);
Actual_condition = allBDataLAK(:,6);

try
KDB_Fitted_StimulusMeans = w(1:5);
catch
end

if  strcmp(Model,'model 1a')
Learning_rates = Allw_1a(iDataset,6);
gamma = Allparam_1a(iDataset);
KDB_Fitted_StimulusMeans = Allw_1a(iDataset,1:5);

elseif strcmp(Model,'model 2')
Learning_rates = Allw_2(iDataset,6:end);
gamma = Allparam_2(iDataset);
KDB_Fitted_StimulusMeans = Allw_2(1:5);

elseif strcmp(Model,'model 3a')
Learning_rates = Allw_3a(iDataset,6);
gamma = Allparam_3a(iDataset);
KDB_Fitted_StimulusMeans = Allw_3a(iDataset,1:5);

elseif strcmp(Model,'model 3a NUNP')
Learning_rates = Allw_3a_NUNP(iDataset,6);
gamma = Allparam_3a_NUNP(iDataset);
KDB_Fitted_StimulusMeans = Allw_3a_NUNP(iDataset,1:5);

elseif strcmpi(Model,'model 3b')
Learning_rates = Allw_3b(iDataset,6:end);
Learning_rates(Learning_rates<0) =0; % Correct learning rates that are below 0 to 0
gamma = Allparam_3b(iDataset);
KDB_Fitted_StimulusMeans = Allw_3b(iDataset,1:5);

elseif strcmpi(Model,'model 3b (red)')
Learning_rates = [Allw_3b_red(iDataset,6), Allw_3b_red(iDataset,7),Allw_3b_red(iDataset,7),...
                  Allw_3b_red(iDataset,7), Allw_3b_red(iDataset,6)];
Learning_rates(Learning_rates<0) =0; % Correct learning rates that are below 0 to 0
gamma = Allparam_3b_red(iDataset);
KDB_Fitted_StimulusMeans = Allw_3b_red(iDataset,1:5);

elseif strcmp(Model,'model 3b NUNP')
Learning_rates = Allw_3b_NUNP(iDataset,6:end);
Learning_rates(Learning_rates<0) =0; % Correct learning rates that are below 0 to 0
gamma = Allparam_3b_NUNP(iDataset);
KDB_Fitted_StimulusMeans = Allw_3b_NUNP(iDataset,1:5);

elseif strcmp(Model,'model 3b (red) NUNP')
Learning_rates = [Allw_3b_red_NUNP(iDataset,6), Allw_3b_red_NUNP(iDataset,7),Allw_3b_red_NUNP(iDataset,7),...
                  Allw_3b_red_NUNP(iDataset,7), Allw_3b_red_NUNP(iDataset,6)];
Learning_rates(Learning_rates<0) =0; % Correct learning rates that are below 0 to 0
gamma = Allparam_3b_red_NUNP(iDataset);
KDB_Fitted_StimulusMeans = Allw_3b_red_NUNP(iDataset,1:5);
end


%criterion trial 1
c(1) = mean(KDB_Fitted_StimulusMeans); % start the criterion position with the mean bias obtained from the mean (stimulus Means)

%criterion rest of trials
for iTrial = 2:length(allBDataLAK(:,1))

    if strcmp(Model,'model 1a')      % punishment learning 1 learning rates
    c(iTrial) = c(iTrial-1)*gamma + Learning_rates *Actual_no_reward(iTrial);
    
   
    elseif strcmp(Model,'model 2')  % income + punishment learning 10 learning rates
        
        if allBDataLAK(iTrial,5)==1  % if this is an unrewarded trials
                c(iTrial) = c(iTrial-1)*gamma + Learning_rates (2)*Actual_no_reward(iTrial);   %error learning rates take values from 6-10
        else
                c(iTrial) = c(iTrial-1)*gamma + Learning_rates (1)*Actual_reward(iTrial);      %income learning rates take values from 1-5
        end
      
    elseif strcmp(Model,'model 3a')  % classic income based
    c(iTrial) = c(iTrial-1)*gamma + Learning_rates *Actual_reward(iTrial);
    
    elseif strcmp(Model,'model 3a NUNP')  % income based + VL + LL
        if allBDataLAK(iTrial,5)==1  % if this is an unrewarded trials
            c(iTrial) = c(iTrial-1);
        else
            c(iTrial) = c(iTrial-1)*gamma + Learning_rates *Actual_reward(iTrial);
        end
   
    elseif strcmp(Model,'model 3b') || strcmp(Model,'model 3b (red)')  % income based + VL + LL
    c(iTrial) = c(iTrial-1)*gamma + Learning_rates (Actual_stimulus(iTrial))*Actual_reward(iTrial);

    elseif strcmp(Model,'model 3b NUNP') || strcmp(Model,'model 3b (red) NUNP')  % income based + VL + LL
        if allBDataLAK(iTrial,5)==1  % if this is an unrewarded trials
            c(iTrial) = c(iTrial-1);
        else
            c(iTrial) = c(iTrial-1)*gamma + Learning_rates (Actual_stimulus(iTrial))*Actual_reward(iTrial);
        end
    end
end

c = c';
for iSess = 1:max(Datamatrix(:,7))
c_sessionKDB(iSess) = mean(c (Datamatrix(:,7)==iSess)) ;
end

%______________________________________________________
% Fit the one-criterion-per-session model

if ProvideStimMeans                         % You provide Stim Means from KDB or you fit them
[RecZSCor,NewZSCor, MLRCoef] = OneCritxSessionModel(Datamatrix,KDB_Fitted_StimulusMeans);
else
[RecZSCor,NewZSCor, MLRCoef] = OneCritxSessionModel(Datamatrix);
StimulusMeans_OCPS = MLRCoef(1:5); 
end

if ~ProvideStimMeans  % if you fit the means dont take the first 5 values because these are the stimulus means
    InferredC_OCPS = MLRCoef(6: end)+mean(StimulusMeans_OCPS) ;
else                  % if you provided them, then take all the means
    InferredC_OCPS = MLRCoef ;
end

InferredC_OCPS = InferredC_OCPS';
%__________________________________________________________________________
% Plotting

% First set the plot with conditions and divisory lines

figure
ConditionLabels = sequence;
block = allBDataLAK(:,7);                 % which block was the trial in

response = allBDataLAK(:,3);              % was the response 1 or 2
responded = not(isnan(response)|response==0);
condition = allBDataLAK(:,6);             % which condition (integer)
condition = condition(responded);
block     = block(responded);
xaxis = 'block';
conditionset = unique(condition);
n = length(conditionset);
i = 1;
conditionEdges = [];
xlabels = {};
lastcondition = condition(1);

for t = 1:length(condition)
    if not(condition(t)==lastcondition) % change detected
        if strcmpi(xaxis,'block')
            xlabels{i} = num2str(block(t)-1); 
        else
            xlabels{i} = num2str(trialnumber(t)-1); 
        end
        conditionEdges(i) = block(t)-1; 
        hold on
        plot([conditionEdges(i),conditionEdges(i)],[-1.75,1.75],'k:')
        hold on
        text(conditionEdges(i)+1,1.25,ConditionLabels{condition(t)})
        i=i+1;
    end
    lastcondition = condition(t);
end

if strcmpi(xaxis,'block')
    xlabels{i}=num2str(block(t));
else
    xlabels{i}=num2str(numtrials);
end

conditionEdges(i)=block(t);
plot([conditionEdges(i),conditionEdges(i)],[-1.75,1.5],'k:')
set(gca,'xtick',sort(conditionEdges));
set(gca,'xticklabel',xlabels);
text(1,1.25, ConditionLabels{condition(1)})

title(convertCharsToStrings(Subjects)',num2str(iDataset)')
Mean_Offset= mean(c_sessionKDB)-mean(InferredC_OCPS);
hold on;p1=plot((1:length(InferredC_OCPS))-0.5,c_sessionKDB(1:length(InferredC_OCPS))-Mean_Offset,'LineWidth',1.5,'color',[255/255 102/255 0/255]);
hold on;p2=plot((1:length(InferredC_OCPS))-0.5,InferredC_OCPS,'-k','LineWidth',1.5);
StimulusMeans_Offset = mean(KDB_Fitted_StimulusMeans)-mean(StimulusMeans_OCPS);

legend([p1,p2],{convertCharsToStrings(Model)','1 crit x Session'})
axis([0, iSess+3, -1.5,+1.5])

end


%% function [RecZSCor,NewZSCor,MLRCoef] = OneCritxSessionModel(UpallBDataKDB,KDB_Fitted_StimulusMeans)
function [RecZSCor,NewZSCor,MLRCoef] = OneCritxSessionModel(UpallBDataKDB,KDB_Fitted_StimulusMeans)

% One Criterion per Session Model   This function uses UpallBDataKDB (already without withdrawals)
%                                   and returns the session Z scored response probabilities (as calculated from
%                                   integrating the area under the curve up to the criterion
%                                   (NewZSCor) and the reconstructed Z scores (using the
%                                   coefficients from a multiple regression
%                                   model (RecZSCor).The regression model encodes the
%                                   stimuli and the session identity in a
%                                   dummy fashion.

%                                   Choose if you want the model to fit the Stimulus Means itself, or you
%                                   want to provide them (in my case from KDB fit, the SDT learning models).

if ~exist('KDB_Fitted_StimulusMeans','var') || isempty(KDB_Fitted_StimulusMeans); KDB_Fitted_StimulusMeans = false; end


CumStimPresentionEachSession = [];
for i=1:max(UpallBDataKDB(:,7))                                             % Loop over sessions

    TrialsxSess(i)=length( find( UpallBDataKDB(:,7)==i));

    for j=1:max(UpallBDataKDB(:,1))                                         % Loop over stimuli. Find how many leftward choices..
                                                                            % there were in each stimulus in each session

        LeftNRxSxS(i,j)= length( find( UpallBDataKDB(:,7)==i & ...          % Number of left responses/Stimulus/session
            UpallBDataKDB(:,1)==j & UpallBDataKDB(:,3)==2));

        SxS(i,j)= length( find( UpallBDataKDB(:,7)==i & ...
            UpallBDataKDB(:,1)==j ));                                       % Stimulus per session
       
        
    end
    StimThisSess = unique(UpallBDataKDB(UpallBDataKDB(:,7)==i,1));
    CumStimPresentionEachSession = [CumStimPresentionEachSession;StimThisSess];
    LeftRxS(i,:) = LeftNRxSxS(i,:)./SxS(i,:);                               % Convert nr. of leftward to relative frequency
end

% Correction for 0 per cent and 100 per cent responses to L

Corr1=find(LeftRxS==1); % 100 per cent
Corr2=find(LeftRxS==0); % 0 per cent

for i=1:length(Corr1)
    LeftRxS(Corr1(i))= 1-(1/(2*SxS(Corr1(i))));
end

for i=1:length(Corr2)
    LeftRxS(Corr2(i))= 1/(2*SxS(Corr2(i)));
end

ZScored_LeftwR=norminv(LeftRxS);
ZScored_LeftwR_vector = [];

for iSess=1:length(ZScored_LeftwR)
    ZScored_LeftwR_vector = [ZScored_LeftwR_vector;ZScored_LeftwR(iSess,:)'];
end

ZScored_LeftwR_vector = ZScored_LeftwR_vector(not(isnan(ZScored_LeftwR_vector)));                        % Get rid of NaN
SmallDummy = zeros(5,5);
DummyVar=[];
Sessions_dummy = zeros (length(SxS)*size(SxS,2),max(UpallBDataKDB(:,7)));

% Build the dummy variable

for i=1:max(UpallBDataKDB(:,7))                                         % loop over sessions

    if i==1
        Sessions_dummy (1:5*i,i)= 1;
    else
        Sessions_dummy ((i-1)*5+1:5*i,i)= 1;
    end

    for j=1:length(1:max(UpallBDataKDB(:,1)))                           % loop over stimuli
        anys = unique(UpallBDataKDB(UpallBDataKDB(:,7)==i,1));

        if anys(anys == j)
            SmallDummy(j,j)=1;
        end

    end
    
    DummyVar=[DummyVar;SmallDummy];
    SmallDummy=[];

end

FinalDummy = [DummyVar,Sessions_dummy];
FinalDummy(FinalDummy(:,1)==0 & FinalDummy(:,2)==0 & ...
FinalDummy(:,3)==0 & FinalDummy(:,4)==0 & FinalDummy(:,5)==0,:)=[];     % get rid of the rows with no stim
% FinalDummy = [FinalDummy,ones(size(FinalDummy,1),1)];                 % if you add intercept ( to capture the general bias over the sessions)

if KDB_Fitted_StimulusMeans                                             % if you provide the Stimulus Means regress its influence out from the dependent variable linearly
                                                                        % you should find it back as a weight. 

    FittedMeans_Contribution = KDB_Fitted_StimulusMeans(CumStimPresentionEachSession);
    ZScoreprep1 = ZScored_LeftwR_vector-FittedMeans_Contribution';
    FinalDummy = FinalDummy(:,6:end);
    [b,bint,r,rint,stats] = regress(ZScoreprep1,FinalDummy); 
else
    [b,bint,r,rint,stats] = regress(ZScored_LeftwR_vector,FinalDummy);
end

Regressors=b;

% Predicted Z-scores (dij = xi + cj) the z-score for stimulus i in
% session j is recontructed as the coef of session j+ coef stim i. As
% in Stüttgen 2011.

for iStim=1:length(1:max(UpallBDataKDB(:,1)))
    for jSess=1:max(UpallBDataKDB(:,7))
        if ~KDB_Fitted_StimulusMeans   
        RecZSCor(jSess,iStim) = Regressors(iStim)+Regressors(5+jSess);
        else
        RecZSCor(jSess,iStim) = KDB_Fitted_StimulusMeans(iStim)+Regressors(jSess);
        end
    end
end
MLRCoef  = b;
RecZSCor = RecZSCor(not(isnan(ZScored_LeftwR_vector)));
NewZSCor = ZScored_LeftwR(not(isnan(ZScored_LeftwR_vector)));
end
