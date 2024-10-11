%% Caller GLM. make some choices: Dataset (1) Response or Stimuluswise (2) and erase influences after withdrawn trial or not(3)

clear all
TrialsInPast = [1,2,3,4];
Dataframe = 0;       % 0 For Rat Long schedule, 2 for Pigeons.

if Dataframe == 0
    load ('Long_struct.mat')
    Dataset = Long_schedule;
elseif Dataframe ==2
    load ('Pigeon_struct.mat')
    Dataset = Pigeon_schedule;
end

for iDataset =1:length(Dataset)
    DataMatrix       = Dataset(iDataset).allBDataKDB ;
    OneHotEncoder    = 1;  % 1 for Session by Session Regressor and 2 for Condition by Condition Regressor
    [Regressors,DependentVar,betas, pvalues, SEWeights, VIF] = fitGLM(DataMatrix,OneHotEncoder,Dataframe);
    FittedStimulusMeans(iDataset,:)       = betas(2:6);
    RewardedWeights(iDataset,:)           = betas(7:(7+length(TrialsInPast)-1) );  %After rewarded trial
    UnrewardedWeights(iDataset,:)         = betas(12:(12+length(TrialsInPast)-1)); %After unrewarded trial
end

% Compute Means of weights and Plot them
xaxis = TrialsInPast;
RewardedWeights(size(RewardedWeights,1)+1,:)     = mean(RewardedWeights,1);
UnrewardedWeights(size(UnrewardedWeights,1)+1,:) = mean(UnrewardedWeights,1);
figure;
axis square

for iDataset=1:length(Dataset)
    p1 = plot(xaxis(~isnan(RewardedWeights(iDataset,:))),RewardedWeights(iDataset,~isnan(RewardedWeights(iDataset,:))),'b--');
    hold on;
    plot(xaxis(~isnan(RewardedWeights(end,:))),RewardedWeights(end,~isnan(RewardedWeights(end,:))),'b-','LineWidth',2)
    hold on;

    p2 = plot(xaxis(~isnan(UnrewardedWeights(iDataset,:))),UnrewardedWeights(iDataset,~isnan(UnrewardedWeights(iDataset,:))),'r--');
    hold on;
    plot(xaxis(~isnan(UnrewardedWeights(end,:))),UnrewardedWeights(end,~isnan(UnrewardedWeights(end,:))),'r-','LineWidth',2)
    hold on
end

xlabel('Trials in the past')
ylabel('Weights')
title('Influence of outcome history on current choice')
legend([p1 p2],{'Reward', 'TimeOut'})
axis square


%% function[Regressors,DependentVar, betas, pvalues, SEWeights, VIF] = fitGLM(DataMatrix,OneHotEncoder, PerseverationRegressor,Dataframe)
function[Regressors,DependentVar, betas, pvalues, SEWeights, VIF] = fitGLM(DataMatrix,OneHotEncoder, Dataframe)
%
% This function fits a Generalized Linear Model to the data formatted as
% DataMatrix (including withdrawals). The dependent variable (DependentVar)
% indicates R1 (0) vs R2 (1) response. The regressors used are 1 Current
% Stimulus Evidence 2a RewardHistoryBiases (t-1,t-2,t-3 & t-4) 2b
% UnrewardedHistoryBiases (t-1,t-2,t-3 & t-4) .

% GLM fitting: We asume that the data comes from a binomial distribution and
% the stretching function (also called link function or 1/link function) is logit
% function. This is also called a logistic regression.
% The regressors are coded from a R2 point of view. Meaning that a 1 means
% a positive tendency towards a R2 response and a -1  a negative one.

% Inputs:

% (1)   DataMatrix in allBDataKDB format
% (2)   Do we include the session by session regressor? 0 =no, 1= yes. If
%       you want to see the effect of the error trials you need 1=yes. This
%       is because in our task error influneces are so small that they you
%       need to clean the stronger session by session correlations that
%       there exist in our dataset so that you can actuall see those.
% (3)   Rats (Dataframe=0) or Pigeons (Dataframe = 2)?

% Book keeping:

% Luis de la Cuesta Ferrer   11/22. Script creation


if Dataframe ==2
    ValidTrials   = find(DataMatrix(:,3)~=0);    % For Pigeons
else
    ValidTrials   = not(isnan(DataMatrix(:,3))); % For Rats
end

[I,~] = find(DataMatrix(:,2)==DataMatrix(:,3));
DataMatrix(I,5)=0;
% Which indices you take

AllIndices          = (1:1:length(DataMatrix));
% Taking all indices
IndicesofInterest   = AllIndices;

IndTrialsRewardedR1 = find( DataMatrix(IndicesofInterest,4)~=0 & DataMatrix(IndicesofInterest,3)==1 ) ;
IndTrialsPunishedR1 = find( DataMatrix(IndicesofInterest,4)==0 & DataMatrix(IndicesofInterest,3)==1 ) ;
IndTrialsRewardedR2 = find( DataMatrix(IndicesofInterest,4)~=0 & DataMatrix(IndicesofInterest,3)==2 ) ;
IndTrialsPunishedR2 = find( DataMatrix(IndicesofInterest,4)==0 & DataMatrix(IndicesofInterest,3)==2 ) ;

TrialsInPast = 1;

[LastTrialRewarded,LastTrialPunished]                 = BuildHistoryRegressors(TrialsInPast,DataMatrix, IndTrialsRewardedR1,IndTrialsRewardedR2,IndTrialsPunishedR1,IndTrialsPunishedR2);
[SecondtoLastTrialRewarded,SecondtoLastTrialPunished] = BuildHistoryRegressors(2,DataMatrix, IndTrialsRewardedR1,IndTrialsRewardedR2,IndTrialsPunishedR1,IndTrialsPunishedR2);
[ThirdtoLastTrialRewarded,ThirdtoLastTrialPunished]   = BuildHistoryRegressors(3,DataMatrix, IndTrialsRewardedR1,IndTrialsRewardedR2,IndTrialsPunishedR1,IndTrialsPunishedR2);
[FifthtoLastTrialRewarded,FifthtoLastTrialPunished]   = BuildHistoryRegressors(4,DataMatrix, IndTrialsRewardedR1,IndTrialsRewardedR2,IndTrialsPunishedR1,IndTrialsPunishedR2);

% One hot encoding (session or condition wise).
% One categorical regressor per session that absorbs the session by session
% biases (it is 1 for that session trials and 0 for the rest)

if OneHotEncoder ==1 % Session
    BlockIdentifyer = zeros( length (DataMatrix), max(DataMatrix(:,7)) );

    for iSess = 1:max(DataMatrix(:,7))
        ThisSessionTrialIndices = find (DataMatrix(:,7)==iSess);
        BlockIdentifyer(ThisSessionTrialIndices,iSess) = 1;
    end
elseif OneHotEncoder ==2 % Condition

    BlockIdentifyer = zeros( length (DataMatrix), max(DataMatrix(:,6)) );

    for iCond = 1:max(DataMatrix(:,6))
        ThisConditionTrialIndices = find (DataMatrix(:,6)==iCond);
        BlockIdentifyer(ThisConditionTrialIndices,iCond) = 1;
    end
end
% Build regressors
% 1Current Stimulus Evidence (Dummy coded)

CurrentStimulus     = zeros(max(DataMatrix(:,1)),size(DataMatrix,1))';
for iStim = 1:max(DataMatrix(:,1))
    CurrentStimulus(DataMatrix(:,1)==iStim,iStim)=1;
end

% 2Y
DependentVar = (DataMatrix(:,3)-1);

% 3 Outcome history
if OneHotEncoder
    History_Regressors = [LastTrialRewarded(:,1),SecondtoLastTrialRewarded(:,1),ThirdtoLastTrialRewarded(:,1),FifthtoLastTrialRewarded(:,1),...
        LastTrialPunished(:,1), SecondtoLastTrialPunished(:,1),ThirdtoLastTrialPunished(:,1),FifthtoLastTrialPunished(:,1),BlockIdentifyer(:,2:end)];
else
    History_Regressors = [LastTrialRewarded(:,1),SecondtoLastTrialRewarded(:,1),ThirdtoLastTrialRewarded(:,1),FifthtoLastTrialRewarded(:,1),...
        LastTrialPunished(:,1), SecondtoLastTrialPunished(:,1),ThirdtoLastTrialPunished(:,1),FifthtoLastTrialPunished(:,1)];
end

% 4a After effect module
% AfterEffect = zeros(1,size(DataMatrix,1));
% AfterEffect(2:end) = CurrentStimulus(1:end-1); % CurrentStimulus used to
% be the identity. Now it is a matrix dummy coding it.

% 5Regressors = [];

Regressors                  = [CurrentStimulus,History_Regressors];
Regressors                  = Regressors(ValidTrials,:);
DependentVar                = DependentVar(ValidTrials);


%[b,dev,stats]               = glmfit(Regressors,DependentVar,'binomial','link','logit','Constant','off');
[b,dev,stats]               = glmfit(Regressors,DependentVar,'binomial','link','logit');

betas     = stats.beta;
SEWeights = stats.se;
pvalues   = stats.p;

% Variance inflation coefficients (to check for colinearity of regressors)

%vif() computes variance inflation coefficients
%VIFs are also the diagonal elements of the inverse of the correlation matrix [1], a convenient result that eliminates the need to set up the various regressions
%[1] Belsley, D. A., E. Kuh, and R. E. Welsch. Regression Diagnostics. Hoboken, NJ: John Wiley & Sons, 1980.

R0 = corrcoef(Regressors); % correlation matrix
VIF  = diag(inv(R0))';
end

%%  function[LastTrialRewarded,LastTrialPunished] = BuildHistoryRegressors(TrialsInPast,allBDataKDB, IndTrialsRewardedR1,IndTrialsRewardedR2,IndTrialsPunishedR1,IndTrialsPunishedR2)

function[LastTrialRewarded,LastTrialPunished] = BuildHistoryRegressors(TrialsInPast,allBDataKDB, IndTrialsRewardedR1,IndTrialsRewardedR2,IndTrialsPunishedR1,IndTrialsPunishedR2)

% This function builds the regressors of the influence of a given
% outcome positive and negative for the response
% specified at trial = TrialInPast to investigate its effect on
% the current decision in the current decision. It needs as input
% IndTrialsRewardedR1/R2 and IndTrialsPunishedR1/R2 .

%                                                  Luis 21/10/2022

LastTrialRewarded = zeros(1,size(allBDataKDB,1)+TrialsInPast);
LastTrialPunished = zeros(1,size(allBDataKDB,1)+TrialsInPast);

LastTrialRewarded(1,IndTrialsRewardedR1+TrialsInPast) =-1; %Response 1 Rewarded trials            %  A note on the regressors signs: If you give the same sign to R1 correct and R1 incorrect
LastTrialRewarded(1,IndTrialsRewardedR2+TrialsInPast) =1; %Response 2 Rewarded trials             %  you would expect negative regressors for the unrewarded outcomes. This visually makes
                                                                                                  %  more sense but clashes with our negative learning rates interpretation of the error-based
LastTrialPunished(1,IndTrialsPunishedR1+TrialsInPast) =-1; %Response 1 Punished trials            %  models.
LastTrialPunished(1,IndTrialsPunishedR2+TrialsInPast) =1; %Response 2 Punished trials

LastTrialRewarded= LastTrialRewarded(:,1:length(LastTrialRewarded)-TrialsInPast)'; %reformatting of the matrix
LastTrialPunished= LastTrialPunished(:,1:length(LastTrialPunished)-TrialsInPast)'; %reformatting of the matrix

end

