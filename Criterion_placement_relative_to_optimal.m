%% Criterion_placement_relative_to_optimal
% This code computes criteria at the onset of the conditions and in the SS 
% relative to that criteria that would maximize reinforcement after having 
% substracted any (sensory or response) baseline biases. Estimates of
% criteria are obtained using a mulitple linear regression as 
% in Stüttgen et al., 2011.

%                                               Luis de la Cuesta Ferrer
%                                               25/09/2024
%% Get data and calculate criteria in the steady state (SS)

clear all
Dataframe =2; % 0 for Rats, 2 for Pigeons

All_InferredC    = [];
All_Reg_Sess_ID  = [];
All_InfCxC_lastN = [];
All_InfCxC_SS_means = [];
All_InfCxC_Base_means = [];

if Dataframe == 0
    SessionsToConsider = 3; 
    SessionsToConsiderforBT = 3; 
    load ('Long_struct')         % Rats datasets are also called Long
    Dataset = Long_schedule;
elseif Dataframe ==2
    SessionsToConsider = 3; 
    SessionsToConsiderforBT = 5; 
    load ('Pigeon_struct')
    Dataset = Pigeon_schedule;
end
StimulusMeans_OCPS = zeros(length(Dataset),5);

for iDataset =1:length(Dataset) % all animals

    allBDataLAK  = Dataset(iDataset).allBDataKDB_LAK;
    sequence     = Dataset(iDataset).sequence (1:max(allBDataLAK(:,6)));
    ChangeOfCond = Dataset(iDataset).ChangeOfCond;
    response  = allBDataLAK(:,3);
    responded = not(isnan(response)|response==0);           % Get rid of NaNs before feeding OCPS model
    UpallBDataLAK = allBDataLAK(responded,:);
    [RecZSCor,NewZSCor, MLRCoef]   = OneCritxSessionModel(allBDataLAK);
    StimulusMeans_OCPS(iDataset,:) = MLRCoef(1:5);
    InferredC_OCPS = MLRCoef(max(UpallBDataLAK(:,1)+1): end)'; 
    InferredC_OCPS = InferredC_OCPS*-1;                     % Change the criteria sign for convinience
    [All_InfCxC_SS_means(iDataset,:),InfCxC_lastN_R, InfCxC_lastN] = SteadyStateCritCondition(Dataframe,InferredC_OCPS,UpallBDataLAK,ChangeOfCond,SessionsToConsider,SessionsToConsiderforBT, sequence);
    [AllInfCxC_Baseline_mean(iDataset,:)]                          = BaselineCritCondition(sequence,InfCxC_lastN);
    
end
    All_InfCxC_SS_means_BT_centered   = All_InfCxC_SS_means - All_InfCxC_SS_means(:,1);
    %All_InfCxC_SS_means_BT_centered   = All_InfCxC_SS_means;
    All_InfCxC_Base_means_BT_centered = AllInfCxC_Baseline_mean;
    All_InfCxC_Base_means_BT_centered = [AllInfCxC_Baseline_mean(:,1),AllInfCxC_Baseline_mean(:,2:end) - All_InfCxC_SS_means(:,1)];
    y=(0.02:0.02:1);
    Nr_cond = 1:length(All_InfCxC_SS_means_BT_centered);
  
% Compute optimal criterion using fitted stimulus Means from OCPS
% Steady state criteria and Stimulus Means are Baseline corrected using the mean
% bias of each subject animal in the baseline testing sessions.

[Pred_c_optimal,c_s,perc_obtained_reinf_SS,perc_obtained_reinf_Baseline,Pred_sequence] = compute_optimal_crit (Dataframe, Dataset, StimulusMeans_OCPS,All_InfCxC_SS_means);  

order = 2;
Updated_sequence = RewriteConditionIdentity(Pred_sequence,order);
C_baselines = zeros (2,7);
[a,b] = sort(Updated_sequence,'ascend');
Pred_c_optimal=Pred_c_optimal(:,b);
offset_SS_from_optimal   = abs(All_InfCxC_SS_means_BT_centered-Pred_c_optimal);
offset_Base_from_optimal = abs(All_InfCxC_Base_means_BT_centered-Pred_c_optimal);

%% Plot (Figure 3b (Rats) and Figure 6c (Pigeons))
ConditionEnd            = [1,3,5];
Distance_Begginig_to_SS = 3;
Distance_Conditions     = 10;
orange                  = [0.8500 0.3250 0.0980];
figure

for jCond = 1:size(offset_Base_from_optimal,2)
xaxis(jCond)         = jCond*Distance_Conditions-1;
Beggining_DataPoints = offset_Base_from_optimal(:,jCond); 
Beggining_Mean_Datapoints(1,jCond) = mean(Beggining_DataPoints);
scatter(xaxis(jCond),Beggining_DataPoints,20,'blue',"x")
hold all
plot(xaxis(jCond),mean(Beggining_DataPoints),'ob','MarkerSize',7,'MarkerFace','blue')
SS_DataPoints               = offset_SS_from_optimal(:,jCond);
SS_Mean_Datapoints(1,jCond) = mean(SS_DataPoints);
scatter((xaxis(jCond)+Distance_Begginig_to_SS),SS_DataPoints,20,orange,"x")
plot(xaxis(jCond)+Distance_Begginig_to_SS,mean(SS_DataPoints),'o','MarkerFace',orange,'MarkerSize',7,'MarkerFace',orange)
end

hold on

for jCond = 1:size(offset_Base_from_optimal,2)
    xaxis1 = xaxis(jCond);
    xaxis2 = xaxis(jCond)+Distance_Begginig_to_SS;
    yaxis1 = Beggining_Mean_Datapoints(jCond);
    yaxis2 = SS_Mean_Datapoints(jCond);
    line ([xaxis1,xaxis2],[yaxis1,yaxis2]);

    if ismember (jCond,ConditionEnd)
        line([xaxis(jCond)+7, xaxis(jCond)+7],[0 2],'LineStyle','--')
    end
    hold all
end
axis ([0 max(xaxis)+7 0 2])
xticks(xaxis+1)
xticklabels({'BT','RL','RR','LL', 'LR','CL','CR'})
axis square
if Dataframe ==0
axis([0,80,0,1.2])
else
axis([0,80,0,1.5])
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
%% function [All_InfCxC_means,InfCxC_lastN_R] = SteadyStateCritCondition(Dataframe,InferredC_OCPS,UpallBDataKDB,ChangeOfCond,SessionsToConsider,SessionsToConsiderforBT,sequence)
function [All_InfCxC_means,InfCxC_lastN_R,InfCxC_lastN] = SteadyStateCritCondition(Dataframe,InferredC_OCPS,UpallBDataKDB,ChangeOfCond,SessionsToConsider,SessionsToConsiderforBT,sequence)
InfCxC = zeros(1,length(ChangeOfCond)); 
All_InfCxC_lastN = [];
All_InfCxC_means = [];

for iSess=1:max(UpallBDataKDB(:,7))                                             % Loop over sessions
    PR2(iSess)  =  length(find(UpallBDataKDB(:,3)==2 & UpallBDataKDB(:,7)==iSess)) / length(find( UpallBDataKDB(:,7)==iSess));
end

    for iCond=1:length(ChangeOfCond)                                            
        PR2xC_lastN(iCond) = mean(PR2(ChangeOfCond(iCond)-(SessionsToConsider-1):ChangeOfCond(iCond)));
        if Dataframe ==2                                                   % this if loop is needed because there is a disagreement
                                                                           % on what the saved variable ChangeOfCond means in Pigeons 
                                                                           % relative to Rats Long & Short. In Rats ChangeOfCond 
                                                                           % refers to the last session of the last conditions. In pigeons
                                                                           % it refers to the first session to the new condition 
                                                                           % (with exception of the last condition because there is not
                                                                           % any new condition following that. In that instance it refers
                                                                           % to the last session of the dataset.
            if iCond == 1
                SessionsToConsider=SessionsToConsiderforBT;
            end
            
            if iCond ~= length(ChangeOfCond)
                InfCxC_lastN(1:SessionsToConsider,iCond) = InferredC_OCPS(ChangeOfCond(iCond)-1-(SessionsToConsider-1):ChangeOfCond(iCond)-1);
            else
                InfCxC_lastN(1:SessionsToConsider,iCond) = InferredC_OCPS(ChangeOfCond(iCond)-(SessionsToConsider-1):ChangeOfCond(iCond));
            end

        else
            InfCxC_lastN(:,iCond) = InferredC_OCPS(ChangeOfCond(iCond)-(SessionsToConsider-1):ChangeOfCond(iCond));
        end
    end
N = 1:length(InfCxC);
sequence = sequence(N);
order = 2;
Updated_sequence = RewriteConditionIdentity(sequence,order);
[a,b] = sort(Updated_sequence,'ascend');
InfCxC_lastN_R = InfCxC_lastN(:,b); %Streamlined order
InfCxC_lastN_R = InfCxC_lastN_R (:,1:7); %Exclude eight session if there is one
    %Streamlined order
    BT = [N(strcmp(sequence,'BT'))];
    Other = setdiff(N,BT);
    for iCond = 1:length( InfCxC_lastN_R)
    All_InfCxC_mean(1,iCond)  = mean( nonzeros( InfCxC_lastN_R(:,iCond) ) ,1);
    end
    All_InfCxC_lastN = vertcat(All_InfCxC_lastN,InfCxC_lastN_R);
    All_InfCxC_means = vertcat(All_InfCxC_means,All_InfCxC_mean);
end
%% function [AllInfCxC_Baseline_mean] = BaselineCritCondition(sequence,InfCxC_lastN)
function [AllInfCxC_Baseline_mean] = BaselineCritCondition(sequence,InfCxC_lastN)

N = 1:12;  % create a vector with the number of initial conditions (including baselines in between real conditions)
sequence = sequence(N);
order = 2;
Updated_sequence = RewriteConditionIdentity(sequence,order);
C_baselines = zeros (2,7);
[a,b] = sort(Updated_sequence,'ascend');
for iCond = 1:7 % cut the baselines in between real conditions
    if iCond ~=1
        C_baselines(:,iCond) = InfCxC_lastN(2:3,find(Updated_sequence==iCond)-1); % select only two last sessions of the last condition
    end
end
    AllInfCxC_Baseline_mean = mean(  C_baselines  ,1);
end
%% function [Updated_sequence] = RewriteConditionIdentity(sequence,order)
function [Updated_sequence] = RewriteConditionIdentity(sequence,order)

%   so that all the conditions have the same identity regardless of the 
%   order in which order they appeared, which is how they are saved by the
%   variable sequence.
%
%   Identities:
%   (First) BT = 1, CL = 2, CL = 3, RPL = 4, RPR = 5, SOL = 6, SOR = 7
%   If there are in between BT, then those BT = 8
                                                           
%                                                       Luis Jan 2022

if order == 1
Updated_sequence = (1:1:length(sequence));
Updated_sequence(strcmp(sequence,'CL'))=2;
Updated_sequence(strcmp(sequence,'CR'))=3;
Updated_sequence(strcmp(sequence,'RPL'))=4; % RPL is called Lean L in the paper
Updated_sequence(strcmp(sequence,'RPR'))=5; % RPR           Lean R
Updated_sequence(strcmp(sequence,'SOL'))=6; % SOL           Rich L
Updated_sequence(strcmp(sequence,'SOR'))=7; % SOR           Rich R
Updated_sequence(strcmp(sequence,'BT'))=8;
Updated_sequence(1) = 1;
elseif order == 2
Updated_sequence = (1:1:length(sequence));
Updated_sequence(strcmp(sequence,'SOL'))=2;
Updated_sequence(strcmp(sequence,'SOR'))=3;
Updated_sequence(strcmp(sequence,'RPL'))=4;
Updated_sequence(strcmp(sequence,'RPR'))=5;
Updated_sequence(strcmp(sequence,'CL'))=6;
Updated_sequence(strcmp(sequence,'CR'))=7;
Updated_sequence(strcmp(sequence,'BT'))=8;
Updated_sequence(1) = 1;
end

end
%% function [Pred_c_optimal,c_s,perc_obtained_reinf_SS,perc_obtained_reinf_Baseline,Pred_sequence] = compute_optimal_crit (Dataframe, Dataset, StimulusMeans_OCPS,All_InfCxC_SS_means)
function [Pred_c_optimal,c_s,perc_obtained_reinf_SS,perc_obtained_reinf_Baseline,Pred_sequence] = compute_optimal_crit (Dataframe, Dataset, StimulusMeans_OCPS,All_InfCxC_SS_means)

% compute optimal criterion using fitted Stimulus Means from OCPS

c_optimal_all_cond = [];
c_rew_all_cond   = [];
c_unrew_all_cond = [];
c_all_all_cond   = [];
Pred_sequence = {'BT','CL','CR','RPL','RPR','SOL','SOR'};
BToffset = All_InfCxC_SS_means(:,1);
for iDataset = 1:size(StimulusMeans_OCPS,1)
    BTcentered_StimulusMeans(iDataset,:) = StimulusMeans_OCPS(iDataset,:)-BToffset(iDataset,1); % we substract the bias (offset at baseline) from all the Stim Means.
end
StimulusMeans  = BTcentered_StimulusMeans;
Pred_c_optimal = zeros(4,7);
conditions     = [1,2,3,4,5,6,7];
c_s = All_InfCxC_SS_means-All_InfCxC_SS_means(:,1);

for  iDataset =1:length(Dataset) % all animals

    AltStimSet   = StimulusMeans(iDataset,:);

    if Dataframe ==2
    d_cond1 =  [AltStimSet(1),AltStimSet(3),AltStimSet(3),AltStimSet(5)]; %BT
    else
    d_cond1 =  [AltStimSet(1),AltStimSet(2),AltStimSet(4),AltStimSet(5)]; %BT
    end
    
    d_cond2 =  [AltStimSet(1),AltStimSet(3),AltStimSet(4),AltStimSet(5)]; %CL
    d_cond3 =  [AltStimSet(1),AltStimSet(2),AltStimSet(3),AltStimSet(5)]; %CR
    d_cond4 =  [AltStimSet(1),AltStimSet(3),AltStimSet(5)]; %RPL
    d_cond5 =  [AltStimSet(1),AltStimSet(3),AltStimSet(5)]; %RPR
    d_cond6 =  [AltStimSet(1),AltStimSet(3),AltStimSet(5)]; %SOL
    d_cond7 =  [AltStimSet(1),AltStimSet(3),AltStimSet(5)]; %SOR

    % main code
    for j = 1:length(conditions)

        switch conditions(j)
            case 1 %BT
                % mu = [-1 0 0 1] * d;
                mu = d_cond1;
                p = 1/4*[1 1 1 1];
                class = [1 1 2 2];
                r = [1.0 1.0 1.0 1.0];

            case 2 %CL
                % mu = [-1 0 1/3 1] * d;
                mu = d_cond2;
                p = 1/4*[1 1 1 1];
                class = [1 2 1 2];
                r = [1.0 1.0 1.0 1.0];

            case 3 %CR
                % mu = [-1 0 1/3 1] * d;
                mu = d_cond3;
                p = 1/4*[1 1 1 1];
                class = [1 2 1 2];
                r = [1.0 1.0 1.0 1.0];

            case 4 %RPL
                % mu = [-1 0 0 1] * d;
                mu = d_cond4;
                p = 1/4*[1 2 1];
                class = [1 2 2];
                r = [1.0 0.25 0.5];

            case 5 %RPR
                % mu = [-1 0 0 1] * d;
                mu = d_cond5;
                p = 1/4*[1 2 1];
                class = [1 1 2];
                r = [0.5 0.25 1];

            case 6 %SOL
                % mu = [-1 0 1] * d;
                mu = d_cond6;
                p = 1/4*[2 1 1];
                class = [1 2 2];
                r = [1.0 1.0 1.0];

            case 7 %SOR
                % mu = [-1 0 1] * d;
                mu = d_cond7;
                p = 1/4*[1 1 2];
                class = [1 1 2];
                r = [1.0 1.0 1.0];
        end

        c = linspace(-2.5+min(mu),max(mu)+2.5,1000);
        [Rf1,Rf2] = reinforcement(c,mu,p,class,r);    % Rf1 and Rf2 are expected reinforcers per trial as a function of criterion values, FOR EACH CATEGORY
        c_optimal = fminbnd(@(c) f_optimal(c,mu,p,class,r),min(c),max(c));
        c_optimal_all_cond = [c_optimal_all_cond,c_optimal];
        Rf_s(iDataset,j) = expected_reward(c_s(j),mu,p,class,r);
        Rf_0(iDataset,j) = expected_reward(0,mu,p,class,r);
        Rf_optimal(iDataset,j) = expected_reward(c_optimal,mu,p,class,r);

    end
    Pred_c_optimal(iDataset,:)   = c_optimal_all_cond;c_optimal_all_cond = [];

end
%performance = Rf_s/Rf_optimal; % how many reinforcers can the subject expect compared to optimal behaviour?
 perc_obtained_reinf_SS =       (mean(Rf_s./Rf_optimal)); %     close all
 perc_obtained_reinf_Baseline = (mean(Rf_0./Rf_optimal)); %
end

%% function y = f_optimal(c,mu,p,class,r)
function y = f_optimal(c,mu,p,class,r)
[Rf1,Rf2] = reinforcement(c,mu,p,class,r);
y = -Rf1-Rf2;
end
 %% function [Rf1,Rf2] = reinforcement(c,mu,p,class,r)
function [Rf1,Rf2] = reinforcement(c,mu,p,class,r)

Rf1 = 0;
for i=find(class==1)
    Rf1 = Rf1 + p(i)*r(i)*normcdf(c,mu(i),1);
end

Rf2 = 0;
for i = find(class==2)
    Rf2 = Rf2 + p(i)*r(i)*(1-normcdf(c,mu(i),1));
end
end

%% function plot_densities(c,mu,p,class,r)
function plot_densities(c,mu,p,class,r)
grey = ones(1,3)*0.7;
s = 1;
y = zeros(size(c));
for i = find(class==1)
    y = y+s*p(i)*r(i)*normpdf(c,mu(i),1);
end
plot(c,y,'Color',grey) % plot sum of all class 1 distributions
y = zeros(size(c));
for i = find(class==2)
    y = y+s*p(i)*r(i)*normpdf(c,mu(i),1);
end
plot(c,y,'Color',grey)
end

%% function Rf = expected_reward(c,mu,p,class,r)
function Rf = expected_reward(c,mu,p,class,r)
[Rf1,Rf2] = reinforcement(c,mu,p,class,r);
Rf = Rf1+Rf2;
end