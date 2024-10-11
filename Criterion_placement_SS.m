%% Criterion_placement_SS
% This code takes the raw data of either Rats (Dataframe = 0) or Pigeons
% (Dataframe = 2), recovers the session-by-session criteria through fitting
% a mulitple linear regression as in Stüttgen et al., 2011, computes steady-state
% criteria condition-wise. It provides and overview of the overall
% performance of subjects in the different experimental conditions.

%                                               Luis de la Cuesta Ferrer
%                                               25/09/2024
%% Multiple Linear regression z scores with Dummy variable 
clear all
Dataframe = 2; % 0 for Long Rats, 1 for Short Rats, 2 for pigeons.

if Dataframe ==0
load("Long_struct.mat")
Dataset  = Long_schedule;
SessionsToConsiderforBT = 3;
SessionsToConsider = 3;

elseif Dataframe ==2
load("Pigeon_struct.mat")
Dataset  = Pigeon_schedule;
Pigeons = 1;
SessionsToConsiderforBT = 5;
SessionsToConsider = 5;
end

All_InferredC    = [];
All_Reg_Sess_ID  = [];
All_InfCxC_lastN = [];
All_InfCxC_means = [];
StimulusMeans_OCPS = zeros(length(Dataset),5);

for iDataset =1:length(Dataset) % all animals 
    
allBDataKDB   = Dataset(iDataset).allBDataKDB;
ChangeOfCond  = Dataset(iDataset).ChangeOfCond;
sequence      = Dataset(iDataset).sequence;
response  = allBDataKDB(:,3);
responded = not(isnan(response)|response==0);
UpallBDataKDB = allBDataKDB(responded,:);
[RecZSCor,NewZSCor, MLRCoef]         = OneCritxSessionModel(UpallBDataKDB);
StimulusMeans_OCPS(iDataset,:) = MLRCoef(1:5);             
InfCxC = zeros(1,length(ChangeOfCond)); % Inferred C from 

for i=1:max(UpallBDataKDB(:,7))                                             % Loop over sessions
    PR2(i)  =  length(find(UpallBDataKDB(:,3)==2 & UpallBDataKDB(:,7)==i)) / length(find( UpallBDataKDB(:,7)==i));
end

InferredC_OCPS = MLRCoef(max(UpallBDataKDB(:,1)+1): end)';
InferredC_OCPS = InferredC_OCPS*-1;
[All_InfCxC_SS_means(iDataset,:),InfCxC_lastN_R, InfCxC_lastN] = SteadyStateCritCondition(Dataframe,InferredC_OCPS,UpallBDataKDB,ChangeOfCond,SessionsToConsider,SessionsToConsiderforBT, sequence);
end
    All_InfCxC_means_SS_BT_centered = All_InfCxC_SS_means - All_InfCxC_SS_means(:,1);
    GrandCxC_mean_SS_BT_centered    = mean( All_InfCxC_means_SS_BT_centered,1);

%% Plot (Figure 2e)
figure('color','w');
colors = [0, 0, 0; 128, 128, 0; 230, 0, 160 ; 255, 179, 0; 128, 0, 64; 0, 255, 128; 48, 219, 233];
colors = colors./255;
xaxis  = (1:1:7);
ax = axes('NextPlot','add','Xlim',[0.0,8],'Ylim',[-1.5,1.5]);

arrayfun( @(i) plot(ax, xaxis(i),All_InfCxC_means_SS_BT_centered(:,i),...
'o','MarkerEdgeColor',colors(i,1:3),'MarkerFaceColor',colors(i,1:3),'MarkerSize',4 ), 1:length( All_InfCxC_means_SS_BT_centered));            % Plot individual sessions       

arrayfun( @(i) plot(ax, xaxis(i),GrandCxC_mean_SS_BT_centered(i),...  
'*','color',colors(i,1:3),'MarkerSize',20 ),...
1:length( All_InfCxC_means_SS_BT_centered) );

hold on
line(ax(1),[0 8],[0 0],'Color','black','LineStyle','--')
xticks(xaxis)
xticklabels({'BT','RL','RR','LL', 'LR','CL','CR'})
ylabel(' Reconstruced c')
xlabel('Conditions')
axis square

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
                InfCxC_lastN(1:SessionsToConsiderforBT,iCond) = InferredC_OCPS(ChangeOfCond(iCond)-1-(SessionsToConsiderforBT-1):ChangeOfCond(iCond)-1);
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