%% Predictions vs observations of the different models based on fitted data
% This code takes the raw data of either Rats (Dataframe = 0) or Pigeons
% (Dataframe = 2), recovers the session-by-session criteria through fitting
% a mulitple linear regression as in Stüttgen et al., 2011, computes steady-state
% criteria condition-wise and plots it against the predicted criteria from
% the diffferent SDT learning-models (also called KDB models) considered in
% the manuscript to generate Figures 3a, 5e, 6d and 6h.

% The code is made to be run entirely. Section 1 needs to preceed section
% 2 and 3 as it loads the data.

%                                               Luis de la Cuesta Ferrer
%                                               25/09/2024
%% Run One-Criterion-per-Session model and retrieve average steady state criteria per condition (Get Data)

clear all
Dataframe = 0; % 0 for Long Rats, 2 for pigeons.

if Dataframe ==0
load("Long_struct.mat")
Dataset  = Long_schedule;
SessionsToConsiderforBT = 3;
SessionsToConsider = 3;
elseif Dataframe ==2
load("Pigeon_struct.mat")
Dataset  = Pigeon_schedule;
SessionsToConsiderforBT = 5;
SessionsToConsider      = 3;
Pigeons = 1;
end

StimulusMeans_OCPS = zeros(length(Dataset),5);

for iDataset  = 1:length(Dataset) % all animals 
allBDataKDB   = Dataset(iDataset).allBDataKDB;
ChangeOfCond  = Dataset(iDataset).ChangeOfCond;
sequence      = Dataset(iDataset).sequence;
response  = allBDataKDB(:,3);
responded = not(isnan(response)|response==0);           % Get rid of NaNs before feeding OCPS model
UpallBDataKDB = allBDataKDB(responded,:);
[RecZSCor,NewZSCor, MLRCoef]   = OneCritxSessionModel(UpallBDataKDB);
StimulusMeans_OCPS(iDataset,:) = MLRCoef(1:5);    
InferredC_OCPS = MLRCoef(max(UpallBDataKDB(:,1)+1): end)'; 
InferredC_OCPS = InferredC_OCPS*-1;                               % change the sign for convinience.
[All_InfCxC_means(iDataset,:)] = SteadyStateCritCondition(Dataframe,InferredC_OCPS,UpallBDataKDB,ChangeOfCond,SessionsToConsider,SessionsToConsiderforBT, sequence);
All_InfCxC_means_BT_centered = All_InfCxC_means - All_InfCxC_means(:,1);
GrandCxC_mean_BT_centered    = mean( All_InfCxC_means_BT_centered,1);
end


%% Plot 1 Predictions from 4 first models - (Figure 3a (Rats) & Figure 6d (Pigeons))

Nr_animals    = 4;
if Dataframe==0
    axes= [-1.25, 1.25, -1.25, 1.25];

elseif Dataframe ==2
    axes= [-2, 2, -2, 2];
end

[Regression_coef, P_vals,Pred_c_opt_BT_centered, Pred_c_IR_BT_centered,...
 Pred_c_IR_SLR_BT_centered, Pred_c_IR_IRO_BT_centered, Pred_c_IRO_BT_centered,...
 Pred_c_IR_SLR_RR_BT_centered,Totals_C ]=compute_correlations_all_models(Dataframe,StimulusMeans_OCPS,All_InfCxC_means_BT_centered);

sequence = {'BT', 'CL','CR','LL','LR','RL','RR'}; %experimental Conditions in plotting order
figure('color','w');
colors = [0, 0, 0; 0, 255, 128; 48, 219, 233; 255, 179, 0; 128, 0, 64; 128, 128, 0;230, 0, 160];
colors = colors./255;
ax = arrayfun( @(i) subplot(1,4,i,'NextPlot','add','Box','off'), [1:4] );

for row=1:Nr_animals

    arrayfun( @(i) plot(ax(1),Pred_c_opt_BT_centered(row,i),Totals_C(row,i),...  % Plot avg last 3 sessions
        'o','MarkerFaceColor',colors(i,1:3),'MarkerEdgeColor',colors(i,1:3),'MarkerSize',5.2  ),...
        1:7);
    line(ax(1),[-3 3],[-3 3],'Color','black','LineStyle','--')
    text(ax(1),axes(1)+0.5,axes(2)-0.0,sprintf('r^2 = %s',num2str(Regression_coef(1)^2, '% 5.2f')) )
    text(ax(1),axes(1)+0.5,axes(2)-0.25,sprintf('p = %s',num2str(P_vals(1), '% 5.3f')) )
    text(ax(1),axes(1)+0.5,axes(2)-0.5,sprintf('N = %s',num2str(Nr_animals)) )
    axis(ax(1), axes )
    xlabel(ax(1),'Predicted c Optimal')
    ylabel(ax(1),'Experimental c')
    pbaspect(ax(1),[1 1 1]);

    arrayfun( @(i) plot(ax(2),Pred_c_IR_BT_centered(row,i),Totals_C(row,i),...  % Plot avg last 3 sessions
        'o','MarkerFaceColor',colors(i,1:3),'MarkerEdgeColor',colors(i,1:3),'MarkerSize',5.2 ),...
        1:7);
    line(ax(2),[-3 3],[-3 3],'Color','black','LineStyle','--')
    text(ax(2),axes(1)+0.5,axes(2)-0.0,sprintf('r^2 = %s',num2str(Regression_coef(2)^2, '% 5.2f')) )
    text(ax(2),axes(1)+0.5,axes(2)-0.25,sprintf('p = %s',num2str(P_vals(2), '% 5.3f')) )
    text(ax(2),axes(1)+0.5,axes(2)-0.5,sprintf('N = %s',num2str(Nr_animals)) )
    axis(ax(2), axes )
    xlabel(ax(2),'Predicted c Integrate Rewards')
    pbaspect(ax(2),[1 1 1]);

    arrayfun( @(i) plot(ax(3),Pred_c_IRO_BT_centered(row,i),Totals_C(row,i),...  % Plot avg last 3 sessions
    'o','MarkerFaceColor',colors(i,1:3),'MarkerEdgeColor',colors(i,1:3),'MarkerSize',5.2 ),...
    1:7);
    line(ax(3),[-3 3],[-3 3],'Color','black','LineStyle','--')
    text(ax(3),axes(1)+0.5,axes(2)-0.0,sprintf('r^2 = %s',num2str(Regression_coef(3)^2, '% 5.2f')) )
    text(ax(3),axes(1)+0.5,axes(2)-0.25,sprintf('p = %s',num2str(P_vals(3), '% 5.2f')) )
    text(ax(3),axes(1)+0.5,axes(2)-0.5,sprintf('N = %s',num2str(Nr_animals)) )
    axis(ax(3), axes )
    xlabel(ax(3),'Predicted c Integrate Reward Omissions')
    pbaspect(ax(3),[1 1 1]);

    arrayfun( @(i) plot(ax(4),Pred_c_IR_IRO_BT_centered(row,i),Totals_C(row,i),...  % Plot avg last 3 sessions
    'o','MarkerFaceColor',colors(i,1:3),'MarkerEdgeColor',colors(i,1:3),'MarkerSize',5.2 ),...
    1:7);
    line(ax(4),[-3 3],[-3 3],'Color','black','LineStyle','--')
    text(ax(4),axes(1)+0.5,axes(2)-0.0,sprintf('r^2 = %s',num2str(Regression_coef(4)^2, '% 5.2f')) )
    text(ax(4),axes(1)+0.5,axes(2)-0.25,sprintf('p = %s',num2str(P_vals(4), '% 5.2f')) )
    text(ax(4),axes(1)+0.5,axes(2)-0.5,sprintf('N = %s',num2str(Nr_animals)) )
    axis(ax(4), axes )
    xlabel(ax(4),'Predicted c Integrate Rewards & Reward Omissions')
    pbaspect(ax(4),[1 1 1]);

end
%% Plot 2: Final Plot including the IR-modulations (Figure 5e (Rats) & Figure 6h (Pigeons))

%Dataframe = 0; % 0 for Rats (Long), 2 for Pigeons 
Nr_animals    = 4;
if Dataframe==0
    axes= [-1.25, 1.25, -1.25, 1.25];

elseif Dataframe ==2
    axes= [-2, 2, -2, 2];
end
[Regression_coef, P_vals,Pred_c_opt_BT_centered, Pred_c_IR_BT_centered,...
 Pred_c_IR_VL_BT_centered, Pred_c_IR_IRO_BT_centered, Pred_c_IRO_BT_centered,...
 Pred_c_IR_VL_NUNP_BT_centered,Totals_C ]=compute_correlations_all_models(Dataframe,StimulusMeans_OCPS,All_InfCxC_means_BT_centered);

sequence = {'BT', 'CL','CR','LL','LR','RL','RR'};
figure('color','w');
colors = [0, 0, 0; 0, 255, 128; 48, 219, 233; 255, 179, 0; 128, 0, 64; 128, 128, 0;230, 0, 160];
colors = colors./255;
ax = arrayfun( @(i) subplot(1,3,i,'NextPlot','add','Box','off'), [1:3] );

for row=1:Nr_animals

    arrayfun( @(i) plot(ax(1),Pred_c_IR_BT_centered(row,i),Totals_C(row,i),...  % Plot avg last 3 sessions
        'o','MarkerFaceColor',colors(i,1:3),'MarkerEdgeColor',colors(i,1:3),'MarkerSize',5.2 ),...
        1:7);
    line(ax(1),[-3 3],[-3 3],'Color','black','LineStyle','--')
    text(ax(1),axes(1)+0.5,axes(2)+0.25,sprintf('r = %s',num2str(Regression_coef(2), '% 5.2f')) )
    text(ax(1),axes(1)+0.5,axes(2)-0.0, sprintf('r^2 = %s',num2str(Regression_coef(2)^2, '% 5.2f')) )
    text(ax(1),axes(1)+0.5,axes(2)-0.25,sprintf('p = %s',num2str(P_vals(2), '% 5.3f')) )
    text(ax(1),axes(1)+0.5,axes(2)-0.5, sprintf('N = %s',num2str(Nr_animals)) )
    axis(ax(1), axes )
    xlabel(ax(1),'Predicted c Integrate Rewards')
    ylabel(ax(1),'Experimental c')
    pbaspect(ax(1),[1 1 1]);

    arrayfun( @(i) plot(ax(2),Pred_c_IR_SLR_BT_centered(row,i),Totals_C(row,i),...  % Plot avg last 3 sessions
    'o','MarkerFaceColor',colors(i,1:3),'MarkerEdgeColor',colors(i,1:3),'MarkerSize',5.2 ),...
    1:7);
    line(ax(2),[-3 3],[-3 3],'Color','black','LineStyle','--')
    text(ax(2),axes(1)+0.5,axes(2)+0.25,sprintf('r = %s',num2str(Regression_coef(5), '% 5.2f')) )
    text(ax(2),axes(1)+0.5,axes(2)-0.0,sprintf('r^2 = %s',num2str(Regression_coef(5)^2, '% 5.2f')) )
    text(ax(2),axes(1)+0.5,axes(2)-0.25,sprintf('p = %s',num2str(P_vals(5), '% 5.3f')) )
    text(ax(2),axes(1)+0.5,axes(2)-0.5,sprintf('N = %s',num2str(Nr_animals)) )
    axis(ax(2), axes )
    xlabel(ax(2),'Predicted c Integrate Rewards - Sensory-dependent LR')
    legend(sequence)
    legend('boxoff')
    pbaspect(ax(2),[1 1 1]);

    arrayfun( @(i) plot(ax(3),Pred_c_IR_SLR_RR_BT_centered(row,i),Totals_C(row,i),...  % Plot avg last 3 sessions
    'o','MarkerFaceColor',colors(i,1:3),'MarkerEdgeColor',colors(i,1:3),'MarkerSize',5.2 ),...
    1:7);
    line(ax(3),[-3 3],[-3 3],'Color','black','LineStyle','--')
    text(ax(3),axes(1)+0.5,axes(2)+0.25,sprintf('r = %s',num2str(Regression_coef(6), '% 5.2f')) )
    text(ax(3),axes(1)+0.5,axes(2)-0.0,sprintf('r^2 = %s',num2str(Regression_coef(6)^2, '% 5.2f')) )
    text(ax(3),axes(1)+0.5,axes(2)-0.25,sprintf('p = %s',num2str(P_vals(6), '% 5.3f')) )
    text(ax(3),axes(1)+0.5,axes(2)-0.5,sprintf('N = %s',num2str(Nr_animals)) )
    axis(ax(3), axes )
    xlabel(ax(3),'Predicted c Integrate Rewards - Sensory-dependent LR - Relative Rew. Diffs')
    legend(sequence)
    legend('boxoff')
    pbaspect(ax(3),[1 1 1]);
end

%% Plot3: Comparison optimal predictions with IR model versions (not ploted in paper, only r is mentioned)

%Dataframe = 0; % 0 for Rats (Long), 2 for Pigeons 

[Regression_coef, P_vals,Pred_c_opt_BT_centered, Pred_c_IR_BT_centered,...
 Pred_c_IR_VL_BT_centered, Pred_c_IR_IRO_BT_centered, Pred_c_IRO_BT_centered,...
 Pred_c_IR_VL_NUNP_BT_centered,Totals_C ]=compute_correlations_all_models(Dataframe,StimulusMeans_OCPS,All_InfCxC_means_BT_centered);

sequence = {'BT', 'CL','CR','LL','LR','RL','RR'};
figure('color','w');
colors = [0, 0, 0; 0, 255, 128; 48, 219, 233; 255, 179, 0; 128, 0, 64; 128, 128, 0;230, 0, 160];
colors = colors./255;
ax = arrayfun( @(i) subplot(1,3,i,'NextPlot','add','Box','off'), [1:3] );

for row=1:Nr_animals

    arrayfun( @(i) plot(ax(1),Pred_c_IR_BT_centered(row,i),Pred_c_opt_BT_centered(row,i),...  % Plot avg last 3 sessions
    'o','MarkerFaceColor',colors(i,1:3),'MarkerEdgeColor',colors(i,1:3),'MarkerSize',5.2 ),...
    1:7);
    line(ax(1),[-3 3],[-3 3],'Color','black','LineStyle','--')
    text(ax(1),axes(1)+0.5,axes(2)+0.25,sprintf('r = %s',num2str(Regression_coef(7), '% 5.2f')) )
    text(ax(1),axes(1)+0.5,axes(2)-0.0, sprintf('r^2 = %s',num2str(Regression_coef(7)^2, '% 5.2f')) )
    text(ax(1),axes(1)+0.5,axes(2)-0.25,sprintf('p = %s',num2str(P_vals(7), '% 5.3f')) )
    text(ax(1),axes(1)+0.5,axes(2)-0.5, sprintf('N = %s',num2str(Nr_animals)) )
    axis(ax(1), axes )
    xlabel(ax(1),'Predicted c Integrate Rewards')
    ylabel(ax(1),'Predicted c Optimal')
    pbaspect(ax(1),[1 1 1]);

    arrayfun( @(i) plot(ax(2),Pred_c_IR_SLR_BT_centered(row,i),Pred_c_opt_BT_centered(row,i),...  % Plot avg last 3 sessions
    'o','MarkerFaceColor',colors(i,1:3),'MarkerEdgeColor',colors(i,1:3),'MarkerSize',5.2 ),...
    1:7);
    line(ax(2),[-3 3],[-3 3],'Color','black','LineStyle','--')
    text(ax(2),axes(1)+0.5,axes(2)+0.25,sprintf('r = %s',num2str(Regression_coef(8), '% 5.2f')) )
    text(ax(2),axes(1)+0.5,axes(2)-0.0,sprintf('r^2 = %s',num2str(Regression_coef(8)^2, '% 5.2f')) )
    text(ax(2),axes(1)+0.5,axes(2)-0.25,sprintf('p = %s',num2str(P_vals(8), '% 5.3f')) )
    text(ax(2),axes(1)+0.5,axes(2)-0.5,sprintf('N = %s',num2str(Nr_animals)) )
    axis(ax(2), axes )
    xlabel(ax(2),'Predicted c Integrate Rewards - Sensory-dependent LR')
    legend(sequence)
    legend('boxoff')
    pbaspect(ax(2),[1 1 1]);


    arrayfun( @(i) plot(ax(3),Pred_c_IR_SLR_RR_BT_centered(row,i),Pred_c_opt_BT_centered(row,i),...  % Plot avg last 3 sessions
    'o','MarkerFaceColor',colors(i,1:3),'MarkerEdgeColor',colors(i,1:3),'MarkerSize',5.2 ),...
    1:7);
    line(ax(3),[-3 3],[-3 3],'Color','black','LineStyle','--')
    text(ax(3),axes(1)+0.5,axes(2)+0.25,sprintf('r = %s',num2str(Regression_coef(9), '% 5.2f')) )
    text(ax(3),axes(1)+0.5,axes(2)-0.0,sprintf('r^2 = %s',num2str(Regression_coef(9)^2, '% 5.2f')) )
    text(ax(3),axes(1)+0.5,axes(2)-0.25,sprintf('p = %s',num2str(P_vals(9), '% 5.3f')) )
    text(ax(3),axes(1)+0.5,axes(2)-0.5,sprintf('N = %s',num2str(Nr_animals)) )
    axis(ax(3),axes )
    xlabel(ax(3),'Predicted c Integrate Rewards - Sensory-dependent LR - Relative Rew. Diffs')
    legend(sequence)
    legend('boxoff')
    pbaspect(ax(3),[1 1 1]);

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

%% function [All_InfCxC_means] = SteadyStateCritCondition(Dataframe,InferredC_OCPS,UpallBDataKDB,ChangeOfCond,SessionsToConsider,SessionsToConsiderforBT,sequence)
function [All_InfCxC_means] = SteadyStateCritCondition(Dataframe,InferredC_OCPS,UpallBDataKDB,ChangeOfCond,SessionsToConsider,SessionsToConsiderforBT,sequence)
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
order = 1;
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
Updated_sequence(strcmp(sequence,'RPL'))=4;
Updated_sequence(strcmp(sequence,'RPR'))=5;
Updated_sequence(strcmp(sequence,'SOL'))=6;
Updated_sequence(strcmp(sequence,'SOR'))=7;
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
%% function [...] = compute_correlations_all_models(Dataframe,StimulusMeans,All_InfCxC_means_BT_centered)
function [Regression_coef, P_vals,Pred_c_opt_BT_centered, Pred_c_IR_BT_centered,...
          Pred_c_IR_VL_BT_centered, Pred_c_IR_IRO_BT_centered, Pred_c_IRO_BT_centered,...
          Pred_c_IR_VL_NUNP_BT_centered,Totals_C ]=compute_correlations_all_models(Dataframe,StimulusMeans,All_InfCxC_means_BT_centered)
if Dataframe == 0
    % load ('Criteria_SS_OCPS_Rats_Long.mat'); %OLD
    % load('Criteria_Cond_OCPS_Rats_New.mat'); %NEW
    load ('Long_struct')
    load ('KDB_Params_for_Predictions_Long.mat')
    Allw     = Allw';
    Allw_all = Allw_all';
    AllwU    = AllwU';
    Allparam_3b_NUNP = Allparam_3bNUNP;
    Allw_3b_NUNP = Allw_3bNUNP;
    Dataset = Long_schedule;
    Allw_opt = StimulusMeans;

elseif Dataframe == 2
   % load ('Criteria_SS_OCPS_Pigeons.mat'); %OLD
   % load('Criteria_Cond_OCPS_Pigeons_New.mat');% NEW
    load ('Pigeon_struct')
    load ('KDB_Params_for_Predictions_Pigeon.mat')
    Dataset = Pigeon_schedule;
    Allw     = Allw;
    Allw_all = Allw_all;
    AllwU    = Allw_U;
    AllparamU= Allparam_U;
    Allparam_3b_NUNP = Allparam_3bNUNP;
    Allw_3b_NUNP = Allw_3bNUNP;
    Dataset = Pigeon_schedule;
    Allw_opt = StimulusMeans;
end

conditions = [1,2,3,4,5,6,7];

for  iDataset =1:length(Dataset) % all animals

    allBDataLAK  = Dataset(iDataset).allBDataKDB_LAK;
    sequence     = Dataset(iDataset).sequence (1:max(allBDataLAK(:,6)));
    ChangeOfCond = Dataset(iDataset).ChangeOfCond;

    % Optimal
    model = 'optimal';
    param = [];
    w_D = Allw_opt(iDataset,:);
    no_update_no_pullback = 0;
    [Pred_c_opt(iDataset,:)] = criterion_prediction_all_conditions(model,[], w_D, no_update_no_pullback,Dataframe);

    % Model 1a (Integrate Rewards Omissions)
    model = 'model 1a';
    param = AllparamU(iDataset);
    w_D     = AllwU(iDataset,:);
    no_update_no_pullback = 0;
    [Pred_c_IRO(iDataset,:)] = criterion_prediction_all_conditions(model,param, w_D, no_update_no_pullback,Dataframe);
    
    % Model 2 (Integrate Rewards & Rewards Omissions)
    model = 'model 2';
    param = Allparam_all(iDataset);
    w_D     = Allw_all(iDataset,:);
    no_update_no_pullback = 0;
    [Pred_c_IR_IRO(iDataset,:)] = criterion_prediction_all_conditions(model,param, w_D, no_update_no_pullback,Dataframe);
    clearvars w_D

    % Model 3a (Integrate Rewards)
    model = 'model 3a';
    param = Allparam(iDataset);
    w_D     = Allw(iDataset,:);
    no_update_no_pullback = 0;
    [Pred_c_IR(iDataset,:)] = criterion_prediction_all_conditions(model,param, w_D, no_update_no_pullback,Dataframe);

  if Dataframe == 2
    % Model 3b Reduced (Integrate Rewards + Var LR ) 
    model = 'model 3b'; % Pigeons only 2 learning rates
    param = Allparam_3b_red(iDataset);
    w_D   = Allw_3b_red(iDataset,:);
    w_D(10)= w_D(6);
    w_D(8)  = w_D(7);   w_D(9)  = w_D(7);
    no_update_no_pullback = 0;
    [Pred_c_IR_VL(iDataset,:)] = criterion_prediction_all_conditions(model,param, w_D, no_update_no_pullback,Dataframe);

    % Model 3b Reduced (Integrate Rewards + Var LR + No-Update No Pull-Back (NUNP) )
    model   = 'model 3b'; % Pigeons only 2 learning rates
    param   = Allparam_3b_red_NUNP(iDataset);
    w_D     = Allw_3b_red_NUNP(iDataset,:);
    w_D(10)= w_D(6);
    w_D(8)  = w_D(7);   w_D(9)  = w_D(7);
    no_update_no_pullback = 1;
    [Pred_c_IR_VL_NUNP(iDataset,:)] = criterion_prediction_all_conditions(model,param, w_D, no_update_no_pullback,Dataframe);
 
  else
   % Model 3b (Integrate Rewards + Var LR )
    model = 'model 3b';
    param = Allparam_3b(iDataset);
    w_D   = Allw_3b(iDataset,:);
    no_update_no_pullback = 0;
    [Pred_c_IR_VL(iDataset,:)] = criterion_prediction_all_conditions(model,param, w_D, no_update_no_pullback,Dataframe);

    % Model 3b (Integrate Rewards + Var LR + No-Update No Pull-Back (NUNP) )
    model = 'model 3b';
    param = Allparam_3b_NUNP(iDataset);
    w_D   = Allw_3b_NUNP(iDataset,:);
    no_update_no_pullback = 1;
    [Pred_c_IR_VL_NUNP(iDataset,:)] = criterion_prediction_all_conditions(model,param, w_D, no_update_no_pullback,Dataframe);

  end
end

% center predictions relative to BT

Pred_c_opt_BT_centered        = Pred_c_opt-Pred_c_opt(:,1);
Pred_c_IR_BT_centered         = Pred_c_IR-Pred_c_IR(:,1);
Pred_c_IR_VL_BT_centered      = Pred_c_IR_VL-Pred_c_IR_VL(:,1);
Pred_c_IR_IRO_BT_centered     = Pred_c_IR_IRO-Pred_c_IR_IRO(:,1);
Pred_c_IRO_BT_centered        = Pred_c_IRO-Pred_c_IRO(:,1);
Pred_c_IR_VL_NUNP_BT_centered = Pred_c_IR_VL_NUNP-Pred_c_IR_VL_NUNP(:,1);

Totals_C        = All_InfCxC_means_BT_centered; % Subtract baseline bias to all criteria
Totals_C_in_line= Totals_C(1:end,1:7);          % Arrange in line


% Compare predictions with actual criteria returned by the OCPS model

[r,p] = corrcoef(Totals_C_in_line(:),Pred_c_opt_BT_centered(:));
Regression_coef(1)=r(2);
P_vals(1) = p(2);
[r,p] = corrcoef(Totals_C_in_line(:),Pred_c_IR_BT_centered(:));
Regression_coef(2)=r(2);
P_vals(2) = p(2);
[r,p] = corrcoef(Totals_C_in_line(:),Pred_c_IRO_BT_centered(:));
Regression_coef(3)=r(2);
P_vals(3) = p(2);
[r,p] = corrcoef(Totals_C_in_line(:),Pred_c_IR_IRO_BT_centered(:));
Regression_coef(4)=r(2);
P_vals(4) = p(2);
[r,p] = corrcoef(Totals_C_in_line(:),Pred_c_IR_VL_BT_centered(:));
Regression_coef(5)=r(2);
P_vals(5) = p(2);
[r,p] = corrcoef(Totals_C_in_line(:),Pred_c_IR_VL_NUNP_BT_centered(:));
Regression_coef(6)=r(2);
P_vals(6) = p(2);

% From here the correlation coefficients indicate optimal vs Integrate
% Rewards

[r,p] = corrcoef(Pred_c_opt_BT_centered(:),Pred_c_IR_BT_centered(:));
Regression_coef(7)=r(2);
P_vals(7) = p(2);
[r,p] = corrcoef(Pred_c_opt_BT_centered(:),Pred_c_IR_VL_BT_centered(:));
Regression_coef(8)=r(2);
P_vals(8) = p(2);
[r,p] = corrcoef(Pred_c_opt_BT_centered(:),Pred_c_IR_VL_NUNP_BT_centered(:));
Regression_coef(9)=r(2);
P_vals(9) = p(2);
end
%% function [c] = criterion_prediction_all_conditions (model,param, w_D, no_update_no_pullback,Dataframe)
function [c] = criterion_prediction_all_conditions(model,param, w_D, no_update_no_pullback,Dataframe)

% This function takes the Stimulus Means and sorts all the inputs necessary
% and runs the function criterion_prediction.m individually for each experimental 
% condition

StimSet = w_D(1:5);
if strcmp(model,'model 3b')
    LRStimSet = w_D(6:10);
end

if Dataframe ==2
    idx1 = [1,3,3,5];
    SMcond1 =  [StimSet(idx1)]; %BT
else
    idx1 = [1,2,4,5];
    SMcond1 =  [StimSet(idx1)]; %BT
end
idx2    =  [1,3,4,5];
SMcond2 =  [StimSet(idx2)]; %CL
idx3    =  [1,2,3,5];
SMcond3 =  [StimSet(idx3)]; %CR
idx4    =  [1,3,5];
SMcond4 =  [StimSet(idx4)]; %LeanL
idx5    =  [1,3,5];
SMcond5 =  [StimSet(idx5)]; %LeanR
idx6    =  [1,3,5];
SMcond6 =  [StimSet(idx6)]; %RichL
idx7    =  [1,3,5];
SMcond7 =  [StimSet(idx7)]; %RichR

if strcmp(model,'model 3b')
    LRcond1 = [LRStimSet(idx1)];
    LRcond2 = [LRStimSet(idx2)];
    LRcond3 = [LRStimSet(idx3)];
    LRcond4 = [LRStimSet(idx4)];
    LRcond5 = [LRStimSet(idx5)];
    LRcond6 = [LRStimSet(idx6)];
    LRcond7 = [LRStimSet(idx7)];
elseif strcmp(model,'model 2')
    LRcond1 = w_D(end)-1:w_D(end);
    LRcond2 = w_D(end)-1:w_D(end);
    LRcond3 = w_D(end)-1:w_D(end);
    LRcond4 = w_D(end)-1:w_D(end);
    LRcond5 = w_D(end)-1:w_D(end);
    LRcond6 = w_D(end)-1:w_D(end);
    LRcond7 = w_D(end)-1:w_D(end);
else
    LRcond1 = [w_D(end)];
    LRcond2 = [w_D(end)];
    LRcond3 = [w_D(end)];
    LRcond4 = [w_D(end)];
    LRcond5 = [w_D(end)];
    LRcond6 = [w_D(end)];
    LRcond7 = [w_D(end)];
end

conditions = [1,2,3,4,5,6,7];

for j = 1:length(conditions)
    switch conditions(j)
        case 1 %BT
            w = [SMcond1,LRcond1];

            p = 1/4*[1 1 1 1];
            class = [1 1 2 2];
            r = [1.0 1.0 1.0 1.0];
            

        case 2 %CL
            w = [SMcond2,LRcond2];
            p = 1/4*[1 1 1 1];
            class = [1 2 1 2];
            r = [1.0 1.0 1.0 1.0];

        case 3 %CR
            w = [SMcond3,LRcond3];
            p = 1/4*[1 1 1 1];
            class = [1 2 1 2];
            r = [1.0 1.0 1.0 1.0];

        case 4 %RPL
            %mu = SMcond4;
            w = [SMcond4,LRcond4];
            p = 1/4*[1 2 1];
            class = [1 2 2];
            r = [1.0 0.25 0.5];

        case 5 %RPR
            w = [SMcond5,LRcond5];
            p = 1/4*[1 2 1];
            class = [1 1 2];
            r = [0.5 0.25 1];

        case 6 %SOL
            %mu = SMcond6;
            w = [SMcond6,LRcond6];
            p = 1/4*[2 1 1];
            class = [1 2 2];
            r = [1.0 1.0 1.0];

        case 7 %SOR
            w = [SMcond7,LRcond7];
            p = 1/4*[1 1 2];
            class = [1 1 2];
            r = [1.0 1.0 1.0];
    end
  
    c(1,j) = criterion_prediction(model, param, w, p, class, r, [], [], no_update_no_pullback);
end
end
%% function [c] = criterion_prediction (model, param, w, p, class, r, c_min, c_max, no_update_no_pullback)
function c = criterion_prediction(model, param, w, p, class, r, c_min, c_max, no_update_no_pullback)
% Computes a prediction of the steady-state criterion for the given model
% and experimental condition.
%
% model can take the following values:
%   ---------------|------------------------------------------------------
%   'optimal'      | optimal SDT criterion that maximizes expected reward
%   'model 1a'     | update on unrewarded trials with the same learning
%                  | rate upsilon for left and right responses
%   'model 1b'     | same as model 1a but with a different learning rate
%                  | for each stimulus
%   'model 2'      | update after all trials with learning rate delta for
%                  | rewarded trials, learning rate upsilon for unrewarded
%                  | trials
%   'model 2b'     | same as model 2 but with a different delta for each
%                  | stimulus and a different upsilon for each stimulus
%   'model 3a      | update on rewarded trials with the same learning
%                  | rate delta for left and right responses
%   'model 3b'     | same as model 3a but with a different learning rate
%                  | for each stimulus
%   ---------------|------------------------------------------------------
% Note that these correspond to the models with the same name in kdbfit2x
% for feedback_mode='treat_unrewarded_as_punished'.
%
% param: leaky-integration parameter in the model
% w:     the other model parameters as fitted by kdbfit2x:
%        - stimulus mean for each stimulus
%        - learning rate(s)
% p:     stimulus presentation probability for each stimulus
% class: stimulus class (1 or 2) for each stimulus
% r:     reward rate for correct responses for each stimulus
%
% mu: Gives the specific stimulus Means of the stimulus present in that condition.
%     This simplifies the preprocessing step.            ADDED BY LUIS.
% w, p, class and r need to be row vectors (1 x number of stimuli).
%
% c_min: lower boundary of the interval in which the criterion should be
% c_max: upper boundary of the interval in which the criterion should be
%
% c_min and c_max can be omitted, the range for the criterion is then
% determined by the stimulus means: (min(mu)-2.5, max(mu)+2.5)
%
% no_update_no_pullback: Defaults to false. When set to true, the
%        criterion does not change on trials without an update step (i.e.
%        the leaky integration factor is 1 for those trials)

%mu = w(1:length(class))';                            %   It was transposed in the original version uploaded by Chrsitina

mu = w(1:length(class));                            

if ~strcmp(model,'optimal')
    gamma = param;
end

if ~exist('c_min','var') || isempty(c_min)
    c_min = min(mu)-2.5;
    c_max = max(mu)+2.5;
end

if ~exist('no_update_no_pullback','var') || isempty(no_update_no_pullback)
    no_update_no_pullback = false;
end

switch model
    case 'optimal'
        c = fminbnd(@(c) -sum(expected_reward(c,mu,p,class,r)),c_min,c_max);
    case 'model 3a'
        delta = repmat(w(end),1,length(class));
        c = fzero(@(c) f_income(c,mu,p,class,r,gamma,delta,no_update_no_pullback),[c_min,c_max]);
    case 'model 3b'
        % delta = w(length(class)+1:end)' ;   % it was transposed in the version uploaded by Christina
        delta = w(length(class)+1:end);
        c = fzero(@(c) f_income(c,mu,p,class,r,gamma,delta,no_update_no_pullback),[c_min,c_max]);
    case 'model 1a'
        upsilon = repmat(w(end),1,length(class));
        c = fzero(@(c) f_unrewarded(c,mu,p,class,r,gamma,upsilon,no_update_no_pullback),[c_min,c_max]);
    case 'model 1b'
        upsilon = w(length(class)+1:end)';
        c = fzero(@(c) f_unrewarded(c,mu,p,class,r,gamma,upsilon,no_update_no_pullback),[c_min,c_max]);
    case 'model 2'
        delta = repmat(w(end-1),1,length(class));
        upsilon = repmat(w(end),1,length(class));
        c = fzero(@(c) f_all(c,mu,p,class,r,gamma,delta,upsilon,no_update_no_pullback),[c_min,c_max]);
    case 'model 2b'
        delta = w(len(class)+1:2*length(class));
        upsilon = w(2*length(class)+1:end);
        c = fzero(@(c) f_all(c,mu,p,class,r,gamma,delta,upsilon,no_update_no_pullback),[c_min,c_max]);
end

end
%% function y = f_income(c,mu,p,class,r,gamma,delta,no_update_no_pullback)
function y = f_income(c,mu,p,class,r,gamma,delta,no_update_no_pullback)
[Rf1,Rf2] = reinforcement(c,mu,p,class,r);
if no_update_no_pullback
    p_update = sum(Rf1+Rf2);
else
    p_update = 1;
end
y = (gamma*c-c)+((Rf1*delta')-(Rf2*delta'))/p_update;
end
%% function y = f_unrewarded(c,mu,p,class,r,gamma,upsilon,no_update_no_pullback)
function y = f_unrewarded(c,mu,p,class,r,gamma,upsilon,no_update_no_pullback)
[Ur1,Ur2] = unrewarded(c,mu,p,class,r);
if no_update_no_pullback
    p_update = sum(Ur1+Ur2);
else
    p_update = 1;
end
y = (gamma*c-c)+((Ur2*upsilon')-(Ur1*upsilon'))/p_update;
end
%% function y = f_all(c,mu,p,class,r,gamma,delta,upsilon,no_update_no_pullback)
function y = f_all(c,mu,p,class,r,gamma,delta,upsilon,no_update_no_pullback)
[Rf1,Rf2] = reinforcement(c,mu,p,class,r);
[Ur1,Ur2] = unrewarded(c,mu,p,class,r);
if no_update_no_pullback
    p_update = sum(Ur1+Ur2+Rf1+Rf2);
else
    p_update = 1;
end
y = (gamma*c-c)+((Rf1*delta')-(Rf2*delta')+(Ur2*upsilon')-(Ur1*upsilon'))/p_update;
end
%% function [Rf1,Rf2] = reinforcement(c,mu,p,class,r)
function [Rf1,Rf2] = reinforcement(c,mu,p,class,r)

Rf1 = zeros(size(class));
for i=find(class==1)
  Rf1(i) = p(i)*r(i)*normcdf(c,mu(i),1);
end

Rf2 = zeros(size(class));
for i = find(class==2)
  Rf2(i) = p(i)*r(i)*(1-normcdf(c,mu(i),1));
end
end
%% function [R1,R2] = response_probabilities(c,mu,p)
function [R1,R2] = response_probabilities(c,mu,p)
R1 = zeros(size(mu));
for i = 1:length(mu)
  R1(i) = p(i) * normcdf(c,mu(i),1);
end
R2 = p-R1;
end
%% function [Ur1,Ur2] = unrewarded(c,mu,p,class,r)
function [Ur1,Ur2] = unrewarded(c,mu,p,class,r)
[Rf1,Rf2] = reinforcement(c,mu,p,class,r);
[R1,R2] = response_probabilities(c,mu,p);
Ur1 = R1-Rf1;
Ur2 = R2-Rf2;
end
%% function Rf = expected_reward(c,mu,p,class,r)
function Rf = expected_reward(c,mu,p,class,r)
[Rf1,Rf2] = reinforcement(c,mu,p,class,r);
Rf = Rf1+Rf2;
end