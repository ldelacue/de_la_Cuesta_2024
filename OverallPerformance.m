%% Overall perfromance (session by session)

% This script takes Rat (Dataframe = 0) or Pigeon (Dataframe=2) data and
% plots session by session performance metrics (p.leftward, p.reward and p.
% withdrawal).

%                                               Luis de la Cuesta Ferrer
%                                               25/09/2024

% Added a second section to generate plots tha track the HR S1
% and FAR S5 of all subjects.                   02/19/2025  
%% PLot reward punishment and, premature responses rates on the top of the data (Figure 2c & Figure S6a)
clear all
Dataframe = 0; % 0 for Long Rats, 1 for Short Rats, 2 for pigeons.

if Dataframe ==0
load("Long_struct.mat")
Dataset  = Long_schedule;
Pigeons = 0;
elseif Dataframe ==2
load("Pigeon_struct.mat")
Dataset  = Pigeon_schedule;
Pigeons = 1;
end

for iDataset=1:length(Dataset)

DataMatrix = Dataset(iDataset).allBDataKDB;
sequence   = Dataset(iDataset).sequence;
[RewxSess] = RewardDensxSess(DataMatrix);
[PremResp] = PremRespxCompute(DataMatrix,'session');
ConditionLabels = sequence(1:max(DataMatrix(:,6)));

xaxis= [1-0.5:1:length(RewxSess(2,:)-0.5)];
l1 = Basic_plotting(DataMatrix,ConditionLabels,[]);
hold all;l2 = plot(xaxis,PremResp(3,:),'color',[0 0.80 1],'Linewidth',2);

if Pigeons
l3 = plot(xaxis,RewxSess(6,:),'color',[0.8 0.8 0],'Linewidth',2);
SmoothRewxSess = [];
else %Rats
    l3 = plot(xaxis,RewxSess(6,:),'color',[0.8 0.8 0],'Linewidth',2);
end

legend([l1 l2 l3], {'Leftward responses','Premature responses','Rewarded trials'},'Location','best')

end
%% Plot change in performance (HR S1 and FAR S5) thoughout exp. sessions (Fig 2d (rats), and FigS6.b (pigeons))

clear all
Dataframe = 2; % 0 for Long Rats, 1 for Short Rats, 2 for pigeons.

if Dataframe ==0
    load("Long_struct.mat")
    Dataset  = Long_schedule;
    Pigeons = 0;
    Nr_sessions = 80;
elseif Dataframe ==2
    load("Pigeon_struct.mat")
    Dataset  = Pigeon_schedule;
    Pigeons = 1;
    Nr_sessions = 120; 
end

for iDataset=1:length(Dataset)
    DataMatrix = Dataset(iDataset).allBDataKDB;
    DataMatrix= DataMatrix(~isnan( DataMatrix(:,1)),:); % Clean aborts in Rats
    DataMatrix= DataMatrix(~DataMatrix(:,3)==0,:);      % Clean aborts in Pigeons

    for iSess = 1:Nr_sessions %number of sessions to plot
        Nr_completed_trials(iDataset,iSess) = size(find(DataMatrix(:,7)==iSess),1);
        Nr_S1_trials (iDataset,iSess)       = size(DataMatrix(DataMatrix(:,1)==1 & DataMatrix(:,7)==iSess,:),1);
        Nr_S5_trials (iDataset,iSess)       = size(DataMatrix(DataMatrix(:,1)==5 & DataMatrix(:,7)==iSess,:),1);
        Hits_S1    (iDataset,iSess)         = size(DataMatrix(DataMatrix(:,1)==1 & DataMatrix(:,3)==1 & DataMatrix(:,7)==iSess,:),1); %count nr of R1 responded S1 stim
        FA_S5     (iDataset,iSess)          = size(DataMatrix(DataMatrix(:,1)==5 & DataMatrix(:,3)==1 & DataMatrix(:,7)==iSess,:),1); %count nr of R1 responded S5 stim

        HR_S1 (iDataset,iSess)  = Hits_S1(iDataset,iSess)/Nr_S1_trials(iDataset,iSess);
        FAR_S5(iDataset,iSess)  = FA_S5(iDataset,iSess)/Nr_S5_trials(iDataset,iSess);
    end
end

% Plot
figure; hold all;

plot(HR_S1','b-','Linewidth',0.5);
l1 = plot(mean(HR_S1,1),'b-','Linewidth',5);
plot(FAR_S5','r-','Linewidth',0.5);
l2 = plot(mean(FAR_S5,1),'r-','Linewidth',5);

legend([l1 l2], {'HR S1','FAR S5'},'Location','best')
axis([0 Nr_sessions+3 0 1])

%% function[RewxSess]=RewardDensxSess(allBDataKDB)
function[RewxSess]=RewardDensxSess(allBDataKDB)
% This function takes allBDataKDB and creates a matrix where columns =
% sessions and

% row 1= Rewards in R1 side
% row 2= Rewards in R2 side
% row 3= Total Rewards
% row 4= Total trials
% row 5= Rewards in R2 vs total rewards
% row 6= Reward Density x session

response  = allBDataKDB(:,3);
responded = not(isnan(response)|response==0);           % Get rid of NaNs before feeding OCPS model
UpallBDataKDB = allBDataKDB(responded,:);


for i=1:max(UpallBDataKDB(:,7))
    
    RewxSess(1,i) = length(find(UpallBDataKDB(:,7)==i & UpallBDataKDB(:,4)==1 & UpallBDataKDB(:,2)==1)); 
    RewxSess(2,i) = length(find(UpallBDataKDB(:,7)==i & UpallBDataKDB(:,4)==1 & UpallBDataKDB(:,2)==2)); 
    RewxSess(3,i) = length(find(UpallBDataKDB(:,7)==i & UpallBDataKDB(:,4)==1)); % 
    RewxSess(4,i) = length(find(UpallBDataKDB(:,7)==i));                         % 
    
end

RewxSess(5,:) = RewxSess(2,:)./RewxSess(3,:); %
RewxSess(6,:) = RewxSess(3,:)./RewxSess(4,:); %
end
%% function [PremResp] = PremRespxCompute(allBDataKDB,unit)
function [PremResp] = PremRespxCompute(allBDataKDB,unit)

%%% This function computes the average premature responses per session or
%%% condition ( unit = 'session' or unit = 'condition' ) using the variable
%%% allBDataKDB.

%%%PremResp(1,:) Absolute number of premature responses
%%%PremResp(2,:) Total number of trials
%%%PremResp(3,:) Rate of Prem Resp

if unit == 'session'
    
    for i = 1: max(allBDataKDB(:,7))

        UpallBDataKDB1 = allBDataKDB(not(isnan(allBDataKDB(:,3))),:); %Clean NaN trials in Rats
        UpallBDataKDB2 = UpallBDataKDB1(~UpallBDataKDB1(:,3)==0,:);         %Clean NaN trials in Rats. These two do not interphere
        
        PremResp(1,i)   = length ( find ( allBDataKDB(:,7)==i ) )-length ( find ( UpallBDataKDB2(:,7)==i ) );
        PremResp(2,i)   = length ( find ( allBDataKDB(:,7)==i ) );
        
    end
    
elseif unit == 'condition'
    
    for i = 1: max(allBDataKDB(:,6))
        PremResp(1,i)   = length ( find ( isnan(allBDataKDB(:,1)) & allBDataKDB(:,6)==i ) );
        PremResp(2,i) = length ( find ( allBDataKDB(:,6)==i ) );
        
    end
end

PremResp(3,:) = PremResp(1,:) ./ PremResp(2,:) ;
end

%% function[h1] =  Basic_plotting(DataMatrix,ConditionLabels,rprobs)
function[h1] =  Basic_plotting(DataMatrix,ConditionLabels,rprobs)
%%% Plots performance (pR2) across conditions without fits or anything. 
%%% Code taken from kdbfit2x.m by Luis de la Cuesta 03/05/22

%   Instructions taken from Kdbfit2x.
%   The DataMatrix has as many rows as there are trials in the experiment.
%   The columns are, in this order:
%   ---|------------------------------------------------------------------
%    1 |    sequence: the stimulus that was shown in each trial
%    2 |    stimulus: the stimulus class (1 or 2)
%    3 |    response: the subject's response (1 or 2)
%    4 | rewardarray: reward if response correct or actual reward
%    5 | punisharray: punishment if resp. incorrect or actual punishment
%    6 |   condition: an integer to identify the condition (for plotting)
%    7 |       block: the block number (for plotting)
%   ---|------------------------------------------------------------------
%   Blocks do not necessarily
%   refer to the blocks as they were done in the experiment! A block is
%   only a set of trials that are averaged to produce a plot! It is assumed
%   that all blocks that belong to one condition were done without any
%   other conditions interleaved. If ConditionLabels is given as a cell
%   array of strings (e.g. {'A','B','C'}) then these labels will be used in
%   the plot.


% Example: 
%       Basic_plotting(allBDataKDB,sequence(1:max(allBDataKDB(:,6))),[]);

xaxis = 'block';  % can be 'block' or 'trials'

% -------------------------------------------------------------------------
% Data preparation
% -------------------------------------------------------------------------
sequence = DataMatrix(:,1);              % stimulus sequence
stimulus = DataMatrix(:,2);              % class 1 or a class 2
response = DataMatrix(:,3);              % was the response 1 or 2
rewardarray = DataMatrix(:,4);           % reward if correct?
punisharray = DataMatrix(:,5);

numtrials = length(sequence);
trialnumber = 1:numtrials;

% get rid of trials in which there was no response (NaN)
responded = not(isnan(response)|response==0);
sequence = sequence(responded);
stimulus = stimulus(responded);
response = response(responded);
rewardarray = rewardarray(responded);
punisharray = punisharray(responded);
trialnumber = trialnumber(responded);

% if there are NaN's in the rewardarray, this probably means that for the
% incorrect trials rewardarray was not recorded---as it's also irrelevant.
% In this case we'll just set the NaN's to 0 (ie no reward). But let's
% check that we only do this for incorrect trials, so no problem arises
% down the road

nans = isnan(rewardarray);
if any(nans & (stimulus==response))
    error('NaN-values in correct trials of rewardarray')
else
    rewardarray(nans)=0;
end

% and similarly for punishments
nans = isnan(punisharray);
if any(nans & not(stimulus==response))
    error('NaN-values in incorrect trials of punisharray')
else
    punisharray(nans)=0;
end

% process variables that are needed for plotting
if size(DataMatrix,2) > 5
    condition = DataMatrix(:,6);         % which condition (integer)
    block = DataMatrix(:,7);             % which block was the trial in
    condition = condition(responded);    % also discard trials with no
    block = block(responded);            % valid response here
    plotting = 1;
    % Blocks should be increasing numbers but only increase by one!
    if any(diff(block)>1)
        warning('There are blocks missing.')
    end
    if any(diff(block)<0)
        error('Blocks are not ordered correctly. They must be increasing!')
    end
    % also, blocks should not include different conditions
    for i=unique(block)'
        if length(unique(condition(block==i))) > 1
            error('Blocks cut across conditions!');
        end
    end
else
    plotting = 0;
end

% construct first part of designmatrix A (see sdtcatfit.m)
% just with stimuli, without income and errors; those are added in the
% generate_predictors function below
sequenceset = unique(sequence);
for i = 1:length(sequenceset)
    % we don't need the actual physical measurement of the stimulus, so we
    % just number them according to their order
    sequence(sequence==sequenceset(i)) = i;
end

% -------------------------------------------------------------------------


h0 = figure;
clf
plot_condition(condition,block,trialnumber,numtrials,ConditionLabels,xaxis);
[h1,raw_data] = plot_data(sequence,response,block,'k');
box off
if strcmpi(xaxis,'block')
    xlabel('block')
else
    xlabel('trials')
end
ylabel('Fraction')
drawnow

% --------------------------------------------------------------------------
function plot_condition(condition,block,trialnumber,numtrials, ...
                        ConditionLabels, xaxis)
conditionset = unique(condition);
n = length(conditionset);
if isempty(ConditionLabels)
    for i=1:n
        ConditionLabels{i}=num2str(i);
    end
end
if not(length(ConditionLabels)==n)
    error('ConditionLabels does not match number of conditions')
end

i = 1;
c = [];
xlabels = {};
lastcondition = condition(1);
for t = 1:length(condition)
    if not(condition(t)==lastcondition) % change detected
        if strcmpi(xaxis,'block')
            xlabels{i} = num2str(block(t)-1); %#ok<AGROW>
        else
            xlabels{i} = num2str(trialnumber(t)-1); %#ok<AGROW>
        end
        c(i) = block(t)-1; %#ok<AGROW>
        plot([c(i),c(i)],[0,1],'k:')
        hold on
        text(c(i)+1,0.95,ConditionLabels{condition(t)})
        i=i+1;
    end
    lastcondition = condition(t);
end
if strcmpi(xaxis,'block')
    xlabels{i}=num2str(block(t));
else
    xlabels{i}=num2str(numtrials);
end
c(i)=block(t);
plot([c(i),c(i)],[0,1],'k:')
set(gca,'xtick',sort(c));
set(gca,'xticklabel',xlabels);
text(1,0.95,ConditionLabels{condition(1)})
end

% --------------------------------------------------------------------------

function [handle,raw_data]=plot_data(sequence,response,block,color,whichStim)
if nargin==4
    whichStim = unique(sequence);
end
q = NaN * ones(max(block),1);
for i=1:(max(block)) % discard last block in plot as wished by Maik! Used to be max(block)-1
    select = (block==i) & ismember(sequence,whichStim);
    n1 = sum(response(select)==1);
    n2 = sum(response(select)==2);
    q(i) = n2./(n1+n2);
end
plot((1:length(q))'-0.5,q,[color, '-'],'Linewidth',2)
handle=plot((1:length(q))'-0.5,q,color);
raw_data = q;
end

end
