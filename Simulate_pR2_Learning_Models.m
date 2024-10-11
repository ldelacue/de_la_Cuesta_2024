%% Simulate_Learning_Models

% This code has two parts that run independnetly. In Part 1 the program
% fits and simulates response probabilities with a SDT-learning-model
% of choice (see possibilities below) to all subjects (Rats, Dataframe = 0 or Pigeons, Dataframe = 2). 
% Additionally it saves all the simulations. The first part requires functions 
% kdbfit2x_Sim.m and sdtfit.m that are not defined within this program.
% Part 2 simply loads simulations (nr = 1000) that have already been run and saved for all
% relevant models in order to plot its spread over the real behavioural data.

% A note about simulations: Simulations are run by using the fitted
% parameters on a series of data (stim sequence, potential reward etc.)
% that is unique for each simulation. In each simulation the trials in a
% condition are randomly shuffled.

% These programs contribute to the plotting of the simulations
% shown in Figures 4a, 5a, 6e, S5.

%% Part1. Fit & Simulate Individual Models on either all Rats or all Pigeon Datasets
% Choose model, dataset (variable Schedule) and number of simulations
% Main possible models included in the manuscript are: 

% Integrate Rewards                                                             Model = 'model 3a'
% Integrate Reward Omissions                                                    Model = 'model 1a'
% Integrate Rewards & Rew. Omissions                                            Model = 'model 2'
% Integrate Rewards w/ Stim. Learning rates                                     Model = 'model 3b'
% Integrate Rewards w/ Stim. Learning rates (reduced)                           Model = 'model 3b (red)'
% Integrate Rewards w/ Relative Rew.Rates                                       Model = 'model 3a NUNP'
% Integrate Rewards w/ Stim. Learning rates & Relative Rew.Rates                Model = 'model 3b NUNP'
% Integrate Rewards w/ Stim. Learning rates (reduced) & Relative Rew.Rates      Model = 'model 3b (red) NUNP'

clear all
Model    = 'model 3b NUNP';
Schedule = 0;       % 0 For Rat Long schedule, 2 for Pigeons 
Nr_sim = 3;

if Schedule == 0
    load ('Long_struct.mat')
    Dataset = Long_schedule;
    Subjects = 'Rats';

elseif Schedule ==2
    load ('Pigeon_struct.mat')
    Dataset = Pigeon_schedule;
    Subjects = 'Pigeons';
end


for iDataset = 1:length(Dataset)
    DataMatrix = Dataset(iDataset).allBDataKDB_LAK;
    sequence    = Dataset(iDataset).sequence (1:max(DataMatrix(:,6)));
    allBDataLAK_Error      = DataMatrix;
    allBDataLAK_Error(:,5) = 0;
    no_update_no_pullback  = 1;
     
    for iSim =1:Nr_sim

        if iSim ==1

            if strcmpi (Model, 'model 1a') || strcmpi (Model, 'model 1b')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(allBDataLAK_Error,'model 1a',sequence(1:max(DataMatrix(:,6))),[],[],[],'treat_unrewarded_as_punished','shuffled');

            elseif strcmpi (Model, 'model 2')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(allBDataLAK_Error,'model 2',sequence(1:max(DataMatrix(:,6))),[],[],[],'treat_unrewarded_as_punished','shuffled');

            elseif strcmpi (Model, 'model 3a')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(DataMatrix,'model 3a',sequence(1:max(DataMatrix(:,6))),[],[],[],[],'shuffled');

            elseif strcmpi (Model, 'model 3a NUNP')
                no_update_no_pullback = 1;
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(DataMatrix,'model 3a',sequence(1:max(DataMatrix(:,6))),[],[],[],[],'shuffled',no_update_no_pullback);

            elseif strcmpi (Model,'model 3b')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(DataMatrix,'model 3b',sequence(1:max(DataMatrix(:,6))),[],[],[],[],'shuffled',no_update_no_pullback,[]);

            elseif strcmpi (Model,'model 3b (red)')
                stimulus_delta_mapping = dictionary([1 2 3 4 5], [1 2 2 2 1]);
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(DataMatrix,'model 3b',sequence(1:max(DataMatrix(:,6))),[],[],[],[],'shuffled',no_update_no_pullback,stimulus_delta_mapping);

            elseif strcmpi (Model,'model 3b NUNP')
                no_update_no_pullback = 1;
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(DataMatrix,'model 3b',sequence(1:max(DataMatrix(:,6))),[],[],[],[],'shuffled',no_update_no_pullback,[]);

            elseif strcmpi (Model,'model 3b NUNP (red)')
                no_update_no_pullback = 1;
                stimulus_delta_mapping = dictionary([1 2 3 4 5], [1 2 2 2 1]);
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(DataMatrix,'model 3b',sequence(1:max(DataMatrix(:,6))),[],[],[],[],'shuffled',no_update_no_pullback,stimulus_delta_mapping);
            end
            FittedGamma = param;
            Allparam(iDataset) = param;
            Allw(iDataset,:) =w;
            AllSimulations_shuffled(:,iSim) = simulated_response;

        else  % for the second simulation we do not need to fit the full model, we already know the parameters (we restrict the grid search of the gamma parameter)

            if strcmpi (Model, 'model 1a') || strcmpi (Model, 'model 1b')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x_Sim(allBDataLAK_Error,'model 1a',sequence(1:max(DataMatrix(:,6))),[],[],[],'treat_unrewarded_as_punished','shuffled');

            elseif strcmpi (Model, 'model 2')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x_Sim(allBDataLAK_Error,'model 2',sequence(1:max(DataMatrix(:,6))),[],[],FittedGamma,[],'treat_unrewarded_as_punished','shuffled');

            elseif strcmpi (Model, 'model 3a')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x_Sim(DataMatrix,'model 3a',sequence(1:max(DataMatrix(:,6))),[],[],FittedGamma,[],[],'shuffled');

            elseif strcmpi (Model, 'model 3a NUNP')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x_Sim(DataMatrix,'model 3a',sequence(1:max(DataMatrix(:,6))),[],[],FittedGamma,[],[],'shuffled',no_update_no_pullback);

            elseif strcmpi (Model,'model 3b')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x_Sim(DataMatrix,'model 3b',sequence(1:max(DataMatrix(:,6))),[],[],FittedGamma,[],[],'shuffled', [], []);

            elseif strcmpi (Model,'model 3b (red)')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x_Sim(DataMatrix,'model 3b',sequence(1:max(DataMatrix(:,6))),[],[],FittedGamma,[],[],'shuffled', [], stimulus_delta_mapping);

            elseif strcmpi (Model,'model 3b NUNP')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x_Sim(DataMatrix,'model 3b',sequence(1:max(DataMatrix(:,6))),[],[],FittedGamma,[],[],'shuffled',no_update_no_pullback,[]);

            elseif strcmpi (Model,'model 3b NUNP (red)')
                [w,nll,A,param,cval,pS2,BIC,raw_data,simulated_data,simulated_response,reward_densities] = kdbfit2x_Sim(DataMatrix,'model 3b',sequence(1:max(DataMatrix(:,6))),[],[],FittedGamma,[],[],'shuffled',no_update_no_pullback,stimulus_delta_mapping);

            end
            AllSimulations_shuffled(:,iSim) = simulated_response;

        end

        Dataset
        iDataset
        iSim

        close all
    end
    Name1 = ['',num2str(Nr_sim)','Sim_','Subject',num2str(iDataset)','_'];
    Name2 = convertCharsToStrings(Model)';
    Name3 = convertCharsToStrings(Subjects)';
    Name  = append (Name1,Name2,Name3);
    save(Name, 'AllSimulations_shuffled','FittedGamma','w');
    clearvars 'AllSimulations_shuffled' 'FittedGamma' 'w' 'no_update_no_pullback' 'stimulus_delta_mapping';
end

%% Part2. PLot already saved Simulations over behavioural data. Subject by subject
clear all
close all
load('1000Sim_All_Models_Pigeon850.mat') %choose dataset
Dataframe = 2; % Is it a Rat (Dataframe = 0) or a Pigeon (Dataframe = 2)?
SubjectNR = 2; % Within Rats or Pigeons, which dataset in ordinal number is it?

if Dataframe == 0
    load ('Long_struct.mat')
    Dataset = Long_schedule;

elseif Dataframe ==2
    load ('Pigeon_struct.mat')
    Dataset = Pigeon_schedule;
end

DataMatrix       = Dataset(SubjectNR).allBDataKDB_LAK;
ConditionLabels  = Pigeon_schedule(SubjectNR).sequence (1:max(DataMatrix(:,6)));
NrSTD            = 1;
Sessions         = 1:max(DataMatrix(:,7));
Type             = 1;


Simulations = AllSimulations1a_shuffled;
Plot_Spread_Simulations(DataMatrix,ConditionLabels,Simulations,Type)
title('Model IRO')

Simulations = AllSimulations2_shuffled;
Plot_Spread_Simulations(DataMatrix,ConditionLabels,Simulations,Type)
title('Model IR&RO')

Simulations = AllSimulations3a_shuffled;
Plot_Spread_Simulations(DataMatrix,ConditionLabels,Simulations,Type)
title('Model IR')

Simulations = AllSimulations3a_NUNP_shuffled;
Plot_Spread_Simulations(DataMatrix,ConditionLabels,Simulations,Type)
title('Model IR-Rel.Rew')

Simulations = AllSimulations3b_shuffled;
Plot_Spread_Simulations(DataMatrix,ConditionLabels,Simulations,Type)
title('Model IR-SLR')

Simulations = AllSimulations3b_NUNP_shuffled;
Plot_Spread_Simulations(DataMatrix,ConditionLabels,Simulations,Type)
title('Model IR-SLR-Rel.Rew')

if Dataframe == 2
Simulations = AllSimulations3b_red_shuffled;
Plot_Spread_Simulations(DataMatrix,ConditionLabels,Simulations,Type)
title('Model IR-SLR(red)')

Simulations = AllSimulations3b_red_NUNP_shuffled;
Plot_Spread_Simulations(DataMatrix,ConditionLabels,Simulations,Type)
title('Model IR-SLR(red)-Rel.Rew')
end

%% function Plot_Spread_Simulations(DataMatrix,ConditionLabels,Simulations,Type)
function Plot_Spread_Simulations(DataMatrix,ConditionLabels,Simulations,Type)

Type  = 1; % 0 = median and 25-75 Interquartile range and 1 = Mean +- 1 STD;
Mean  = mean(Simulations,2);
STD   = std(Simulations,0,2);
NrSTD = 1;

Median = median(Simulations,2);
Range= prctile(Simulations(:,:)',[25,75]); % returns interquartile range. It is best because with STD we get negative values sometimes
%% Plot behavioral Data (pR2)

Basic_plotting(DataMatrix,ConditionLabels)
axis([0 max(DataMatrix(:,7))+2 0.15 0.85])                       % Here we edit the axes
hold all;

if Type ==0
Make_shade(Median,Range,Type)
else
Make_shade(Mean,STD,Type,NrSTD)
end

%--------------------------------------------

    function Make_shade(varargin)
        if Type ==0
        y = Median; % your median vector;
        x = (1:numel(y))-0.5;
        curve1 = Range(1,:)';
        curve2 = Range(2,:)';

        else
        y = Mean; % your mean vector;
        x = (1:numel(y))-0.5;
        std_dev = NrSTD*STD;
        curve1 = y + std_dev;
        curve2 = y - std_dev;
        end

        x2 = [x, fliplr(x)];
        inBetween = [curve2', fliplr(curve1')];
        h = fill(x2, inBetween, 'b');
        h.FaceAlpha = 0.2 ; % for 20% transparent
        hold on;
        plot(x, y, 'color','b', 'LineWidth', 2);

    end


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
