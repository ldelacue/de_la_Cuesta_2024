%% Model parameters & Model performance
% In this code we plot model fitting performance through BIC comparision
% and additionally we plot diverse set of parameters that result from
% fitting the different SDT-learning models. 

%% Load data (the different parameters fitted from each model and the BICs)
clear all
close all
Dataframe = 2;

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


%% Model comparision
% (a) Basic models (Figure S4.b(Rats) and S6.b(Pigeons))
figure (1)
AllBIC        = [AllBIC_3a,AllBIC_1a,AllBIC_2];
boxplot(AllBIC,'Labels',{'IR','IRO','IR&RO'},'PlotStyle','traditional')
hold on

for jCond = 1:3
    ThisConditionDataPoints = AllBIC(:,jCond);
scatter(jCond*ones(size(ThisConditionDataPoints)).*(1+(rand(size(ThisConditionDataPoints))-0.5)/10),ThisConditionDataPoints,20,'blue','filled')
hold on
end

axis square
if Dataframe  ==0
    ylim([12000,32000])
else
    ylim([20000,28000])
end

% (b) New IR models with learning and steady-state modulations (Figure 5d(Rats) and 6g(Pigeons))
figure (2)

if Dataframe ==0
    AllBIC  = [AllBIC_3a,AllBIC_3a_NUNP,AllBIC_3b,AllBIC_3b_NUNP];
elseif Dataframe ==2
    AllBIC  = [AllBIC_3a,AllBIC_3a_NUNP,AllBIC_3b_red,AllBIC_3b_red_NUNP]; % In pigeon we take the reduced models (the ones used in the paper)
end

RelativeBICs = AllBIC(:,:) - AllBIC(:,4);
boxplot(RelativeBICs(:,1:3),'Labels',{'IR - Full','IR-RD - Full','IR-SDL - Full'},'PlotStyle','traditional')
hold on

for jCond = 1:3
    ThisConditionDataPoints = RelativeBICs(:,jCond);
    scatter(jCond*ones(size(ThisConditionDataPoints)).*(1+(rand(size(ThisConditionDataPoints))-0.5)/10),ThisConditionDataPoints,20,'blue','filled')
    hold on
end

axis square

%% Model parameters
% a) Learning rates and leak terms for initial models IR, IRO and IR&RO
% models (Figure 4b (Rats) & 6f (Pigeons))

addpath('\\uni-mainz.de\dfs\Profiles\Settings\ldelacue\Desktop\code\Project 1 - Schedule to be fitted\Model Fitting')
clear all
Dataframe = 2;
load('StimulusMeans_1crit_sess')

if Dataframe == 0
    load ('AllFittedParams_BICs_4Rats.mat')
elseif Dataframe ==2
    load ('AllFittedParams_BICs_4Pigeons.mat')
end

ParameterstoPlot = [Allw_3a(:,6), 1-Allparam_3a' ,Allw_1a(:,6), 1-Allparam_1a ,Allw_2(:,6), Allw_2(:,7), 1-Allparam_2];
figure;
boxplot(ParameterstoPlot,'Labels',{'d (IR)','1-gamma (IR)','u (IRO)','1-gamma (IRO)','d (IR&RO)','u (IR&RO)','1-gamma (IR&RO)'})
title('Parameters Initial Learning models')
hold on
for jCond = 1:size(ParameterstoPlot,2)
    ThisConditionDataPoints = ParameterstoPlot(:,jCond);
scatter(jCond*ones(size(ThisConditionDataPoints)).*(1+(rand(size(ThisConditionDataPoints))-0.5)/10),ThisConditionDataPoints,20,'blue','filled')
hold on
end
hold on
line([0 10],[0 0],'Color','black','LineStyle','--')
axis square

% b) Learning rates as a function of the Stimulus Means of the best model
% ((Figure 5c (Rats)).

if Dataframe == 0 % Rats
    figure % unscaled version
    for i = 1:size(Allw_3b_NUNP,1)
        plot(Allw_3b_NUNP(i,1:5),Allw_3b_NUNP(i,6:10),'.','MarkerSize',16)
        hold on
    end
    ylabel('Raw Learning rates')
    xlabel('Distance to category boundary (KDB)')
    title('Stimulus-dependent learning rates (IR-SLR-RD)')
    axis square

    figure % scaled version
    for i = 1:size(Allw_3b_NUNP,1)
        ScaledLearnRates = Allw_3b_NUNP(i,6:10)./max(Allw_3b_NUNP(i,6:10));
        plot(Allw_3b_NUNP(i,1:5),ScaledLearnRates,'.','MarkerSize',16)
        hold on
    end
    ylabel('Scaled Learning rates')
    xlabel('Distance to category boundary (KDB)')
    title('Discrimination-dependent-dependent learning rates (IR-SLR-RR)')
    axis square

elseif Dataframe ==2 % Pigeons
    figure % unscaled version
    for i = 1:size(Allw_3b_NUNP,1)
        plot(Allw_3b_red_NUNP(i,6:end),'.','MarkerSize',16)
        hold on
    end
    ylabel('Raw Learning rates')
    title('Discrimination-dependent learning rates (IR-SLR(red)-RD)')
    axis square
    xticks ([1 2]);
    xticklabels({'Easy = Stim 1&5','Diff = Stim 2,3,4' })
    axis([0 3 -0.02 0.16])
end

% c) Plot parameters all model versions that are to be discarded because of negative LR (Figure S4.a & S6.d)

if Dataframe ==0
    load('AllFittedParams_BICs_4Rats.mat')
elseif Dataframe ==2
    load('AllFittedParams_BICs_4Pigeons.mat')
end

% We show the learning rates of the Reward Omission components of all the models
Nr_points = [size(Allw_1a,2)-5 + size(Allw_1b,2)-5 + size(Allw_1a_NUNP,2)-5 + size(Allw_1b_NUNP,2)-5 ...
             + size(Allw_2,2)-6 + size(Allw_2b,2)-10];
xaxis = (1:1:Nr_points);
col   = [1,0,1; 1,0,1; 1,0,1; 1,0,1; 1,0,1; 1,0,1; 1,0,1; 1,0,1; 1,0,1; 1,0,1; 1,0,1; 1,0,1;...
         1,1,0; 1,1,0; 1,1,0; 1,1,0; 1,1,0; 1,1,0];
yaxis = [Allw_1a(:,6), Allw_1b(:,6:10), Allw_1a_NUNP(:,6), Allw_1b_NUNP(:,6:10),...
         Allw_2(:,7) , Allw_2b(:,11:15)] ;

yaxis(5,:) = mean (yaxis,1);

figure;% scatter(xaxis(i),yaxis(:,:),'x');hold on
ax = axes('NextPlot','add','Xlim',[0, length(xaxis)+2],'Ylim',[min(yaxis(:))-0.02, max(yaxis(:))+0.02]);
arrayfun ( @(i) plot(ax,xaxis(i),yaxis(1:4,i),'x','MarkerEdgeColor',col(i,1:3),'MarkerSize',6), [1:length(xaxis)] )
arrayfun ( @(i) plot(ax,xaxis(i),yaxis(5,i),'x','MarkerEdgeColor',col(i,1:3),'MarkerSize',12), [1:length(xaxis)] )
line(ax, [0 length(xaxis)+2],[0 0],'Color','black','LineStyle','--')
hold on
xticks(xaxis);
xticklabels({'IRO \upsilon','IRO-SLR \upsilon1','\upsilon2','\upsilon3','\upsilon4','\upsilon5'...
            ,'IRO-RD \upsilon','IRO-SLR-RD\upsilon1','\upsilon2','\upsilon3','\upsilon4','\upsilon5'...
            ,'IR&RO\upsilon','IR&RO-SLR\upsilon1','\upsilon2','\upsilon3','\upsilon4','\upsilon5'});
title ('Unrewarded-learning parameters ')
