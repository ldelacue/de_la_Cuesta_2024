%% Callers Overall Performance - Extended
% THis plots (1) Rew. Dens and (2) Movement/Reaction Times across
% conditions. It also (3) reports a negative contrast
% effect observed in Rats when they enterd the firs Lean condition (Lean
% L).

%                                                   Luis 25/02/2025
%% Caller Rew.Dens across conditions
clear all
load("Long_prime_struct.mat")      % this struct is different because it contains not allBDataKDB but allBData matrixes
Dataset = Long_prime_schedule;
SaveOrder = {'BT','SOL','SOR','RPL','RPR','CL','CR'};
Mean_RD   = zeros(length(Dataset),length(SaveOrder));
LowB_95CI_RD    = zeros(length(Dataset),length(SaveOrder));
HighB_95CI_RD   = zeros(length(Dataset),length(SaveOrder));

for iDataset = 1:length(Dataset)
    DataMatrix   = Dataset(iDataset).allBDataKDB;
    DataMatrix1  = DataMatrix(~isnan(DataMatrix(:,3)),:); % remove Withdrawals
    [RewxSess]   = RewardDensxSess(DataMatrix1);
    RewDensxSess = RewxSess(6,:);
    sequence     = Dataset(iDataset).sequence;
    for iCond = 1:max(DataMatrix1(:,6))
        Trial_ind_in_Cond       = find(DataMatrix1(:,6)==iCond);
        Nr_Rew_in_Cond          = length(find(DataMatrix1(Trial_ind_in_Cond,4)==1) );
        Mean_RD(iDataset,iCond) = Nr_Rew_in_Cond/ length(Trial_ind_in_Cond);
        [phat,pci]              = binofit(Nr_Rew_in_Cond,length(Trial_ind_in_Cond));
        LowB_95CI_RD_2          = Mean_RD(iDataset,iCond)-1.96*msep(Mean_RD(iDataset,iCond),length(Trial_ind_in_Cond));HighB_95CI_RD_2 = Mean_RD(iDataset,iCond)+1.96*msep(Mean_RD(iDataset,iCond),length(Trial_ind_in_Cond));  % normal approximation - should be highly similar to pci
        LowB_95CI_RD (iDataset,iCond) = pci(1); HighB_95CI_RD (iDataset,iCond) = pci(2);

    end

    for iOrder = 1:length(SaveOrder)
        WhichCond   = find ( strcmp(sequence,SaveOrder{iOrder})==1);
        WhichCond   = WhichCond(1);
        CondName    = SaveOrder{iOrder};
        Ordered_Mean_RD(iDataset,iOrder)       = Mean_RD(iDataset,WhichCond);
        Ordered_LowB_95CI_RD(iDataset,iOrder)  = LowB_95CI_RD(iDataset,WhichCond);
        Ordered_HighB_95CI_RD(iDataset,iOrder) = HighB_95CI_RD(iDataset,WhichCond);


    end
end

MeanRD_across_animals = mean(Ordered_Mean_RD,1);
% PLotting

ConditionEnd = [1,3,5];
Distance_Begginig_to_SS = 3;
Distance_Conditions     = 10;
orange                  = [0.8500 0.3250 0.0980];
figure;
xaxis = ones(size(Ordered_Mean_RD));
xaxisMeans = (12.5:10:72.5);
for jCond = 1:size(Ordered_Mean_RD,2)
    xaxis(:,jCond)           = xaxis(:,jCond)*Distance_Conditions*jCond;
end

for iDataset = 1:size(Ordered_Mean_RD,1)
    xaxis(iDataset,:)        = xaxis(iDataset,:)+iDataset;
end

%uncomment if you want errorbars
% errorbar(xaxis,Ordered_Mean_RD,abs( Ordered_LowB_95CI_RD-Ordered_Mean_RD),abs( Ordered_HighB_95CI_RD-Ordered_Mean_RD),"ko")
scatter (xaxis,Ordered_Mean_RD,[],"k.");hold on
scatter (xaxisMeans,MeanRD_across_animals,160,"kx")
hold on

for jCond = 1:size(Ordered_Mean_RD,2)
    if ismember (jCond,ConditionEnd)
        line([xaxis(4,jCond)+3.5, xaxis(4,jCond)+3.5],[0 1.25],'LineStyle','--','color',[.7, .7, .7])
    end
    hold all
end
xticks(xaxis(2,:)+0.5)
xticklabels({'B','RichL','RichR','LeanL', 'LeanR','CL','CR'})
axis square
axis( [ 5 80 0 1] )
ylabel('Reward Density')

%% Caller Reaction (RT) and Movement Times (MT)
clear all
load("Long_prime_struct.mat")
Dataset = Long_prime_schedule;
SaveOrder = {'BT','SOL','SOR','RPL','RPR','CL','CR'}; % This the name in which conditions were initially called.
                                                      % Name on paper is listed below.

for iDataset = 1:length(Dataset)
    DataMatrix  = Dataset(iDataset).allBData;
    DataMatrix1 = DataMatrix(~isnan(DataMatrix(:,9)),:); % remove Withdrawals
    RT          = DataMatrix1(:,7) - DataMatrix1(:,4);
    MT          = DataMatrix1(:,8);
    sequence    = Dataset(iDataset).sequence;
    for iCond = 1:max(DataMatrix1(:,16))
        RT_cond(iDataset,iCond)   =  mean(  RT (find( DataMatrix1(:,16)==iCond)) );
        MT_cond(iDataset,iCond)   =  mean(  MT (find( DataMatrix1(:,16)==iCond)) );
    end

    for iOrder = 1:length(SaveOrder)
        WhichCond   = find ( strcmp(sequence,SaveOrder{iOrder})==1);
        WhichCond   = WhichCond(1); % that way we select the first BT and not the rest. For the rest of conditions it does not apply
        CondName    = SaveOrder{iOrder};
        Ordered_RT(iDataset,iOrder) = RT_cond(iDataset,WhichCond);
        Ordered_MT(iDataset,iOrder) = MT_cond(iDataset,WhichCond);

    end
end

figure;
subplot(2,1,1)
boxplot(Ordered_RT,'Labels',{'BT', 'Rich L', 'Rich R', 'Lean L', 'Lean R', 'CL', 'CR'})
ylabel('seconds')
title('Reaction Times')
subplot(2,1,2)
boxplot(Ordered_MT,'Labels',{'BT', 'Rich L', 'Rich R', 'Lean L', 'Lean R', 'CL', 'CR'})
ylabel('seconds')
title('Movement Times')
%% Caller (negative) contrast effect

% Calculate change in pR2 after a reward as a function of the stimulus identity for Rich L vs Lean L at begginig and end of conditions

% Rats sensitivity to reinforcement was increased upon entering
% Lean L condition, in which all animals consistently experimented
% a negative contrast effect, namely less reward density compared to previous
% conditions. Importantly this effect is not present when entering Rich L/R, Lean R (which always
% followed Lean L) and crucially, is even gone in the latter second half of Lean L.
% So this is really specific to the moment of entering Lean L. This
% shows that this is an adaptive effect and not a constant one (otherwise to the
% inverse correlation between reward densitiy of the condition
% and RT). Literature has usually shown that negative contrast effects
% (going from higher to lower reward densitites regimes)
% drive animals to initiate less trials, run less etc (see "From Positive
% and Negative Contrast Effects Using Delayed Reinforcements", by
% R. L. Mellgren, 1972). That is, in terms of reduced vigor.
% Here we report a qualtiatively different effect whereby, negative contrast effect seem to
% drive faster adaptability (learning?) to condition changes.

clear all
close all

Dataframe = 0; % It does not work for pigeons among other reasons because they experienced random experimental sequences.
if Dataframe == 0
    load("Long_struct.mat")
    Dataset = Long_schedule;
elseif Dataframe == 2
    load("Pigeon_struct.mat")
    Dataset = Pigeon_schedule;
end

CondOfInterest = {'SOL', 'RPL', 'SOR', 'RPR'} ; %correspond to Rich L, Lean L, Rich R and Lean R
SessionstoConsider = 3;
for iDataset = 1:length(Dataset)

    DataMatrix  = Dataset(iDataset).allBDataKDB;
    DataMatrix1 = DataMatrix(~isnan(DataMatrix(:,2)),:); % remove Withdrawals
    sequence    = Dataset(iDataset).sequence;

    % select Trials of interest

    for iCond = 1:max(DataMatrix1(:,6))
        TrialsIndxCond(:,iCond)   = DataMatrix1(:,6)==iCond;
    end

    for iOrder = 1:length(CondOfInterest)
        CondPosinDataset(iOrder)   = find(strcmp(sequence,CondOfInterest{iOrder})==1);
    end

    for iOrder = 1:length(CondOfInterest)
        Datachunk1 = DataMatrix1( TrialsIndxCond(:,CondPosinDataset(iOrder)),:);
        FirstSessions = min(Datachunk1(:,7)) : min(Datachunk1(:,7)) + SessionstoConsider-1;
        LastSessions  = max(Datachunk1(:,7)) - SessionstoConsider+1 : max(Datachunk1(:,7));

        %FirstDatachunk (first sessions in Rich/ Lean)
        FirstDataChunk = Datachunk1( logical( sum(Datachunk1(:,7)==FirstSessions,2)),: );
        [pR1overall, pR1_after_Rew, pR1_after_Stim1Rew, pR1_after_Stim3Rew, pR1_after_Stim5Rew] = compute_change_in_Pr1(FirstDataChunk);

        AllModulation_index_Stim1_pR1_First(iDataset,iOrder) = pR1_after_Stim1Rew-pR1_after_Rew;
        AllModulation_index_Stim3_pR1_First(iDataset,iOrder) = pR1_after_Stim3Rew-pR1_after_Rew;
        AllModulation_index_Stim5_pR1_First(iDataset,iOrder) = pR1_after_Stim5Rew-pR1_after_Rew;

        %LastDatachunk (last sessions in Rich/ Lean)
        LastDataChunk  = Datachunk1( logical( sum(Datachunk1(:,7)== LastSessions,2)),: ); %I stop here, here we have already the trials selectedfor first and last blocks
        [pR1overall, pR1_after_Rew, pR1_after_Stim1Rew, pR1_after_Stim3Rew, pR1_after_Stim5Rew] = compute_change_in_Pr1(LastDataChunk);

        AllModulation_index_Stim1_pR1_Last(iDataset,iOrder) = pR1_after_Stim1Rew-pR1_after_Rew;
        AllModulation_index_Stim3_pR1_Last(iDataset,iOrder) = pR1_after_Stim3Rew-pR1_after_Rew;
        AllModulation_index_Stim5_pR1_Last(iDataset,iOrder) = pR1_after_Stim5Rew-pR1_after_Rew;

    end
    clearvars TrialsIndxCond pR1_after_Stim1Rew pR1_after_Stim3Rew pR1_after_Stim5Rew pR1overall

end
Colorvalues = [0     0.447 0.741; 0.85  0.325 0.098;...
               0.929 0.694 0.125; 0.494 0.184 0.556 ]	;
% stats (Paired t-test)

[h_First_Rich_vs_Lean(1,1),p_First_Rich_vs_Lean(1,1)] = ttest(AllModulation_index_Stim1_pR1_First(:,1), AllModulation_index_Stim1_pR1_First(:,2));
[h_First_Rich_vs_Lean(1,2),p_First_Rich_vs_Lean(1,2)] = ttest(AllModulation_index_Stim1_pR1_First(:,3), AllModulation_index_Stim1_pR1_First(:,4));
[h_First_Rich_vs_Lean(2,1),p_First_Rich_vs_Lean(2,2)] = ttest(AllModulation_index_Stim3_pR1_First(:,1), AllModulation_index_Stim3_pR1_First(:,2));
[h_First_Rich_vs_Lean(2,2),p_First_Rich_vs_Lean(2,2)] = ttest(AllModulation_index_Stim3_pR1_First(:,3), AllModulation_index_Stim3_pR1_First(:,4));
[h_First_Rich_vs_Lean(3,1),p_First_Rich_vs_Lean(3,2)] = ttest(AllModulation_index_Stim5_pR1_First(:,1), AllModulation_index_Stim5_pR1_First(:,2));
[h_First_Rich_vs_Lean(3,2),p_First_Rich_vs_Lean(3,2)] = ttest(AllModulation_index_Stim5_pR1_First(:,3), AllModulation_index_Stim5_pR1_First(:,4));

[h_Last_Rich_vs_Lean(1,1),p_Last_Rich_vs_Lean(1,1)] = ttest(AllModulation_index_Stim1_pR1_Last(:,1), AllModulation_index_Stim1_pR1_Last(:,2));
[h_Last_Rich_vs_Lean(1,2),p_Last_Rich_vs_Lean(1,2)] = ttest(AllModulation_index_Stim1_pR1_Last(:,3), AllModulation_index_Stim1_pR1_Last(:,4));
[h_Last_Rich_vs_Lean(2,1),p_Last_Rich_vs_Lean(2,2)] = ttest(AllModulation_index_Stim3_pR1_Last(:,1), AllModulation_index_Stim3_pR1_Last(:,2));
[h_Last_Rich_vs_Lean(2,2),p_Last_Rich_vs_Lean(2,2)] = ttest(AllModulation_index_Stim3_pR1_Last(:,3), AllModulation_index_Stim3_pR1_Last(:,4));
[h_Last_Rich_vs_Lean(3,1),p_Last_Rich_vs_Lean(3,2)] = ttest(AllModulation_index_Stim5_pR1_Last(:,1), AllModulation_index_Stim5_pR1_Last(:,2));
[h_Last_Rich_vs_Lean(3,2),p_Last_Rich_vs_Lean(3,2)] = ttest(AllModulation_index_Stim5_pR1_Last(:,3), AllModulation_index_Stim5_pR1_Last(:,4));

% PLotting
figure;
ax = arrayfun( @(i) subplot(2,3,i,'NextPlot','add','Box','off'), [1:6] );
xaxis = [2,4,6,8];
row2 = { 'Rich L', 'Lean L', 'Rich R', 'Lean R' };

row1 = {'Stim1R1', 'Stim1R1', 'Stim1R1', 'Stim1R1'};
labelArray = [row1; row2];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
arrayfun( @(i) plot( ax(1), xaxis(1:2),AllModulation_index_Stim1_pR1_First(i,1:2)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(1), xaxis(3:4),AllModulation_index_Stim1_pR1_First(i,3:4)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(1), xaxis(i,:),mean(AllModulation_index_Stim1_pR1_First,1),'k*','MarkerSize',10), 1);
ylabel  (ax(1), 'change in pR1 following a reward - First Three sessions')
xticks  (ax(1), xaxis)
set     (ax(1),'XTickLabel',tickLabels);
axis    (ax(1), [0 10 -0.2 0.2])
axis    (ax(1), 'square')

row1 = {'Stim3R2', 'Stim3R2', 'Stim3R1', 'Stim3R1'};
labelArray = [row1; row2];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
arrayfun( @(i) plot( ax(2), xaxis(1:2),AllModulation_index_Stim3_pR1_First(i,1:2)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(2), xaxis(3:4),AllModulation_index_Stim3_pR1_First(i,3:4)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(2), xaxis(i,:),mean(AllModulation_index_Stim3_pR1_First,1),'k*','MarkerSize',10), 1);
xticks  (ax(2), xaxis)
set     (ax(2),'XTickLabel',tickLabels);
axis    (ax(2), [0 10 -0.2 0.2])
axis    (ax(2), 'square')

row1 = {'Stim5R2', 'Stim5R2', 'Stim5R2', 'Stim5R2'};
labelArray = [row1; row2];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
arrayfun( @(i) plot( ax(3), xaxis(1:2),AllModulation_index_Stim5_pR1_First(i,1:2)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(3), xaxis(3:4),AllModulation_index_Stim5_pR1_First(i,3:4)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(3), xaxis(i,:),mean(AllModulation_index_Stim5_pR1_First,1),'k*','MarkerSize',10), 1);
xticks  (ax(3), xaxis)
set     (ax(3),'XTickLabel',tickLabels);
axis    (ax(3), [0 10 -0.2 0.2])
axis    (ax(3), 'square')

row1 = {'Stim1R1', 'Stim1R1', 'Stim1R1', 'Stim1R1'};
labelArray = [row1; row2];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
arrayfun( @(i) plot( ax(4), xaxis(1:2),AllModulation_index_Stim1_pR1_Last(i,1:2)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(4), xaxis(3:4),AllModulation_index_Stim1_pR1_Last(i,3:4)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(4), xaxis(i,:),mean(AllModulation_index_Stim1_pR1_Last,1),'k*','MarkerSize',10), 1);
ylabel  (ax(4), 'change in pR1 following a reward - Last Three sessions')
xticks  (ax(4), xaxis)
set     (ax(4),'XTickLabel',tickLabels);
axis    (ax(4), [0 10 -0.2 0.2])
axis    (ax(4), 'square')

row1 = {'Stim3R2', 'Stim3R2', 'Stim3R1', 'Stim3R1'};
labelArray = [row1; row2];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
arrayfun( @(i) plot( ax(5), xaxis(1:2),AllModulation_index_Stim3_pR1_Last(i,1:2)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(5), xaxis(3:4),AllModulation_index_Stim3_pR1_Last(i,3:4)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(5), xaxis(i,:),mean(AllModulation_index_Stim3_pR1_Last,1),'k*','MarkerSize',10), 1);
xticks  (ax(5), xaxis)
set     (ax(5),'XTickLabel',tickLabels);
axis    (ax(5), [0 10 -0.2 0.2])
axis    (ax(5), 'square')

row1 = {'Stim3R2', 'Stim3R2', 'Stim3R1', 'Stim3R1'};
labelArray = [row1; row2];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
arrayfun( @(i) plot( ax(6), xaxis(1:2),AllModulation_index_Stim5_pR1_Last(i,1:2)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(6), xaxis(3:4),AllModulation_index_Stim5_pR1_Last(i,3:4)', '.-','Color',Colorvalues(i,:),'LineWidth',1.5), [1:length(Dataset)] );
arrayfun( @(i) plot( ax(6), xaxis(i,:),mean(AllModulation_index_Stim5_pR1_Last,1),'k*','MarkerSize',10), 1);
xticks  (ax(6), xaxis)
set     (ax(6),'XTickLabel',tickLabels);
axis    (ax(6), [0 10 -0.2 0.2])
axis    (ax(6), 'square')

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



responded     = not(isnan(allBDataKDB(:,3)));                               % Get rid of NaNs for Rats
UpallBDataKDB = allBDataKDB(responded,:);
UpallBDataKDB = UpallBDataKDB(~allBDataKDB(:,3)==0,:);                      % Get rid of NaNs for Pigeons. These two do not interfere


for i=1:max(UpallBDataKDB(:,7))
    
    RewxSess(1,i) = length(find(UpallBDataKDB(:,7)==i & UpallBDataKDB(:,4)==1 & UpallBDataKDB(:,2)==1)); 
    RewxSess(2,i) = length(find(UpallBDataKDB(:,7)==i & UpallBDataKDB(:,4)==1 & UpallBDataKDB(:,2)==2)); 
    RewxSess(3,i) = length(find(UpallBDataKDB(:,7)==i & UpallBDataKDB(:,4)==1)); % 
    RewxSess(4,i) = length(find(UpallBDataKDB(:,7)==i));                         % 
    
end

RewxSess(5,:) = RewxSess(2,:)./RewxSess(3,:); %
RewxSess(6,:) = RewxSess(3,:)./RewxSess(4,:); %
end
%% function [sep] = msep(p,n)

function [sep] = msep(p,n)
% function [sep] = msep(p,n)
%
% computes the standard error for proportions p for sample size n
% this should not be confused with the standard deviation of the binomial distribution
% the variance of the binomial distribution is n*p*(1-p) (for a series of Bernoulli trials)
% the variance of a proportion is pq, the standard deviation accordingly sqrt(pq)
%
% if p and n are vectors of the same size, sep will be computed for p(1) and n(1), p(2) and n(2) etc.
%
% INPUT
% p     proportion
% n     total number of samples
%
% Maik C. Stuettgen, July 2014
% the (little) works
if numel(p)~=numel(n)
    error('Vector inputs do not match in size.')
end
sep = sqrt((p.*(1-p))./n);
end
%% function [pR1overall, pR1_after_Rew, pR1_after_Stim1Rew, pR1_after_Stim3Rew, pR1_after_Stim5Rew] = compute_change_in_Pr1(Datamatrix)

function [pR1overall, pR1_after_Rew, pR1_after_Stim1Rew, pR1_after_Stim3Rew, pR1_after_Stim5Rew] = compute_change_in_Pr1(Datamatrix)


RewTrialsInd   = find( Datamatrix(:,4)==1);
RewOmTrialsInd = find( Datamatrix(:,4)==0);

Stim1RewTrialsInd   = find( Datamatrix(:,1)==1 & Datamatrix(:,4)==1);
Stim1RewOmTrialsInd = find( Datamatrix(:,1)==1 & Datamatrix(:,4)==0);
Stim3RewTrialsInd   = find( Datamatrix(:,1)==3 & Datamatrix(:,4)==1);
Stim3RewOmTrialsInd = find( Datamatrix(:,1)==3 & Datamatrix(:,4)==0);
Stim5RewTrialsInd   = find( Datamatrix(:,1)==5 & Datamatrix(:,4)==1);
Stim5RewOmTrialsInd = find( Datamatrix(:,1)==5 & Datamatrix(:,4)==0);

pR1overall             = length( find( Datamatrix(:,3)==1))/size( Datamatrix,1);
pR1_after_Rew          = length( find( Datamatrix(RewTrialsInd(1:end-1)+1,3)==1))       / (length(RewTrialsInd)-1);
pR1_after_RewOm        = length( find( Datamatrix(RewOmTrialsInd(1:end-1)+1,3)==1))     / (length(RewOmTrialsInd)-1);
pR1_after_Stim1Rew     = length( find(Datamatrix(  Stim1RewTrialsInd(1:end-1)+1,3)==1)) / (length(Stim1RewTrialsInd)-1);
pR1_after_Stim1RewOm   = length( find(Datamatrix(Stim1RewOmTrialsInd(1:end-1)+1,3)==1)) / (length(Stim1RewOmTrialsInd)-1);
pR1_after_Stim3Rew     = length( find(Datamatrix(  Stim3RewTrialsInd(1:end-1)+1,3)==1)) / (length(Stim3RewTrialsInd)-1);
pR1_after_Stim3RewOm   = length( find(Datamatrix(Stim3RewOmTrialsInd(1:end-1)+1,3)==1)) / (length(Stim3RewOmTrialsInd)-1);
pR1_after_Stim5Rew     = length( find(Datamatrix(  Stim5RewTrialsInd(1:end-1)+1,3)==1)) / (length(Stim5RewTrialsInd)-1);
pR1_after_Stim5RewOm   = length( find(Datamatrix(Stim5RewOmTrialsInd(1:end-1)+1,3)==1)) / (length(Stim5RewOmTrialsInd)-1);
end
