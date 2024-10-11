%% Fit_Learning_Models
% This code fits response probabilities with a SDT-learning model of choice to all
% subjects (Rats, Dataframe = 0 or Pigeons, Dataframe = 2). Additionally it
% saves all the fitted parameters under KDB_Fitted_Parameters. It requires
% functions kdbfit2x.m and sdtfit.m that are not defined in this program.

% Main possible models included in the manuscript are: 

% Integrate Rewards                                                             Model = 'model 3a'
% Integrate Reward Omissions                                                    Model = 'model 1a'
% Integrate Rewards & Rew. Omissions                                            Model = 'model 2'
% Integrate Rewards w/ Stim. Learning rates                                     Model = 'model 3b'
% Integrate Rewards w/ Stim. Learning rates (reduced)                           Model = 'model 3b (red)'
% Integrate Rewards w/ Relative Rew.Rates                                       Model = 'model 3a NUNP'
% Integrate Rewards w/ Stim. Learning rates & Relative Rew.Rates                Model = 'model 3b NUNP'
% Integrate Rewards w/ Stim. Learning rates (reduced) & Relative Rew.Rates      Model = 'model 3b (red) NUNP'

%% Fit a model to a whole dataset
clear all
Model    = 'model 3b NUNP';
Dataframe = 0;       % 0 For Rat Long schedule, 2 for Pigeons 
[KDB_Fitted_Parameters] = BatchPlotFits (Dataframe,Model);

%% [KDB_Fitted_Parameters] = BatchPlotFits (Dataframe,Model)
function [KDB_Fitted_Parameters] = BatchPlotFits (Dataframe,Model)
% This function takes the Schedule of choice (0 for Rats Long, , 2 for Pigeons to fit the desired Model by the kdbfit2x.m. It requires
% therefore this program)
                                                                                                %Luis to fit his codes 01/09/23
if Dataframe == 0
    load ('Long_struct.mat')
    Dataset = Long_schedule;

elseif Dataframe ==2
    load ('Pigeon_struct.mat')
    Dataset = Pigeon_schedule;
end

% Reconstruct criterion by fittingt KDB model and save the data
for iDataset = 1:length(Dataset)

    allBDataLAK = Dataset(iDataset).allBDataKDB_LAK;
    sequence    = Dataset(iDataset).sequence (1:max(allBDataLAK(:,6)));

    if strcmpi (Model, 'model 1a') || strcmpi (Model, 'model 1b')
        allBDataLAK(:,5) = 0; % Needed to run unrewarded learning model
        [w,nll,A,param,cval,pS2,BIC,raw_data,fitted_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(allBDataLAK,Model,sequence(1:max(allBDataLAK(:,6))),[],[],[],'treat_unrewarded_as_punished','shuffled');

    elseif strcmpi (Model, 'model 2')
        allBDataLAK(:,5) = 0; % Needed to run unrewarded learning model
        [w,nll,A,param,cval,pS2,BIC,raw_data,fitted_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(allBDataLAK,Model,sequence(1:max(allBDataLAK(:,6))),[],[],[],'treat_unrewarded_as_punished','shuffled');

     elseif strcmpi (Model, 'model 3a')
        [w,nll,A,param,cval,pS2,BIC,raw_data,fitted_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(allBDataLAK,Model,sequence(1:max(allBDataLAK(:,6))),[],[],[],[],'shuffled');

    elseif strcmpi (Model, 'model 3a NUNP')
        no_update_no_pullback = 1;
        [w,nll,A,param,cval,pS2,BIC,raw_data,fitted_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(allBDataLAK,'model 3a',sequence(1:max(allBDataLAK(:,6))),[],[],[],[],'shuffled',no_update_no_pullback);

    elseif strcmpi (Model,'model 3b')
        [w,nll,A,param,cval,pS2,BIC,raw_data,fitted_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(allBDataLAK,Model,sequence(1:max(allBDataLAK(:,6))),[],[],[],[],'shuffled');
    
    elseif strcmpi (Model,'model 3b (red)')
        stimulus_delta_mapping = dictionary([1 2 3 4 5], [1 2 2 2 1]);
        [w,nll,A,param,cval,pS2,BIC,raw_data,fitted_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(allBDataLAK,Model,sequence(1:max(allBDataLAK(:,6))),[],[],[],[],'shuffled',[],stimulus_delta_mapping);

    elseif strcmpi (Model,'model 3b NUNP')
        no_update_no_pullback = 1;
        [w,nll,A,param,cval,pS2,BIC,raw_data,fitted_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(allBDataLAK,'model 3b',sequence(1:max(allBDataLAK(:,6))),[],[],[],[],'shuffled',no_update_no_pullback);
    
    elseif strcmpi (Model,'model 3b NUNP (red)')
        stimulus_delta_mapping = dictionary([1 2 3 4 5], [1 2 2 2 1]);
        no_update_no_pullback = 1;
        [w,nll,A,param,cval,pS2,BIC,raw_data,fitted_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(allBDataLAK,'model 3b',sequence(1:max(allBDataLAK(:,6))),[],[],[],[],'shuffled',no_update_no_pullback,stimulus_delta_mapping);
    end
    ylim([0.15, 0.85])
    KDB_Fitted_Parameters(iDataset,:) = [w',param];
end
end

