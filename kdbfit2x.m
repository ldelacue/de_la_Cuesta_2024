function [w,nll,A,param,cval,pS2,BIC,raw_data,fitted_data,simulated_data,simulated_response,reward_densities] = kdbfit2x(DataMatrix,model,ConditionLabels,rprobs,StimulusMeans,scale_reward,feedback_mode,simulation_type,no_update_no_pullback,stimulus_delta_mapping)
% function [w,nll,A,param,cval,pS2,BIC,raw_data,fitted_data,simulated_data,simulated_response] = kdbfit2x(DataMatrix,model,ConditionLabels,rprobs)
%
% same as kdbfit, but engineered to extract optimal criterion

% SAME AS KDBFIT2, BUT USES ALL BLOCKS - KDBFIT2 DISCARDS THE LAST BLOCK (CODE SAYS I WANTED THIS TO BE)
% CHANGES APPLY TO NESTED FUNCTIONS PLOT_DATA AND PLOT_FIT
% CODE HAS BEEN MODIFIED TO TOLERATE OMISSION OF RPROBS AND TO TOLERATE NON-CONTIGUOUS CONDITION LABELS
%
%
% FRANK: applies only when stimuli are equiprobable throughout - would e.g. not work in Nature Behavior project
% FRANK: Leider hast Du fr?her oft die gleiche Condition mehrmals dem Tier gezeigt, so dass vielleicht die Conditions
% eine Abfolge 4 3 1 2 4 hatten, mit der Baseline am Anfang und am Ende. Dann funktioniert Dein Code nicht.
% Dein Code geht nur, wenn jede Bedingung auch nur einmal vorkommt. Deshalb ist mein Code in plot_condition so h??lich.
%
% input argument rprobs should contain as many rows as experimental conditions (incl. baseline)
% and as many columns as there were stimuli; importantly, rprobs can be omitted for the last condition if rprobs
% are identical to the first condition (usually given when testing began and ended with baseline testing)
% new output arguments cval and pS2 (optimal criterion and pS2 values for each condition)
%
% [w,nll,A,param] = kdbfit2(DataMatrix)
% [w,nll,A,param] = kdbfit2(DataMatrix,model)
% [w,nll,A,param] = kdbfit2(DataMatrix,model,ConditionLabels)
% [w,nll,A,param] = kdbfit2(DataMatrix,model,ConditionLabels,rprobs,StimulusMeans)
%
%   Fit the Kac-Dorfman-Biderman model (kdb) with the leaky integrator
%   modification described by Stuettgen et al. (2013) to single trial data.
%   Instead of learning from correct trials / errors as in the original KDB
%   model, this model learns only from rewarded/punished trials.
%
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
%   The last two are only used for plotting and if they are dropped from
%   the DataMatrix no plot will be generated. Otherwise all responses in
%   one block will be averaged for the plot. Blocks do not necessarily
%   refer to the blocks as they were done in the experiment! A block is
%   only a set of trials that are averaged to produce a plot! It is assumed
%   that all blocks that belong to one condition were done without any
%   other conditions interleaved. If ConditionLabels is given as a cell
%   array of strings (e.g. {'A','B','C'}) then these labels will be used in
%   the plot.
%
%   Possible options for model are as the same as described in Dorfman &
%   Biderman (1971):
%   ---------------|------------------------------------------------------
%   'general model'| general case of kdb; learn on rewarded and punished
%                  | trials with no constraints on the learning parameters
%   'model 1'      | constrain the general model so that we only learn on
%                  | punished trials
%   'model 1a'     | same as model 1 but force same learning rate for left
%                  | and right responses (same as model by Kac)
%   'model 1b'     | same as model 1a but with a different learning rate
%                  | for each stimulus
%   'model 2'      | make sure that rewarded trials have same learning
%                  | rate and that punished have the same learning
%                  | rate; but rewarded and punished can be
%                  | different
%   'model 2b'     | same as model 2 but with a different learning rate for
%                  | each stimulus
%   'model 3'      | only learn on rewarded trials; in the original kdb
%                  | model this does not make sense because this will lead
%                  | to exclusive choice---but we have a forgetting term
%                  | that pulls the criterion back to zero
%   'model 3a'     | same as model 3 but constrain learning rate to be the
%                  | same for the two response options
%   'model 3b'     | same as model 3a but with a different learning rate
%                  | for each stimulus
%   ---------------|------------------------------------------------------
%   If model is not specified then the model from Stuettgen et al. (2013)
%   will be used, which is 'model 3a'. Note also that we split up the
%   trials in rewarded and punished according to the two vectors in the
%   DataMatrix (columns 4 and 5). Hence, nothing happens on trials where
%   there was no reward and no punishment. The code assumes that only
%   correct responses can be rewarded and only errors can be punished.
%
%   If StimulusMeans is provided, those will be used as the means for the
%   stimuli instead of fitting the means as free parameters. Another
%   parameter for the initial criterion will be fitted in this case.
%
%   Possible options for scale_reward are 'received' and 'programmed'. If
%   either of those is given, before fitting the model the rewards in each
%   condition will be scaled by 1/reward_density, where reward_density is
%   the reward received (for 'received') or possible to gain (for
%   'programmed') averaged over all trials in that condition.
%
%   feedback_mode can have three values (defaults to 'standard'):
%       'standard': The criterion is not updated in trials where there is 
%           neither reward nor punishment.
%       'treat_unrewarded_as_punished': All trials where no reward was
%           given are treated as if they were punished. Do not provide any
%           values in punisharray (set them all to 0 or NaN).
%       'treat_unpunished_as_rewarded': All trials where no punishment was
%           given are treated as if they were rewarded. Do not provide any
%           values in rewardarray (set them all to 0 or NaN).
%
%   simulation_type:
%       'experienced': Simulates the responses for an agent with the
%        same sequence as encountered in DataMatrix.
%       'shuffled'   : Simulates the responses with a Condition-wise
%        shuffled trial order.
%
%   no_update_no_pullback: Defaults to false. When set to true, the
%       criterion does not change on trials without an update step (i.e.
%       the leaky integration factor is 1 for those trials).
%
%   stimulus_delta_mapping: For models 1b, 2b and 3b, it is possible to
%       group several stimuli together to use the same learning rate delta.
%       To do so, a dictionary with stimulus numbers as keys and delta
%       indices as values needs to be provided here.
%       E.g. stimulus_delta_mapping = dictionary([1 2 3 4 5], [1 2 2 2 1])
%
%   The function returns the best fitting parameters for the stimulus means
%   in w, the negative log likelihood in nll, the designmatrix that was
%   used in the glm as A, and the temporal integration parameter param
%   (gamma in the leaky integrator of St?ttgen et al.).
%   If scale_reward was 'received' or 'programmed', additionally the
%   respective reward densities for all conditions are returned.
%
%   References:
%       Dorfman, D. D. & Biderman, M. (1971). A Learning Model for a
%         Continuum of Sensory States. Journal of Mathematical Psychology,
%         8, 264-284.
%       Stüttgen et al. (2013). Suboptimal Criterion Setting in a
%         Perceptual Choice Task with Asymmetric Reinforcement. Behavioural
%         Processes, 96, 59-70.
%
%   Version 4.6 --- May 2024, Frank Jäkel, frank.jaekel@tu-darmstadt.de
%                             Christina Koß, christina.koss@tu-darmstadt.de
%
DEBUG = 0; % show a plot for fitting the leaky integration parameter
% In Stüttgen et al. (2013) we used 'model 3a' of Dorfman & Biderman (1971)
% so let's use this as a default if the model is not specified
if nargin == 1; model = 'model 3a'; end
if nargin < 3; ConditionLabels = {}; end
if ~exist('StimulusMeans','var') || isempty(StimulusMeans); StimulusMeans = false; end
if ~exist('scale_reward','var') || isempty(scale_reward); scale_reward = ''; end
if ~exist('feedback_mode','var') || isempty(feedback_mode); feedback_mode = 'standard'; end
if ~exist('simulation_type','var') || isempty(simulation_type); simulation_type = 'experienced'; end
if ~exist('no_update_no_pullback','var') || isempty(no_update_no_pullback); no_update_no_pullback = false; end
if ~exist('stimulus_delta_mapping','var') || isempty(stimulus_delta_mapping)
    sequenceset = unique(DataMatrix(:,1));
    stimulus_delta_mapping = dictionary(sequenceset, sequenceset);
end

integrator = 'leaky integrator';   % as in Stüttgen et al. (2013)
%integrator = 'moving average';    % for comparison
%integrator = 'cumulative sum';    % the standard kdb model

xaxis = 'block';  % can be 'block' or 'trials'

if no_update_no_pullback
    if strcmp(integrator,'leaky integrator')
        integrator = 'leaky integrator on update';
    else
        error('no_update_no_pullback only makes sense when the leaky integrator is used.')
    end
end

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

% If feedback mode is treat_unrewarded_as_punished
% (treat_unpunished_as_rewarded), no punishments (rewards) should be
% explicitly provided.
if strcmpi(feedback_mode,'treat_unrewarded_as_punished')
    if all(punisharray==0)
        punisharray(:) = nan;
    end
    if ~all(isnan(punisharray))
        error('Explicit punishment was provided. It does not make sense to treat unrewarded trials as punished.')
    end
end
if strcmpi(feedback_mode,'treat_unpunished_as_rewarded')
    if all(rewardarray==0)
        rewardarray(:) = nan;
    end
    if ~all(isnan(rewardarray))
        error('Explicit rewards were provided. It does not make sense to treat unpunished trials as rewarded.')
    end
end

% if there are NaN's in the rewardarray, this probably means that for the
% incorrect trials rewardarray was not recorded---as it's also irrelevant.
% In this case we'll just set the NaN's to 0 (ie no reward). But let's
% check that we only do this for incorrect trials, so no problem arises
% down the road
if strcmpi(feedback_mode,'standard') || strcmpi(feedback_mode,'treat_unrewarded_as_punished')
    nans = isnan(rewardarray);
    if any(nans & (stimulus==response))
        error('NaN-values in correct trials of rewardarray')
    else
        rewardarray(nans)=0;
    end
end
    
% and similarly for punishments
if strcmpi(feedback_mode,'standard') || strcmpi(feedback_mode,'treat_unpunished_as_rewarded')
    nans = isnan(punisharray);
    if any(nans & not(stimulus==response))
        error('NaN-values in incorrect trials of punisharray')
    else
        punisharray(nans)=0;
    end
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

% scale reward with 1/reward density
if strcmp(scale_reward,'programmed') || strcmp(scale_reward,'received')
    [rewardarray, reward_densities] = scale_reward_with_density(stimulus, response, rewardarray, condition, scale_reward);
else
    reward_densities = [];
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
if ~StimulusMeans
    A = zeros(length(sequence),length(sequenceset));
    for i = 1:length(sequenceset)
        % In the designmatrix we set the corresponding column to have a -1
        % whenever the stimulus was shown
        A(sequence==i,i) = -1;
    end
    % do NOT add a constant bias term for all trials since we already have a
    % free parameter for every stimulus
    % A = [A ones(length(sequence),1)];
    b = zeros(length(sequence),1);
else
    A = ones(length(sequence),1);
    b = nan(length(sequence),1);
    for i = 1:length(sequenceset)
        b(sequence==i) = - StimulusMeans(i);
    end
end

disp('------------------------------------------------------------------')
switch lower(model)
    case 'general model'
        disp('Fitting general model: reward (2 params) and punishment (2 params)')
        if all(rewardarray == 0)
            error('There was no reward given. It does not make sense to fit the general model.')
        end
        if all(punisharray == 0)
            error('There was no punishment given. It does not make sense to fit the general model.')
        end
        for i = unique(response)'
            if all(rewardarray(response==i)==0)
                error('There was no reward given for one of the responses. It does not make sense to have 2 params for reward learning.')
            end
            if all(punisharray(response==i)==0)
                error('There was no punishment given for one of the responses. It does not make sense to have 2 params for punishment learning.')
            end
        end
    case 'model 1'
        disp('Fitting model 1: only punishment (2 params), reward column ignored')
        if all(punisharray == 0)
            error('There was no punishment given. It does not make sense to fit model 1.')
        end
        for i = unique(response)'
            if all(punisharray(response==i)==0)
                error('There was no punishment given for one of the responses. It does not make sense to have 2 params for punishment learning.')
            end
        end
    case 'model 1a'
        disp('Fitting model 1a: only punishment (1 param), reward column ignored')
        if all(punisharray == 0)
            error('There was no punishment given. It does not make sense to fit model 1a.')
        end
    case 'model 1b'
        disp('Fitting model 1b: only punishment, one per stimulus (n_stim params), reward column ignored')
        if all(punisharray == 0)
            error('There was no punishment given. It does not make sense to fit model 1b.')
        end
        for i = unique(sequence)'
            if all(punisharray(sequence==i)==0)
                error('There was no punishment given for one of the stimuli. It does not make sense to have a param for each stimulus.')
            end
        end
    case 'model 2'
        disp('Fitting model 2: reward (1 param) and punishment (1 param)')
        if all(rewardarray == 0)
            error('There was no reward given. It does not make sense to fit model 2.')
        end
        if all(punisharray == 0)
            error('There was no punishment given. It does not make sense to fit model 2.')
        end
    case 'model 2b'
        disp('Fitting model 2b: reward (1 param per stimulus) and punishment (1 param per stimulus)')
        if all(rewardarray == 0)
            error('There was no reward given. It does not make sense to fit model 2b.')
        end
        if all(punisharray == 0)
            error('There was no punishment given. It does not make sense to fit model 2b.')
        end
    case 'model 3'
        disp('Fitting model 3: only reward (2 params), punishment column ignored')
        if all(rewardarray == 0)
            error('There was no reward given. It does not make sense to fit model 3.')
        end
        for i = unique(response)'
            if all(rewardarray(response==i)==0)
                error('There was no reward given for one of the responses. It does not make sense to have 2 params for reward learning.')
            end
        end
    case 'model 3a'
        disp('Fitting model 3a: only reward (1 param), punishment column ignored')
        if all(rewardarray == 0)
            error('There was no reward given. It does not make sense to fit model 3a.')
        end
    case 'model 3b'
        disp('Fitting model 3b: only reward, one per stimulus (n_stim params), punishment column ignored')
        if all(rewardarray == 0)
            error('There was no reward given. It does not make sense to fit model 3b.')
        end
        for i = unique(sequence)'
            if all(rewardarray(sequence==i)==0)
                error('There was no reward given for one of the stimuli. It does not make sense to have a param for each stimulus.')
            end
        end
    otherwise
        error('Unknown model restriction')
end
disp('------------------------------------------------------------------')

% plot the data right away so we have something to look at while we're
% waiting for the fitting to finish
if plotting
    h0 = figure;
    clf
    plot_condition(condition,block,trialnumber,numtrials,ConditionLabels,xaxis);
    [h1,raw_data] = plot_data(sequence,response,block,'b');
    box off
    if strcmpi(xaxis,'block')
        xlabel('block')
    else
        xlabel('trials')
    end
    ylabel('fraction S_2 responses')
    drawnow
end

% -------------------------------------------------------------------------
% Model fitting
% -------------------------------------------------------------------------
tic
switch integrator
    case {'leaky integrator', 'leaky integrator on update'}
        if DEBUG
            params = linspace(0.01,0.99,100);
            params = params(2:end-1);
            params = [logspace(-2,-10,10) params (1-logspace(-2,-10,10))];
            params = sort(params);
            nlls = nan(size(params));
            figure
            clf
            for i = 1:length(params)
                param = params(i);
                nll = kdb_fit(param,stimulus,sequence,response,...
                    rewardarray,punisharray,A,b,model,integrator,...
                    stimulus_delta_mapping);
                nlls(i) = nll;
                plot(params,nlls)
                ylabel('negative log likelihood')
                xlabel('Leaky integration parameter')
                title(lower(model))
                box off
            end
        end
        drawnow
        opt = optimset;
        param = fminbnd(@kdb_fit,0.01,1,opt,...
            stimulus,sequence,response,rewardarray,punisharray,...
            A,b,model,integrator,stimulus_delta_mapping);
        [nll,w,A] = kdb_fit(param,stimulus,sequence,response,...
            rewardarray,punisharray,A,b,model,integrator,...
            stimulus_delta_mapping);
        if DEBUG
            hold on
            plot([param param],ylim,'r')
        end
        k = length(w)+1; % number of free params, +1 for leaky integrator
        fprintf('Negative log likelihood: %6.8g\n', nll)
        fprintf('                    BIC: %6.8g\n\n', ...
            2*nll+k*log(length(sequence)))
        fprintf('Leaky integration: %.3g\n', param)
        fprintf('        half time: %d trials\n', ...
            round(log(0.5)/log(param)))
        fprintf('       100th time: %d trials\n\n', ...
            round(log(0.01)/log(param)))
    case 'moving average'
        params = 10:5:200;
        nlls = nan(size(params));
        figure
        clf
        for i = 1:length(params)
            param = params(i);
            nll = kdb_fit(param,stimulus,sequence, response,...
                rewardarray,punisharray,A,b,model,integrator,...
                stimulus_delta_mapping);
            nlls(i) = nll;
            plot(params,nlls)
            ylabel('negative log likelihood')
            xlabel('moving average window size (trials)')
            drawnow
        end
        [~,i] = min(nlls);
        param = params(i);
        [nll,w,A] = kdb_fit(param,stimulus,sequence,response,...
            rewardarray,punisharray,A,b,model,integrator,...
            stimulus_delta_mapping);
        k = length(w)+1; % number of free parameters, add 1 window size
        fprintf('Negative log likelihood: %6.8g\n', nll)
        fprintf('                    BIC: %6.8g\n\n', ...
            2*nll+k*log(length(sequence)))
        fprintf('Integration window: %d trials\n\n', param)
    case 'cumulative sum'
        param = Inf;
        [nll,w,A] = kdb_fit(param,stimulus,sequence,response,...
            rewardarray,punisharray,A,b,model,integrator,...
            stimulus_delta_mapping);
        k = length(w);
        fprintf('Negative log likelihood: %6.8g\n', nll)
        fprintf('                    BIC: %6.8g\n\n', ...
            2*nll+k*log(length(sequence)))
    otherwise
        error('unknown integrator')
end
% and output parameter vector
if ~StimulusMeans
    fprintf('Parameter vector (%d stimuli, %d learning rates):\n\n',...
        length(sequenceset),length(w)-length(sequenceset))
else
    fprintf('Parameter vector (initial criterion, %d learning rates):\n\n',...
        length(w)-1)
end
disp(w)
BIC = 2*nll+k*log(length(sequence));  % just so that BIC can be output
toc

tic
disp('-------------------------------------------')
disp('Generating simulated data from fitted model')
disp('-------------------------------------------')
if strcmpi(simulation_type,'shuffled')
    [sim_stimulus, sim_sequence, sim_rewardarray, sim_punisharray, ~] = shuffle_data_conditionwise(stimulus,sequence,rewardarray,punisharray,condition);
else
    sim_stimulus = stimulus;
    sim_sequence = sequence;
    sim_rewardarray = rewardarray;
    sim_punisharray = punisharray;
end
[simulated_response, ~, pred1] = simulate_data(sim_stimulus,sim_sequence,sim_rewardarray,sim_punisharray,...
    w,model,integrator,param,StimulusMeans,stimulus_delta_mapping);
% check that predictors calculated in simulate_data are really the same as
% calculated with generate_predictors
pred2 = generate_predictors(sim_stimulus,sim_sequence,simulated_response,sim_rewardarray,sim_punisharray,model,integrator,param,stimulus_delta_mapping);
if not(all(abs(pred1(:)-pred2(:))<10^-11))
    figure
    plot(pred1(:)-pred2(:))
    error('arrgh')
end
toc

if ~isempty(rprobs)
  disp('------------------------------------------------------------------')
  disp('Maik''s code snippet for determining the optimal criterion')
  disp('------------------------------------------------------------------')
  if size(rprobs,2)~=numel(sequenceset),error('size(rprobs,2)~=numel(sequenceset)'),end
  if size(rprobs,1)==max(condition)-1
    rprobs(size(rprobs,1)+1,:)=rprobs(1,:);
    disp('---> appended rprobs')
  end
  for i = 1:size(rprobs,1)
    if StimulusMeans
        mu = StimulusMeans;
    else
        mu = w(1:length(sequenceset));
    end
    cval(i) = fminbnd(@(x) getvORF(x,mu,rprobs(i,:)),min(mu)-3,max(mu)+3); %#ok<AGROW>
    pS2(i)  = 1-mean(normcdf(cval(i),mu,1)); %#ok<AGROW>
  end
  % comment by Frank: pS2(i)... ist es richtig den mean zu nehmen? Das stimmt nur, wenn in jedem Block die Stimuli alle gleich
  % wahrscheinlich sind dranzukommen. Das stimmt nicht immer. Insbesondere hatten wir oft Designs wo einige Stimuli in mancher Condition gar nicht dran kamen.
  disp('Optimal criterion values: ')
  disp(cval)
  disp('Optimal p(S2) values: ')
  disp(pS2)
else
  cval = nan(1,size(rprobs,1));
  pS2  = nan(1,size(rprobs,1));
end

% finally, plot the fit on top of the data
if plotting
    figure(h0),hold on
    [h2,fitted_data] = plot_fit(sequence,A,b,w,block,'r');
    [h3,simulated_data] = plot_data(sim_sequence,simulated_response,block,'g');
    % legend([h1, h2],'data','fit','Location','SouthOutside')
    legend boxoff
    title(lower(model))

    % added by Maik: plot optimal bias for each condition
%     if ~isempty(rprobs)
%       for i = 1:max(condition)
%         blocksInCondition = unique(block(condition==i));  % conditions may be repeated, e.g. baseline, therefore plot optimum for each block
%         for j = 1:numel(blocksInCondition)
%           h3 = plot([blocksInCondition(j)-1,blocksInCondition(j)],[pS2(i) pS2(i)],'k-');
%         end
%       end
%       legend([h1,h2,h3],'data','fit','sim','Location','best')
%     end
    legend([h1,h2,h3],'data','fit','sim','Location','best')


end

% =========================================================================
% Helper functions
% =========================================================================

% -------------------------------------------------------------------------
function [nll,w,A] = kdb_fit(param,stimulus,sequence,response,...
                             rewardarray,punisharray,A,b,model,integrator,...
                             stimulus_delta_mapping)

% Generate additional predictor variables that depend on the temporal
% integration parameter param
B = generate_predictors(stimulus,sequence,response,rewardarray,punisharray, ...
                        model,integrator,param,stimulus_delta_mapping);
A = [A B];

% Fit model using sdtcatfit.m
r1 = response==1;
r2 = response==2;
[w,nll] = sdtfit(A,b,r1,r2);%,1);

% -------------------------------------------------------------------------
function next = moving_average(last,x,i,lag)
% incrementally calculate a moving average, average at time point i (next)
% is returned, based on the last value at time point i-1 (last), a history
% vector x and a lag value (the size of the moving average window).
next=last;
fn = fieldnames(x);
for k = 1:length(fn)
    K = fn{k};
    j = i-lag;
    if j>0
        next.(K) = last.(K)-x(i-lag).(K)/lag;
    else
        next.(K) = last.(K);
    end
    next.(K) = next.(K)+x(i).(K)/lag;
end

% -------------------------------------------------------------------------
function next = leaky_integrator(last,x,i,param)
next=last;
fn = fieldnames(x);
for k = 1:length(fn)
    K = fn{k};
    next.(K) = param * last.(K) + x(i).(K);
end

% -------------------------------------------------------------------------
function next = cumulative_sum(last,x,i)
next=last;
fn = fieldnames(x);
for k = 1:length(fn)
    K = fn{k};
    next.(K) = last.(K) + x(i).(K);
end

% -------------------------------------------------------------------------
function x = trial_stats(stimulus,sequence,response,reward,punishment,num_stimuli)
x = struct([]);
for i = 1:length(stimulus)

    x(i).i1 = zeros(1,num_stimuli);
    x(i).i2 = zeros(1,num_stimuli);
    x(i).j1 = zeros(1,num_stimuli);
    x(i).j2 = zeros(1,num_stimuli);
    if sequence(i)  % false when initializing trial_stats with all 0s
        x(i).i1(sequence(i)) = (response(i)==1) * reward(i);
        x(i).i2(sequence(i)) = (response(i)==2) * reward(i);
        x(i).j1(sequence(i)) = (response(i)==1) * punishment(i);
        x(i).j2(sequence(i)) = (response(i)==2) * punishment(i);
    end
end

% -------------------------------------------------------------------------

function [i1_selected, i2_selected, j1_selected, j2_selected] = select_trials(i1, i2, j1, j2, idx)

i1_selected = i1(idx,:);
i2_selected = i2(idx,:);
j1_selected = j1(idx,:);
j2_selected = j2(idx,:);

% -------------------------------------------------------------------------

function array = rebuild_one(array_selected, idx)

array = nan(length(idx), size(array_selected,2));
array(idx,:) = array_selected;
% first row should contain zeros if there is no update
if idx(1) == 0
    array(1,:) = zeros(1, size(array,2));
end
% fill all NaN-rows with the values from the row above
while any(isnan(array),'all')
    nan_idx = isnan(array);
    previous = logical([nan_idx(2:end,:); zeros(1, size(array,2))]);
    array(nan_idx) = array(previous);
end

% -------------------------------------------------------------------------

function [i1, i2, j1, j2] = rebuild(i1_selected, i2_selected, j1_selected, j2_selected, idx)

i1 = rebuild_one(i1_selected, idx);
i2 = rebuild_one(i2_selected, idx);
j1 = rebuild_one(j1_selected, idx);
j2 = rebuild_one(j2_selected, idx);


% -------------------------------------------------------------------------

function pred = generate_predictors(stimulus,sequence,response,...
                                    rewardarray,punisharray,...
                                    model,integrator,param,...
                                    stimulus_delta_mapping)

% For this part of the design matrix, the exact stimulus identity does not
% matter, only which learning rate the stimulus is assigned to. So we
% replace sequence accordingly.
sequenceset = unique(sequence);
sequence_mod = sequence;
for s = sequenceset'
    sequence_mod(sequence_mod==s) = stimulus_delta_mapping(s);
end
sequenceset_mod = unique(sequence_mod);

% The learning part of the model depends on the following variables:
% i1: response 1 and reward
% i2: response 2 and reward
% j1: response 1 and punishment
% j2: response 2 and punishment

% this part should also work with rewardarray being the actual rewards
% and/or punisharray being the actual punishments rather than reward
% (punishment) in case the response was correct (wrong)

[rewards, punishments] = get_rewards_and_punishments(stimulus,response,rewardarray,punisharray);

i1 = [];
i2 = [];
for i = 1:length(sequenceset_mod)
    i1 = [i1, ((response==1) & (sequence_mod==i)) .* rewards];
    i2 = [i2, ((response==2) & (sequence_mod==i)) .* rewards];
end

j1 = [];
j2 = [];
for i = 1:length(sequenceset_mod)
    j1 = [j1, ((response==1) & (sequence_mod==i)) .* punishments];
    j2 = [j2, ((response==2) & (sequence_mod==i)) .* punishments];
end

% select trials where criterion should be updated
% If no_update_no_pullback is true, only the trials where a learning step
% happens are selected. Otherwise, all trials are selected.
% The unselected trials are ignored when computing the next criterion, and
% will later be re-inserted with the criterion value of the preceding trial.
if strcmp(integrator, 'leaky integrator on update')
    idx_punished = any(j1,2) | any(j2,2);
    idx_rewarded = any(i1,2) | any(i2,2);
    if ismember(model, {'model 1', 'model 1a', 'model 1b'})
        idx_update = idx_punished;
    elseif ismember(model, {'model 3', 'model 3a', 'model 3b'})
        idx_update = idx_rewarded;
    else
        idx_update = idx_punished | idx_rewarded;
    end
    i1 = i1(idx_update,:);
    i2 = i2(idx_update,:);
    j1 = j1(idx_update,:);
    j2 = j2(idx_update,:);
end

% integrate values
if strcmpi(integrator,'cumulative sum')
    % this will give the standard kdb model
    i1 = cumsum(i1);
    i2 = cumsum(i2);
    j1 = cumsum(j1);
    j2 = cumsum(j2);
else
    switch lower(integrator)
        case 'moving average'
            F = ones(param,1)/param;
        case {'leaky integrator', 'leaky integrator on update'}
            n = size(i1,1);
            F = cumprod([1; ones(n-1,1)*param]);
        otherwise
            error('unknown integrator')
    end
    i1 = filter(F,1,i1);
    i2 = filter(F,1,i2);
    j1 = filter(F,1,j1);
    j2 = filter(F,1,j2);
end

% re-insert trials where the criterion is not updated
% the criterion gets the same value as in the preceding trial
if strcmp(integrator, 'leaky integrator on update')
    [i1, i2, j1, j2] = rebuild(i1, i2, j1, j2, idx_update);
end


% shift by one so that the reward from trial i only enters the prediction
% of the i+1 trial
y.i1 = [zeros(1,length(sequenceset_mod)); i1(1:end-1,:)];
y.i2 = [zeros(1,length(sequenceset_mod)); i2(1:end-1,:)];
y.j1 = [zeros(1,length(sequenceset_mod)); j1(1:end-1,:)];
y.j2 = [zeros(1,length(sequenceset_mod)); j2(1:end-1,:)];
pred = predictors(y,model);

function pred = predictors(y,model)
i1 = cat(1,y.i1);  % why do we need cat? seems to do the same with i1=y.i1
i2 = cat(1,y.i2);
j1 = cat(1,y.j1);
j2 = cat(1,y.j2);
% depending on the model that is used we need different combinations of
% i1, i2, j1, and j2
switch lower(model)
    case 'general model'
        % +i1: if response is 1 and I get reinforced I want to shift to the
        % positive side and hence the learning rate should be positive
        % -i2: if response is 2 and I get reinforced I want to shift the
        % criterion to the negative side and hence the learning rate should
        % be postive, too
        % -j1: if response is 1 and I get punished I want to shift
        % to the negative side and hence the learning rate should be
        % positive, too
        % +j2: if response is 2 and I get punished I want to shift
        % to the positive side and hence the learning rate should be
        % positive, too
        pred = [+sum(i1,2), -sum(i2,2), -sum(j1,2), +sum(j2,2)];
    case 'model 1'
        % constrain the general model so that we only learn on
        % punished trials
        pred = [-sum(j1,2), +sum(j2,2)];
    case 'model 1a'
        % same as model 1 but force same learning rate for left and right
        % responses
        pred = -sum(j1,2)+sum(j2,2);
    case 'model 1b'
        % same as model 1a but allow a different learning rate for each
        % stimulus
        pred = -j1+j2;
    case 'model 2'
        % make sure that reinforced trials have same learning rate and that
        % punished trials have the same learning rate... but reinforced and
        % punished can be different
        pred = [(+sum(i1,2)-sum(i2,2)), (-sum(j1,2)+sum(j2,2))];
    case 'model 2b'
        % same as model 2 but allow a different learning rate for each
        % stimulus
        pred = [+i1-i2, -j1+j2];
    case 'model 3'
        % only learn on reinforced trials, in the original kdb model this
        % does not make sense because this will lead to exclusive
        % choice---but we have a forgetting term that pulls the criterion
        % back to zero
        pred = [+sum(i1,2), -sum(i2,2)];
    case 'model 3a'
        % same as model 3 but constrain learning rate to be the same for
        % one and two responses
        pred = +sum(i1,2)-sum(i2,2);
    case 'model 3b'
        % same as model 3a but allow a different learning rate for each
        % stimulus
        pred = +i1-i2;
    otherwise
        error('Unknown model restriction')
end


% -------------------------------------------------------------------------
function [rewards,punishments]=get_rewards_and_punishments(stimulus,...
    response,rewardarray,punisharray)

if isnan(rewardarray(1))  % 'treat_unpunished_as_rewarded'
    punishments = (stimulus~=response) .* punisharray;
    rewards = ~punishments;
elseif isnan(punisharray(1))  % 'treat_unrewarded_as_punished'
    rewards = (stimulus==response) .* rewardarray;
    punishments = ~rewards;
else  % 'standard
    rewards = (stimulus==response) .* rewardarray;
    punishments = (stimulus~=response) .* punisharray;
end


% -------------------------------------------------------------------------
function [stimulus,sequence,rewardarray,punisharray,original_trial]=...
    shuffle_data_conditionwise(stimulus,sequence,rewardarray,punisharray,...
    condition)

% This piece of codes generates the new conditionwise shuffled
% data. Maybe put this as another helper function instead of
% being here?

condition_start = [1; find(diff(condition))+1]; % index of first trial in each condition
condition_end = [find(diff(condition)); length(condition)]; % index of last trial in each condition
original_trial = zeros(size(condition));

for iCond=1:length(condition_start)
    This_Condition_Trial_idx          = condition_start(iCond):condition_end(iCond);            % Input Data
    Shuffled_idx_order                = randperm(size(This_Condition_Trial_idx, 2));            % Shuffle order
    This_Condition_Shuffled_Trial_idx = This_Condition_Trial_idx(Shuffled_idx_order);           % Output shuffled indices
    
    % shuffle variables according to the shuffled indices
    sequence(This_Condition_Trial_idx) = sequence(This_Condition_Shuffled_Trial_idx);
    stimulus(This_Condition_Trial_idx) = stimulus(This_Condition_Shuffled_Trial_idx);
    rewardarray(This_Condition_Trial_idx) = rewardarray(This_Condition_Shuffled_Trial_idx);
    punisharray(This_Condition_Trial_idx) = punisharray(This_Condition_Shuffled_Trial_idx);
    original_trial(This_Condition_Trial_idx) = This_Condition_Shuffled_Trial_idx;
end


% -------------------------------------------------------------------------
function [response,p,pred]=simulate_data(stimulus,sequence,rewardarray,...
    punisharray,w,model,integrator,param,StimulusMeans,stimulus_delta_mapping)

num_stimuli = numel(unique(sequence));

p = zeros(size(stimulus));          % probability of responding 1 on trial
c = zeros(size(stimulus));          % criterion
response = zeros(size(stimulus));   % the response (1 or 2)

if ~StimulusMeans
    mu = w(1:num_stimuli);
    delta = w(num_stimuli+1:end);
else
    mu = StimulusMeans;
    delta = w(2:end);
    c(1) = w(1);
end

% generate response for first trial (assuming criterion of 0)
i = 1;
p(i) = normcdf(-mu(sequence(i))+c(i));  % only here stimulus identity matters
response(i) = (rand>p(i))+1;

% For the rest of the simuluation, the exact stimulus identity does not
% matter, only which learning rate the stimulus is assigned to. So we
% replace sequence accordingly.
sequenceset = unique(sequence);
sequence_mod = sequence;
for s = sequenceset'
    sequence_mod(sequence_mod==s) = stimulus_delta_mapping(s);
end
num_stimuli_mod = numel(unique(sequence_mod));

[reward, punishment] = get_rewards_and_punishments(stimulus(i),response(i),rewardarray(i),punisharray(i));

% collect all the trial statistics
x(i) = trial_stats(stimulus(i),sequence_mod(i),response(i),reward,punishment,num_stimuli_mod);
x(length(stimulus))=trial_stats(0,0,0,0,0,num_stimuli_mod);
% integrated/averaged versions of all the trial stats
y(i) = trial_stats(0,0,0,0,0,num_stimuli_mod);
y(length(stimulus)) = trial_stats(0,0,0,0,0,num_stimuli_mod);

% now go through all the remaining trials
for i=2:length(stimulus)

    % integrate quantities from previous trials
    switch integrator
        case 'moving average'
            y(i) = moving_average(y(i-1),x,i-1,param);
        case 'leaky integrator'
            y(i) = leaky_integrator(y(i-1),x,i-1,param);
        case 'leaky integrator on update'
            punished = any(x(i-1).j1) || any(x(i-1).j2);
            rewarded = any(x(i-1).i1) || any(x(i-1).i2);
            if ismember(model, {'model 1', 'model 1a', 'model 1b'})
                update = punished;
            elseif ismember(model, {'model 3', 'model 3a', 'model 3b'})
                update = rewarded;
            else
                update = punished || idx_rewarded;
            end
            if update
                y(i) = leaky_integrator(y(i-1),x,i-1,param);
            else
                y(i) = y(i-1);
            end
        case 'cumulative sum'
            y(i) = cumulative_sum(y(i-1),x,i-1);
        otherwise
            error('unknown integrator')
    end

    % this is the decisive line: How to calculate/update the criterion?
    c(i) = predictors(y(i),model)*delta;

    % now generate a response for this trial
    p(i) = normcdf(-mu(sequence(i))+c(i));
    response(i) = (rand>p(i))+1;

    % and update the stats that we're interested in
    [reward, punishment] = get_rewards_and_punishments(stimulus(i),response(i),rewardarray(i),punisharray(i));
    x(i) = trial_stats(stimulus(i),sequence_mod(i),response(i),reward,punishment,num_stimuli_mod);

end
pred = predictors(y,model);

% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
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
plot((1:length(q))'-0.5,q,[color, '.'])
handle=plot((1:length(q))'-0.5,q,color);
raw_data = q;

% -------------------------------------------------------------------------
function [handle,fitted_data]=plot_fit(sequence,A,b,w,block,color,whichStim)
if ~exist('whichStim','var') || isempty(whichStim)
    whichStim = unique(sequence);
end
q = NaN * ones(max(block),1);
for i=1:(max(block)) % discard last block in plot as wished by Maik! Used to be max(block)-1
    select = block==i & ismember(sequence,whichStim);
    q(i) = mean(normcdf(A(select,:)*w+b(select)));
end
handle=plot((1:length(q))'-0.5,1-q,color);
fitted_data = 1-q;

% -------------------------------------------------------------------------
function [vORF] = getvORF(crit,mu,rprobs)
if size(mu,1)>size(mu,2),mu=mu';end
if size(rprobs,1)>size(rprobs,2),rprobs=rprobs';end
pRew = normcdf(crit,mu,1);
pRew(:,numel(mu)/2+1:end) = 1-pRew(:,numel(mu)/2+1:end);
pRew = pRew.*rprobs;
vORF = 1-mean(pRew);

% -------------------------------------------------------------------------
function [rewardarray, reward_density] = scale_reward_with_density(stimulus, response, rewardarray, condition, received_or_programmed)
if ~exist('received_or_programmed','var')
    received_or_programmed = 'received';
end
received_reward = (stimulus==response) & rewardarray>0;
start_idx = [1; find(diff(condition))+1; length(condition)+1];
reward_density = zeros(length(start_idx)-1,1);
for i = 1:length(start_idx)-1
    if strcmp(received_or_programmed, 'programmed')
        condition_reward = rewardarray(start_idx(i):start_idx(i+1)-1);
    elseif strcmp(received_or_programmed, 'received')
        condition_reward = received_reward(start_idx(i):start_idx(i+1)-1);
    end
    reward_density(i) = sum(condition_reward) / length(condition_reward);
    rewardarray(start_idx(i):start_idx(i+1)-1) = rewardarray(start_idx(i):start_idx(i+1)-1) / reward_density(i);
end
