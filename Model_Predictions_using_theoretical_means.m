%% Model_Predictions_using_theoretical_means
% This function plots the theoretical criteria in the steady state of the
% different condition asuming the Stimulus Means = [-1.5, -0.5, 0, 0.5,
% 1.5].It does so for the initially comnsidered SDT learning models, namely:
% optimal, income (in the paper IR), unrewarded- (IRO in the paper), error-
% and  income & unrewarded-based models (IR&RO) .  These plots are 
% are responsible for the Figures 1d, Figures S1.1 and S1.2.
% 
%                                                       Luis 09 / 02 / 2024

clc,grey=ones(1,3)*0.7;
figure('units','normalized','position',[.1,.1,.8,.8])

d     = 1.5;    % multiplication factor
gamma = 0.99;   % whether the design works is very sensitive to gamma
delta = 0.04;
conditions = [1,2,3,4,5,6,7]; %plotted from top (1) to bottom (7)

for j = 1:length(conditions)

  switch conditions(j)
            case 1 %BT  
                mu = [-1 0 0 1] * d;
                p = 1/4*[1 1 1 1];
                class = [1 1 2 2];
                r = [1.0 1.0 1.0 1.0];
            case 2 %CL
                mu = [-1 0 1/3 1] * d;
                p = 1/4*[1 1 1 1];
                class = [1 2 1 2];
                r = [1.0 1.0 1.0 1.0];

            case 3 %CR
                mu = [-1 -1/3 0 1] * d;
                p = 1/4*[1 1 1 1];
                class = [1 2 1 2];
                r = [1.0 1.0 1.0 1.0];

            case 4 %RPL (Lean L)
                mu = [-1 0 0 1] * d;
                p = 1/4*[1 1 1 1];
                class = [1 1 2 2];
                r = [1.0 0.0 0.5 0.5];

            case 5 %RPR (Lean R)
                mu = [-1 0 0 1] * d;
                p = 1/4*[1 1 1 1];
                class = [1 1 2 2];
                r = [0.5 0.5 0 1];

            case 6 %SOL (Rich L)
                mu = [-1 0 1] * d;
                p = 1/4*[2 1 1];
                class = [1 2 2];
                r = [1.0 1.0 1.0];

            case 7 %SOR (Rich R)
                mu = [-1 0 1] * d;
                p = 1/4*[1 1 2];
                class = [1 1 2];
                r = [1.0 1.0 1.0];

  end

  c = linspace(-2.5+min(mu),max(mu)+2.5,1000);
  [Rf1,Rf2] = reinforcement(c,mu,p,class,r);    % Rf1 and Rf2 are expected reinforcers per trial as a function of criterion values, FOR EACH CATEGORY
  [Er1,Er2] = errors(c,mu,p,class,r);           % Er1 and Er2 are expected errors per trial, FOR EACH CATEGORY
  [Ur1,Ur2] = unrewarded(c,mu,p,class,r);       % Ur1 and Ur2 are expected nonrewards per trial, FOR EACH CATEGORY

  % optimal model
  subplot(length(conditions),4,(j-1)*4+1),hold on
  plot_densities(c,mu,p,class,r);
  plot(c,Rf1+Rf2,'g');
  set(gca,'ylim',[0,0.85],'xlim',[min(c),max(c)]);
  if j == 1; title('Optimal'); end
  ylabel('E(Rf_1)+E(Rf_2)')
  if j == max(conditions); xlabel('criterion'); end
  c_optimal = fminbnd(@(c) f_optimal(c,mu,p,class,r),min(c),max(c));
  plot([c_optimal,c_optimal],ylim,'g')

  % income-based model (IR)
  subplot(length(conditions),4,(j-1)*4+2),hold on
  plot(c,Rf1,'k')
  plot(c,-Rf2,'Color',grey)
  plot(c,Rf1-Rf2,'b');
  plot(c,-(gamma*c-c)/delta,'b:')
  set(gca,'ylim',[-0.75,+0.75],'xlim',[min(c),max(c)],'ytick',[-3/4,-1/2,-1/4,0,1/4,1/2,3/4]);
  if j == 1; title('Rewarded-based'); end
  ylabel('E(Rf_1)-E(Rf_2)')
  if j == max(conditions); xlabel('criterion'); end
  c_income = fzero(@(c) f_income(c,mu,p,class,r,gamma,delta),[min(c),max(c)]);
  plot([c_income,c_income],ylim,'b')

%   % error-based model
%   subplot(length(conditions),4,(j-1)*4+3),hold on
%   plot(c,-Er1,'k')
%   plot(c,Er2,'Color',grey)
%   plot(c,Er2-Er1,'b');
%   plot(c,-(gamma*c-c)/delta,'b:')
%   set(gca,'ylim',[-0.75,+0.75],'xlim',[min(c),max(c)],'ytick',[-3/4,-1/2,-1/4,0,1/4,1/2,3/4]);if j == 1; title('Error-based criterion shift'); end
%   ylabel('E(Er_2)-E(Er_1)')
%   if j == max(conditions); xlabel('criterion'); end
%   c_errors = fzero(@(c) f_errors(c,mu,p,class,r,gamma,delta),[min(c),max(c)]);
%   plot([c_errors,c_errors],ylim,'b')

  % unrewarded-based model (IRO)
  subplot(length(conditions),4,(j-1)*4+3),hold on
  plot(c,-Ur1,'k')
  plot(c,Ur2,'Color',grey)
  plot(c,Ur2-Ur1,'r');
  plot(c,-(gamma*c-c)/delta,'r:')
  set(gca,'ylim',[-1,+1],'xlim',[min(c),max(c)],'ytick',[-1,-3/4,-1/2,-1/4,0,1/4,1/2,3/4,1]);
  if j == 1; title('Unrewarded-based'); end
  ylabel('E(noRf_2)-E(noRf_1)')
  if j == max(conditions); xlabel('criterion');end
  c_unrew = fzero(@(c) f_unrewarded(c,mu,p,class,r,gamma,delta),[min(c),max(c)]);
  plot([c_unrew,c_unrew],ylim,'r')

    % model that learns after all trials (IR&IRO)
  subplot(length(conditions),4,(j-1)*4+4),hold on
  plot(c,Ur2-Ur1,'k')
  plot(c,Rf1-Rf2,'Color',grey)
  plot(c,Rf1-Rf2+Ur2-Ur1,'b');
  plot(c,-(gamma*c-c)/delta,'b:')
  set(gca,'ylim',[-1,+1],'xlim',[min(c),max(c)],'ytick',[-1,-3/4,-1/2,-1/4,0,1/4,1/2,3/4,1]);
  if j == 1; title('All-trial-based criterion shift'); end
  ylabel('E(Rf_1)-E(Rf_2)+E(noRf_2)-E(noRf_1)')
  if j == max(conditions); xlabel('criterion');end
  c_all = fzero(@(c) f_all(c,mu,p,class,r,gamma,delta),[min(c),max(c)]);
  plot([c_all,c_all],ylim,'b')

  if j == 6
    Ur2(1)
    Ur1(end)
  end

end


%% function y = f_optimal(c,mu,p,class,r)
function y = f_optimal(c,mu,p,class,r)
[Rf1,Rf2] = reinforcement(c,mu,p,class,r);
y = -Rf1-Rf2;
end

%% function y = f_income(c,mu,p,class,r,gamma,delta)
function y = f_income(c,mu,p,class,r,gamma,delta)
[Rf1,Rf2] = reinforcement(c,mu,p,class,r);
y = (Rf1-Rf2+(gamma*c-c)/delta);
end

%% function y = f_errors(c,mu,p,class,r,gamma,delta)
function y = f_errors(c,mu,p,class,r,gamma,delta)
[Er1,Er2] = errors(c,mu,p,class,r);
y = (Er2-Er1+(gamma*c-c)/delta);
end

%% function y = f_unrewarded(c,mu,p,class,r,gamma,delta)
function y = f_unrewarded(c,mu,p,class,r,gamma,delta)
[Ur1,Ur2] = unrewarded(c,mu,p,class,r);
y = (Ur2-Ur1+(gamma*c-c)/delta);
end

%% function y = f_all(c,mu,p,class,r,gamma,delta)
function y = f_all(c,mu,p,class,r,gamma,delta)
[Rf1,Rf2] = reinforcement(c,mu,p,class,r);
[Ur1,Ur2] = unrewarded(c,mu,p,class,r);
y = Rf1-Rf2+Ur2-Ur1+(gamma*c-c)/delta;
end

%% function plot_densities(c,mu,p,class,r)
function plot_densities(c,mu,p,class,r)
grey = ones(1,3)*0.7;
%s = length(class)/2;
s = 1;
y = zeros(size(c));
for i = find(class==1)
  %plot(c,s*p(i)*r(i)*normpdf(c,mu(i),1),'k:') % plot each class 1 distribution (dotted lines)
  y = y+s*p(i)*r(i)*normpdf(c,mu(i),1);
end
plot(c,y,'k') % plot sum of all class 1 distributions
y = zeros(size(c));
for i = find(class==2)
  %plot(c,s*p(i)*r(i)*normpdf(c,mu(i),1),':','Color',grey)
  y = y+s*p(i)*r(i)*normpdf(c,mu(i),1);
end
plot(c,y,'Color',grey)
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

%% function [Er1,Er2] = errors(c,mu,p,class,r)
function [Er1,Er2] = errors(c,mu,p,class,r)
% This function is different from that in new_experiment4.m of the same name.
Er2 = 0;
for i=find(class==1)
  Er2 = Er2 + p(i)*(1-normcdf(c,mu(i),1));
end

Er1 = 0;
for i = find(class==2)
  Er1 = Er1 + p(i)*normcdf(c,mu(i),1);
end
end

%% function [R1,R2] = response_probabilities(c,mu,p)
function [R1,R2] = response_probabilities(c,mu,p)
R1 = 0;
for i = 1:length(mu)
  R1 = R1 + p(i) * normcdf(c,mu(i),1);
end
R2 = 1-R1;
end

%% function [Ur1,Ur2] = unrewarded(c,mu,p,class,r)
function [Ur1,Ur2] = unrewarded(c,mu,p,class,r)
[Rf1,Rf2] = reinforcement(c,mu,p,class,r);
[R1,R2] = response_probabilities(c,mu,p);
Ur1 = R1-Rf1;
Ur2 = R2-Rf2;
end
