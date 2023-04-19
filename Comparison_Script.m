%% Script to generate results presented in "Exact solutions for diffusive
%% transport on heterogeneous growing domains" by Johnston and Simpson, 2023.
%% This script will provide results for either a single domain or three
%% domains, symmetric around zero, as in the manuscript. Code is set up to
%% generate the results in Figure 5(c).

close all
clear all

%% Define options for model simulations

growthType = 'linear';                      % Domain growth options: either 'linear' 'exponential' or 'oscillating'.
includeThreeTimes = 'yes';                  % Option that will present results at 3 different time points: 'no' or 'yes'.
onlySplitting = 'no';                       % Option that will only return splitting/survival probabilities at the final time and not.

%% Define parameters used in the model

L0 = [50; 50];                              % Vector of initial domain sizes.
beta = [0.02; 0.01];                        % Growth rate for each domain. Note that beta values for exponential growth are negative, positive for exponential contraction.
tEnd = 5000;                                % Final simulation time.
x0 = [0];                                   % Initial location of individuals. Can be a vector.
Pm = [1; 0.5];                              % Probability of movement for each domain.
D = Pm/2;                                   % Diffusivity for each domain.
minL = [30; 30];                            % Minimum domain size (only used for exponential with beta < 0 and oscillating).
nIndividuals = 1000;                        % Number of agents in an individual realisation of the random walk.
nRepeats = 1000;                            % Number of repeats of the random walk.
nSimulationSavePoints = min(501,tEnd);      % Number of time points to save results.
nTerms = 10;                                % Number of terms in the infinite summations evaluated.
nMOISavePoints = 300;                       % Number of time points in the PDE saved.
nMOISpacePoints = 1001;                     % Number of space points in the PDE in each domain.
startingPoints = x0;                        % Location of Dirac delta initial condition.
nDomains = numel(L0);                       % Number of domains.
nStartingPoints = numel(startingPoints);    % Number of starting points.
T1Time = ceil(tEnd/5);                      % Time to plot first solution.
T2Time = ceil(tEnd/2);                      % Time to plot second solution.

%% Initialise summary statistics

finalRightSplit = zeros(nStartingPoints,1); % Final right splitting probability.
finalLeftSplit = zeros(nStartingPoints,1);  % Final left splitting probability.
finalSurvive = zeros(nStartingPoints,1);    % Final survival probability.

relativeD = sqrt(flipud(D))/sum(sqrt(D));   % Calculate weights from diffusivity.

defineGrowthFunction;                       % Define growth function and associated transforms.

if strcmpi(onlySplitting,'no')
    MultidomainRandomWalk;                  % Perform realisations of the random walk.
end
if numel(L0) == 2
    TwoDomainExactSolution;                 % Enumerate exact solution for multiple domains.
else
    OneDomainExactSolution;                 % Enumerate exact solution for a single domain.
end

PlotResults;                                % Plot results.

% Generate solution profiles for multiple times
if strcmpi(includeThreeTimes,'yes')
    
    locationsEnd_Store = locationsEnd;      % Store old final locations.
    locationsEnd = locationsT2Time;         % New "final" locations.    
    tEnd_Store = tEnd;                      % Store old final time.
    tEnd = T2Time;                          % New "final" time.
    
    % Enumerate exact solutions.
    if numel(L0) == 2
        TwoDomainExactSolution;             % Enumerate exact solution for multiple domains.
    else
        OneDomainExactSolution;             % Enumerate exact solution for a single domain.
    end
    
    PlotResults                             % Plot results.
    
    locationsEnd = locationsT1Time;         % New "final" locations.
    tEnd = T1Time;                          % New "final" time.
    
    % Enumerate exact solutions.
    if numel(L0) == 2
        TwoDomainExactSolution;             % Enumerate exact solution for multiple domains.
    else
        OneDomainExactSolution;             % Enumerate exact solution for a single domain.
    end
    
    PlotResults                             % Plot results.
    
end

% Store final splitting and survival probabilites.
finalRightSplit = rightSplit(:,end);
finalLeftSplit = leftSplit(:,end);
finalSurvive = 1 - leftSplit(:,end) - rightSplit(:,end);