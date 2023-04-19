%% Script to plot results for "Exact solutions for diffusive
%% transport on heterogeneous growing domains" by Johnston and Simpson, 2023.

figCount = 0;                                                               % Counting variable for figure.

% Plot non-splitting probabilty results if desired.
if strcmpi(onlySplitting,'no')
    % Loop over all possible starting points
    for iIC = 1:nStartingPoints
        
        figCount = figCount + 1;
        figure(figCount);
        
        % Calculate normalised counts of agent locations.
        counts = histcounts(locationsEnd{iIC},floor(allBoundaries(1)):floor(allBoundaries(end))+1, ...
            'Normalization','pdf')*numel(locationsEnd{1})/(nIndividuals*nRepeats);
        
        % Plot counts and solution profiles in each domain.
        for i = 1:2*nDomains
            plotPoints = floor(allBoundaries(i))+0.5:floor(allBoundaries(i+1))-0.5;                     % Calculate plot locations in each domain.
            plot(plotPoints,counts(plotPoints-floor(allBoundaries(1))+0.5),'linewidth',3); hold on;
            index = (i-1)*nMOISpacePoints+1:i*nMOISpacePoints;                                          % Indices for relevant solution.
            plot(currentValidX(index),C_finalStore(iIC,index),'linewidth',4)
        end
        
        % Plotting splitting and survival probabilities.
        figCount = figCount + 1;
        figure(figCount); hold on; plot(savePoints,nLeft(iIC,:)/nIndividuals,'linewidth',3); hold on
        plot(savePoints,nRight(iIC,:)/nIndividuals,'linewidth',3);
        plot(savePoints,nSurvive(iIC,:)/nIndividuals,'linewidth',3);
        plot(timeValues,1-leftSplit(iIC,:)-rightSplit(iIC,:),'--','linewidth',3);
        plot(timeValues,rightSplit(iIC,:),'--','linewidth',3);
        plot(timeValues,leftSplit(iIC,:),'--','linewidth',3);
    end
end

% Plot final splitting and survival probabilities if more than one starting
% location is considered.
if nStartingPoints > 1
    % Plot both simulation and exact solutions.
    if strcmpi(onlySplitting,'no')
        figure; hold on;
        plot(x0,nSurvive(:,end)/nIndividuals,'linewidth',3)
        plot(x0,nLeft(:,end)/nIndividuals,'linewidth',3)
        plot(x0,nRight(:,end)/nIndividuals,'linewidth',3)
        plot(x0,1-leftSplit(:,end)-rightSplit(:,end),'--','linewidth',3)
        plot(x0,rightSplit(:,end),'--','linewidth',3)
        plot(x0,leftSplit(:,end),'--','linewidth',3)
    % Plot only exact solutions.
    else
        figure; hold on;
        plot(x0,1-leftSplit(:,end)-rightSplit(:,end),'--','linewidth',3)
        plot(x0,rightSplit(:,end),'--','linewidth',3)
        plot(x0,leftSplit(:,end),'--','linewidth',3)
    end
end