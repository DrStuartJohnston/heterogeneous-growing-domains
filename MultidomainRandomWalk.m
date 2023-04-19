%% Script to perform random walk simulations for "Exact solutions for diffusive
%% transport on heterogeneous growing domains" by Johnston and Simpson, 2023.

%% Initialise the random walk

latticeSizeStart = sum(L0);                                 % Overall number of initial lattice sites in the domain.
latticeSizesStart = L0;                                     % Number of initial lattice sites in each domain.

nLeft = zeros(nStartingPoints,nSimulationSavePoints);       % Number of agents that cross the left external boundary.
nRight = zeros(nStartingPoints,nSimulationSavePoints);      % Number of agents that cross the right external boundary.
nSurvive = zeros(nStartingPoints,nSimulationSavePoints);    % Number of agents that do not cross the boundaries.        
savePoints = linspace(0,tEnd,nSimulationSavePoints);        % Times at which to save simulation results.

locationsEnd = cell(nStartingPoints,1);                     % Cell containing the final locations of each agent.
for i = 1:nStartingPoints
    locationsEnd{i} = 1e6*ones(nRepeats*nIndividuals,1);    % Initialise locations.
end

locationsT1Time = cell(nStartingPoints,1);                  % Cell containing the locations of each agent at first solution time.
for i = 1:nStartingPoints
    locationsT1Time{i} = 1e6*ones(nRepeats*nIndividuals,1); % Initialise locations.
end

locationsT2Time = cell(nStartingPoints,1);                  % Cell containing the locations of each agent at second solution time.
for i = 1:nStartingPoints
    locationsT2Time{i} = 1e6*ones(nRepeats*nIndividuals,1); % Initialise locations.
end

% Main random walk loop. Loops over the vector of starting points to allow
% for ease of analysis of impact of initial locations.

for iStart = 1:nStartingPoints
    iStart                                                      % Print number of starting point
    for iRepeat = 1:nRepeats
        
        iRepeat                                                 % Print realisation number
        startingPoint = startingPoints(iStart);                 % Initial location of agents.
        counter = 1;                                            % Counting variable for save points.
        nLeftTmp = 0;                                           % Temporary variable for agents that cross the left boundary.
        nRightTmp = 0;                                          % Temporary variable for agents that cross the right boundary.
        latticeSize = latticeSizeStart;                         % Current number of total lattice site.
        latticeSizes = latticeSizesStart;                       % Current number of lattice sites in each domain.
        domainSizes = latticeSizes;                             % Current domain sizes.
        initialLocation = startingPoints(iStart);               % Initial location of agents.
        locations = initialLocation*ones(nIndividuals,1);       % Current location of each agent.
        nRemaining = nIndividuals;                              % Current number of agents remaining in the simulation.
        
        for iTime = 1:tEnd
            
            % Save relevant information if first solution time is exceeded.
            if iTime == T1Time+1
                locationsT1Time{iStart}((iRepeat-1)*nIndividuals+1:(iRepeat-1)*nIndividuals+nRemaining) = locations;
                latticeSizeT1 = latticeSize;
            end
            
            % Save relevant information if second solution time is
            % exceeded.
            if iTime == T2Time+1
                locationsT2Time{iStart}((iRepeat-1)*nIndividuals+1:(iRepeat-1)*nIndividuals+nRemaining) = locations;
                latticeSizeT2 = latticeSize;
            end
            
            % Save summary statistics as each save point is exceeded.
            if iTime >= savePoints(counter)
                nLeft(iStart,counter) = nLeft(iStart,counter) + nLeftTmp/nRepeats;
                nRight(iStart,counter) = nRight(iStart,counter) + nRightTmp/nRepeats;
                nSurvive(iStart,counter) = nSurvive(iStart,counter) + nRemaining/nRepeats;
                counter = counter + 1;
            end
            
            oldLatticeSize = latticeSize;                       % Temporary copy of number of total lattice sites.
            oldLatticeSizes = latticeSizes;                     % Temporary copy of number of lattice sites in each domain.
            
            % Update the domain size according to the domain growth
            % function.
            for iDomain = 1:nDomains
                domainSizes(iDomain) = domainGrowthFunction(iTime,iDomain);
            end
            
            % Update the number of lattice sites according to new domain
            % sizes.
            latticeSizes = floor(domainSizes);
            latticeSize = sum(latticeSizes);
            
            % Check if a new lattice site must be added due to domain
            % growth.
            if any(latticeSizes > oldLatticeSizes)
                
                addSizesPos = latticeSizes-oldLatticeSizes;                 % Number of sites to be added in the positive half-domain.
                addSizesNeg = latticeSizes-oldLatticeSizes;                 % Number of sites to be added in the negative half-domain.
                newSitesPos = ceil(oldLatticeSizes.*rand(numel(L0),1));     % Locations of new sites (within a domain) to be added in the positive half-domain.
                newSitesNeg = -ceil(oldLatticeSizes.*rand(numel(L0),1));    % Locations of new sites (within a domain) to be added in the negative half-domain.
                
                % Loop over each domain to add new sites.
                for iDomain = 1:numel(L0)
                    if latticeSizes(iDomain) > oldLatticeSizes(iDomain)
                        actualNewSitePos = newSitesPos(iDomain) + sum(latticeSizes(1:iDomain-1));               % Locations of new sites (overall) to be added in the positive half-domain.
                        actualNewSiteNeg = newSitesNeg(iDomain) - sum(latticeSizes(1:iDomain-1));               % Locations of new sites (overall) to be added in the negative half-domain.
                        locations(locations>=actualNewSitePos) = locations(locations>=actualNewSitePos) + ...
                            addSizesPos(iDomain);                                                               % Shift all agents to the right of the new sites.
                        locations(locations<=actualNewSiteNeg) = locations(locations<=actualNewSiteNeg) - ...
                            addSizesNeg(iDomain);                                                               % Shift all agents to the left of the new sites.
                    end
                end
            end
            
            % Check if a lattice site must be removed due to new domain
            % sizes.
            if any(latticeSizes < oldLatticeSizes)
                
                removeSizesPos = oldLatticeSizes-latticeSizes;              % Number of sites to be removed in the positive half-domain.
                removeSizesNeg = oldLatticeSizes-latticeSizes;              % Number of sites to be removed in the negative half-domain.
                deleteSitesPos = ceil(oldLatticeSizes.*rand(numel(L0),1));  % Locations of new sites (within a domain) to be removed in the positive half-domain.
                deleteSitesNeg = -ceil(oldLatticeSizes.*rand(numel(L0),1)); % Locations of new sites (within a domain) to be removed in the negative half-domain.
                
                % Loop over each domain to add new sites.
                for iDomain = 1:numel(L0)
                    if latticeSizes(iDomain) < oldLatticeSizes(iDomain)
                        actualDeleteSitePos = deleteSitesPos(iDomain) + sum(latticeSizes(1:iDomain-1));             % Locations of new sites (overall) to be removed in the positive half-domain.
                        actualDeleteSiteNeg = deleteSitesNeg(iDomain) - sum(latticeSizes(1:iDomain-1));             % Locations of new sites (overall) to be removed in the negative half-domain.
                        locations(locations>actualDeleteSitePos) = locations(locations>actualDeleteSitePos) - ...   
                            removeSizesPos(iDomain);                                                                % Shift all agents to the right of the removed sites.   
                        locations(locations<actualDeleteSiteNeg) = locations(locations<actualDeleteSiteNeg) + ...   
                            removeSizesNeg(iDomain);                                                                % Shift all agents to the left of the removed sites.
                    end
                end
            end
            
            % Check if all agents have left the domain and terminate after
            % saving the relevant date if so.
            if nRemaining == 0
                nRight(iStart,counter:end) = nRight(iStart,counter:end) + nRightTmp/nRepeats;
                nLeft(iStart,counter:end) = nLeft(iStart,counter:end) + nLeftTmp/nRepeats;
                nSurvive(iStart,counter:end) = nSurvive(iStart,counter:end) + nRemaining/nRepeats;
                break
            end
            
            latticePoints = [0;cumsum(latticeSizes)];   % Location of internal boundaries.
            
            % Select number of individuals for movement.
            for iInd = 1:nRemaining
                indiv = ceil(rand*nRemaining);          % Number of individual.   
                direction = rand;                       % Random number for movement direction.        
                moveCheck = 0;
                absLoc = abs(locations(indiv));         % Location of individual for determining probability of movement.
               
                % Determine probability of movement based on relevant
                % domain.
                if locations(indiv) >= 0
                    for i = 1:numel(L0)
                        if absLoc < latticePoints(i+1) && ...
                                absLoc >= latticePoints(i)
                            moveCheck = Pm(i);          % Set probability of movement.
                        end
                    end
                else
                    for i = 1:numel(L0)
                        if absLoc <= latticePoints(i+1) && ...
                                absLoc > latticePoints(i)
                            moveCheck = Pm(i);          % Set probability of movement.
                        end
                    end
                end
                
                % Perform movement event.
                if direction < moveCheck/2
                    % Left movement event.
                    locations(indiv) = locations(indiv) - 1;
                    % Check if individual crosses left external boundary.
                    if locations(indiv) <= -latticeSize
                        locations(indiv) = [];
                        nLeftTmp = nLeftTmp + 1;
                        nRemaining = nRemaining-1;
                    end
                elseif direction > moveCheck/2 && direction < moveCheck
                    % Right movement event.
                    locations(indiv) = locations(indiv) + 1;
                    % Check if individual crosses right external boundary.
                    if locations(indiv) >= latticeSize
                        locations(indiv) = [];
                        nRightTmp = nRightTmp + 1;
                        nRemaining = nRemaining-1;
                    end
                end
            end
        end
        
        % Save final locations of individuals
        locationsEnd{iStart}((iRepeat-1)*nIndividuals+1:(iRepeat-1)*nIndividuals+nRemaining) = locations;
    end
    
    % Remove individuals that have crossed the external boundaries.
    locationsEnd{iStart}(locationsEnd{iStart}==1e6) = [];
    locationsT1Time{iStart}(locationsT1Time{iStart}==1e6) = [];
    locationsT2Time{iStart}(locationsT2Time{iStart}==1e6) = [];
end