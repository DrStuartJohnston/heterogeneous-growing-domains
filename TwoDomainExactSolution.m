%% Script to enumerate exact multidomain solutions for "Exact solutions for diffusive
%% transport on heterogeneous growing domains" by Johnston and Simpson, 2023.

%% Initialise variables

totalMass = zeros(nStartingPoints,nMOISavePoints);                  % Total mass across all domains.
middleMass = zeros(nStartingPoints,nMOISavePoints);                 % Mass in the inner domain.
leftMass = zeros(nStartingPoints,nMOISavePoints);                   % Mass in the outer left domain.
rightMass = zeros(nStartingPoints,nMOISavePoints);                  % Mass in the outer right domain.
leftSplit = zeros(nStartingPoints,nMOISavePoints);                  % Left splitting probability.
rightSplit = zeros(nStartingPoints,nMOISavePoints);                 % Right splitting probability.
leftInternalSplit = zeros(nStartingPoints,nMOISavePoints);          % Left internal splitting probability.
rightInternalSplit = zeros(nStartingPoints,nMOISavePoints);         % Right internal splitting probability.
C_finalStore = zeros(nStartingPoints,4*nMOISpacePoints);            % Final solution profile.

% Loop over each starting location.
for iStart = 1:nStartingPoints
    
    C = zeros(nMOISavePoints,4*nMOISpacePoints);                    % All solution profiles for a single starting point.
    convolutionTerms = zeros(nMOISavePoints,4*nMOISpacePoints);     % All contributions from convolution terms.
    timeValues = linspace(tEnd/nMOISavePoints,tEnd,nMOISavePoints); % Time points where solution is enumerated.
    domainSizes = zeros(2,1);                                       % Current domain sizes.
    rightFlux = zeros(1,nMOISavePoints);                            % Flux through the internal right boundary.
    leftFlux = zeros(1,nMOISavePoints);                             % Flux through the internal left boundary.
    dx = zeros(4,nMOISavePoints);                                   % Current node spacing (space step).
    dt = timeValues(2)-timeValues(1);                               % Time step.
    
    % Loop over each time point.
    for iTime = 1:nMOISavePoints
        
        t = timeValues(iTime);                                                  % Current time point.
        
        % Calculate current domain sizes according to growth functions.
        for iDomain = 1:nDomains
            domainSizes(iDomain) = domainGrowthFunction(t,iDomain);
        end
        allBoundaries = [-flipud(cumsum(domainSizes));0;cumsum(domainSizes)];   % Current boundary locations.
        timeTransform = transformFunction(t,1);                                 % Perform time transformation.
        currentValidX = zeros(1,4*nMOISpacePoints);                             % Current transformed x values.    
        
        %% Solution for the -L_1(t) < x < 0 domain.
        
        xScale = L0(1)/domainSizes(1);                                          % Rescaling of x.
        Dn = D(1);                                                              % Relevant diffusivity.
        validX = linspace(allBoundaries(2),allBoundaries(3),nMOISpacePoints);   % Relevant x values.
        indices = nMOISpacePoints+1:2*nMOISpacePoints;                          % Indices for storage of solution.
        
        % Enumerate solution in the domain.
        C(iTime,indices) = 1/sqrt(4*pi*Dn*timeTransform)*xScale*exp(-(validX*xScale-x0(iStart)).^2/(4*Dn*timeTransform));
        for iTerm = 1:nTerms
            C(iTime,indices) = C(iTime,indices) + 1/sqrt(4*pi*Dn*timeTransform)*xScale*( ...
                (-1)^iTerm*exp(-(validX*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) + ...
                (-1)^iTerm*exp(-(validX*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)));
        end
        
        % Calculate the flux through the boundary at x = -L_1(t).
        leftFlux(iTime) = -2*Dn*(allBoundaries(2)*xScale-x0(iStart))/(sqrt(pi)*(4*Dn*timeTransform)^(3/2))*xScale^2* ...
            exp(-(allBoundaries(2)*xScale-x0(iStart)).^2/(4*Dn*timeTransform));
        for iTerm = 1:nTerms
            leftFlux(iTime) = leftFlux(iTime) - 1/(sqrt(pi)*(4*Dn*timeTransform)^(3/2))*xScale^2*( ...
                (-1)^iTerm*exp(-(allBoundaries(2)*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) * ...
                2*Dn*(allBoundaries(2)*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))) + ...
                (-1)^iTerm*exp(-(allBoundaries(2)*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) * ...
                2*Dn*(allBoundaries(2)*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))));
        end
        
        %% Solution for the 0 < x < L_1(t) Domain.
        
        xScale = L0(1)/domainSizes(1);                                          % Rescaling of x.
        Dn = D(1);                                                              % Relevant diffusivity.
        validX = linspace(allBoundaries(3),allBoundaries(4),nMOISpacePoints);   % Relevant x values.
        indices = 2*nMOISpacePoints+1:3*nMOISpacePoints;                        % Indices for storage of solution.
        
        % Enumerate solution in the domain.
        C(iTime,indices) = 1/sqrt(4*pi*Dn*timeTransform)*xScale*exp(-(validX*xScale-x0(iStart)).^2/(4*Dn*timeTransform));
        for iTerm = 1:nTerms
            C(iTime,indices) = C(iTime,indices) + 1/sqrt(4*pi*Dn*timeTransform)*xScale*( ...
                (-1)^iTerm*exp(-(validX*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) + ...
                (-1)^iTerm*exp(-(validX*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)));
        end
        
        % Calculate the flux through the internal boundary at x = L_1(t).
        rightFlux(iTime) = 2*Dn*(allBoundaries(4)*xScale-x0(iStart))/(sqrt(pi)*(4*Dn*timeTransform)^(3/2))*xScale^2* ...
            exp(-(allBoundaries(4)*xScale-x0(iStart)).^2/(4*Dn*timeTransform));
        for iTerm = 1:nTerms
            rightFlux(iTime) = rightFlux(iTime) + 1/(sqrt(pi)*(4*Dn*timeTransform)^(3/2))*xScale^2*( ...
                (-1)^iTerm*exp(-(allBoundaries(4)*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) * ...
                2*Dn*(allBoundaries(4)*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))) + ...
                (-1)^iTerm*exp(-(allBoundaries(4)*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) * ...
                2*Dn*(allBoundaries(4)*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))));
        end
        
        %% Introduction of the convolution terms.
        
        %% Solution for the -L_1(t) < x < 0 domain.
        
        Dn = D(1);                                                              % Relevant diffusivity.      
        validX = linspace(allBoundaries(2),allBoundaries(3),nMOISpacePoints);   % Relevant x values.
        dx(2,iTime) = validX(2)-validX(1);                                      % Current space step.
        indices = nMOISpacePoints+1:2*nMOISpacePoints;                          % Indices for storage of solution.
        currentValidX(indices) = validX;                                        % Relevant x values.
        
        % Do not calculate solution if only the final splitting
        % probabilities are wanted.
        if strcmpi(onlySplitting,'no')
            % Enumerate convolution terms by integrating over time.        
            for jTime = 1:iTime-1
                tmpDomainSizes = zeros(nDomains,1);                             % Temporary domain size.    
                tau = timeValues(jTime);                                        % Temporary time.
                tauTimeTransform = relativeTransformFunction(tau,t,1);          % Temporary time transform.
                
                % Calculate temporary domain size.
                for iDomain = 1:nDomains
                    tmpDomainSizes(iDomain) = domainGrowthFunction(tau,iDomain);
                end
                
                tmpAllBoundaries = [-flipud(cumsum(tmpDomainSizes));0;cumsum(tmpDomainSizes)];  % Temporary boundary locations.
                relevantX = validX*domainGrowthFunction(tau,1)/domainSizes(1);                  % Temporary relevant x values.
                relScale =  domainGrowthFunction(tau,1)/domainGrowthFunction(t,1);              % Ratio of current domain size to temporary domain size.
                
                % Calculate the appropriate weights.
                fluxScale = relativeD(1)*sqrt(relativeTransformFunction(tau,t,1))*domainGrowthFunction(tau,2)/domainGrowthFunction(t,2);
                fluxScale = fluxScale/(relativeD(1)*sqrt(relativeTransformFunction(tau,t,1))*domainGrowthFunction(tau,1)/domainGrowthFunction(t,1) + ...
                    relativeD(2)*sqrt(relativeTransformFunction(tau,t,2))*domainGrowthFunction(tau,2)/domainGrowthFunction(t,2));
                
                % Enumerate convolution terms.
                convolutionTerms(iTime,indices) = convolutionTerms(iTime,indices) + ...
                    leftFlux(jTime) * dt * ...
                    1/sqrt(pi*Dn*tauTimeTransform) * relScale * fluxScale * ...
                    exp(-(relevantX-tmpAllBoundaries(2)).^2/(4*Dn*tauTimeTransform)) + ...
                    rightFlux(jTime) * dt * ...
                    1/sqrt(pi*Dn*tauTimeTransform) * relScale * fluxScale * ...
                    exp(-(relevantX-tmpAllBoundaries(4)).^2/(4*Dn*tauTimeTransform));
                
                % Include terms from the summation.
                for iTerm = 1:nTerms
                    convolutionTerms(iTime,indices) = convolutionTerms(iTime,indices) + ...
                        (-1)^iTerm * leftFlux(jTime) * dt * ...
                        1/sqrt(pi*Dn*tauTimeTransform)* relScale * fluxScale * ( ...
                        exp(-(relevantX-(-2*iTerm*tmpAllBoundaries(1) + ...
                        (-1)^iTerm*tmpAllBoundaries(2))).^2/(4*Dn*tauTimeTransform)) + ...
                        exp(-(relevantX-(2*iTerm*tmpAllBoundaries(1) + ...
                        (-1)^iTerm*tmpAllBoundaries(2))).^2/(4*Dn*tauTimeTransform))) + ...
                        (-1)^iTerm * rightFlux(jTime) * dt * ...
                        1/sqrt(pi*Dn*tauTimeTransform)* relScale * fluxScale * ( ...
                        exp(-(relevantX-(-2*iTerm*tmpAllBoundaries(5) + ...
                        (-1)^iTerm*tmpAllBoundaries(4))).^2/(4*Dn*tauTimeTransform)) + ...
                        exp(-(relevantX-(2*iTerm*tmpAllBoundaries(5) + ...
                        (-1)^iTerm*tmpAllBoundaries(4))).^2/(4*Dn*tauTimeTransform)));
                end
                
            end
        end
        
        % Update solution profile, total mass, and mass in the inner
        % domain.
        C(iTime,indices) = C(iTime,indices) + convolutionTerms(iTime,indices);
        totalMass(iStart,iTime) = totalMass(iStart,iTime) + trapz(C(iTime,indices))*dx(2,iTime);
        middleMass(iStart,iTime) = middleMass(iStart,iTime) + trapz(C(iTime,indices))*dx(2,iTime);
        
        %% Solution for the -L_2(t)-L_1(t) < x < -L_1(t) domain.
        
        Dn = D(2);                                                              % Relevant diffusivity.
        validX = linspace(allBoundaries(1),allBoundaries(2),nMOISpacePoints);   % Relevant x values.
        dx(1,iTime) = validX(2)-validX(1);                                      % Current space step.
        indices = 1:nMOISpacePoints;                                            % Indices for storage of solution.
        currentValidX(indices) = validX;                                        % Relevant x values.
        
        % Do not calculate solution if only the final splitting
        % probabilities are wanted.
        if strcmpi(onlySplitting,'no')
            % Enumerate convolution terms by integrating over time. 
            for jTime = 1:iTime-1
                tmpDomainSizes = zeros(nDomains,1);                             % Temporary domain sizes.
                tau = timeValues(jTime);                                        % Temporary time.
                tauTimeTransform = relativeTransformFunction(tau,t,2);          % Temporary time transformation.
                
                % Calculate temporary domain size.
                for iDomain = 1:nDomains
                    tmpDomainSizes(iDomain) = domainGrowthFunction(tau,iDomain);
                end
                
                tmpAllBoundaries = [-flipud(cumsum(tmpDomainSizes));0;cumsum(tmpDomainSizes)];  % Temporary boundary locations.
                % Temporary relevant x values.
                relevantX = (validX + domainSizes(1))/domainSizes(2)*domainGrowthFunction(tau,2) - domainGrowthFunction(tau,1);
                relScale = domainGrowthFunction(tau,2)/domainGrowthFunction(t,2);               % Ratio of temporary domain size to current domain size.
                % Calculate the appropriate weights.
                fluxScale = relativeD(2)*sqrt(relativeTransformFunction(tau,t,2))*domainGrowthFunction(tau,1)/domainGrowthFunction(t,1);
                fluxScale = fluxScale/(relativeD(1)*sqrt(relativeTransformFunction(tau,t,1))*domainGrowthFunction(tau,1)/domainGrowthFunction(t,1) + ...
                    relativeD(2)*sqrt(relativeTransformFunction(tau,t,2))*domainGrowthFunction(tau,2)/domainGrowthFunction(t,2));
                
                % Enumerate convolution terms.
                convolutionTerms(iTime,indices) = convolutionTerms(iTime,indices) + ...
                    leftFlux(jTime) * dt * ...
                    1/sqrt(pi*Dn*tauTimeTransform) * relScale * fluxScale  * ...
                    exp(-(relevantX-tmpAllBoundaries(2)).^2/(4*Dn*tauTimeTransform));
                
                % Include terms from the summation.
                for iTerm = 1:nTerms
                    convolutionTerms(iTime,indices) = convolutionTerms(iTime,indices) + ...
                        (-1)^iTerm * leftFlux(jTime) * dt * ...
                        1/sqrt(pi*Dn*tauTimeTransform) * relScale * fluxScale * ( ...
                        exp(-(relevantX-(-2*iTerm*tmpAllBoundaries(1) + ...
                        (-1)^iTerm*tmpAllBoundaries(2))).^2/(4*Dn*tauTimeTransform)) + ...
                        exp(-(relevantX-(2*iTerm*tmpAllBoundaries(1) + ...
                        (-1)^iTerm*tmpAllBoundaries(2))).^2/(4*Dn*tauTimeTransform)));
                end
                
            end
        end
        
        % Update solution profile, total mass, and mass in the left outer
        % domain.
        C(iTime,indices) = C(iTime,indices) + convolutionTerms(iTime,indices);
        totalMass(iStart,iTime) = totalMass(iStart,iTime) + trapz(C(iTime,indices))*dx(1,iTime);
        leftMass(iStart,iTime) = leftMass(iStart,iTime) + trapz(C(iTime,indices))*dx(1,iTime);
        
        %% Solution for the 0 < x < L_1(t) domain.
        
        Dn = D(1);                                                              % Relevant diffusivity.
        validX = linspace(allBoundaries(3),allBoundaries(4),nMOISpacePoints);   % Relevant x values.
        dx(3,iTime) = validX(2)-validX(1);                                      % Current space step.
        indices = 2*nMOISpacePoints+1:3*nMOISpacePoints;                        % Indices for storage of solution.
        currentValidX(indices) = validX;                                        % Relevant x values.
        
        % Do not calculate solution if only the final splitting
        % probabilities are wanted.
        if strcmpi(onlySplitting,'no')
            % Enumerate convolution terms by integrating over time. 
            for jTime = 1:iTime-1
                tmpDomainSizes = zeros(nDomains,1);                             % Temporary domain size.
                tau = timeValues(jTime);                                        % Temporary time.
                tauTimeTransform = relativeTransformFunction(tau,t,1);          % Temporary time transformation.
                
                % Calculate temporary domain size.
                for iDomain = 1:nDomains
                    tmpDomainSizes(iDomain) = domainGrowthFunction(tau,iDomain);
                end
                
                tmpAllBoundaries = [-flipud(cumsum(tmpDomainSizes));0;cumsum(tmpDomainSizes)];  % Temporary boundary locations.
                relevantX = validX*domainGrowthFunction(tau,1)/domainSizes(1);                  % Temporary relevant x values.
                relScale =  domainGrowthFunction(tau,1)/domainGrowthFunction(t,1);              % Ratio of temporary domain size to current domain size.
                % Calculate the appropriate weights.
                fluxScale = relativeD(1)*sqrt(relativeTransformFunction(tau,t,1))*domainGrowthFunction(tau,2)/domainGrowthFunction(t,2);
                fluxScale = fluxScale/(relativeD(1)*sqrt(relativeTransformFunction(tau,t,1))*domainGrowthFunction(tau,1)/domainGrowthFunction(t,1) + ...
                    relativeD(2)*sqrt(relativeTransformFunction(tau,t,2))*domainGrowthFunction(tau,2)/domainGrowthFunction(t,2));
                
                % Enumerate convolution terms.
                convolutionTerms(iTime,indices) = convolutionTerms(iTime,indices) + ...
                    rightFlux(jTime) * dt * ...
                    1/sqrt(pi*Dn*tauTimeTransform) * relScale * fluxScale * ...
                    exp(-(relevantX-tmpAllBoundaries(4)).^2/(4*Dn*tauTimeTransform)) + ...
                    leftFlux(jTime) * dt * ...
                    1/sqrt(pi*Dn*tauTimeTransform) * relScale * fluxScale * ...
                    exp(-(relevantX-tmpAllBoundaries(2)).^2/(4*Dn*tauTimeTransform));
                
                % Include terms from the summation.
                for iTerm = 1:nTerms
                    convolutionTerms(iTime,indices) = convolutionTerms(iTime,indices) + ...
                        (-1)^iTerm * rightFlux(jTime) * dt * ...
                        1/sqrt(pi*Dn*tauTimeTransform)* relScale * fluxScale * ( ...
                        exp(-(relevantX-(-2*iTerm*tmpAllBoundaries(5) + ...
                        (-1)^iTerm*tmpAllBoundaries(4))).^2/(4*Dn*tauTimeTransform)) + ...
                        exp(-(relevantX-(2*iTerm*tmpAllBoundaries(5) + ...
                        (-1)^iTerm*tmpAllBoundaries(4))).^2/(4*Dn*tauTimeTransform))) + ...
                        (-1)^iTerm * leftFlux(jTime) * dt * ...
                        1/sqrt(pi*Dn*tauTimeTransform)* relScale * fluxScale * ( ...
                        exp(-(relevantX-(-2*iTerm*tmpAllBoundaries(1) + ...
                        (-1)^iTerm*tmpAllBoundaries(2))).^2/(4*Dn*tauTimeTransform)) + ...
                        exp(-(relevantX-(2*iTerm*tmpAllBoundaries(1) + ...
                        (-1)^iTerm*tmpAllBoundaries(2))).^2/(4*Dn*tauTimeTransform)));
                end
            end
        end
        
        % Update solution profile, total mass, and mass in the left outer
        % domain.
        C(iTime,indices) = C(iTime,indices) + convolutionTerms(iTime,indices);
        totalMass(iStart,iTime) = totalMass(iStart,iTime) + trapz(C(iTime,indices))*dx(3,iTime);
        middleMass(iStart,iTime) = middleMass(iStart,iTime) + trapz(C(iTime,indices))*dx(3,iTime);
        
        %% Solution for the L_1(t) < x L_1(t)+L_2(t) domain.
        
        Dn = D(2);                                                              % Relevant diffusivity.
        validX = linspace(allBoundaries(4),allBoundaries(5),nMOISpacePoints);   % Relevant x values.
        dx(4,iTime) = validX(2)-validX(1);                                      % Current space step.
        indices = 3*nMOISpacePoints+1:4*nMOISpacePoints;                        % Indices for storage of solution.
        currentValidX(indices) = validX;                                        % Relevant x values.
        
        % Do not calculate solution if only the final splitting
        % probabilities are wanted.
        if strcmpi(onlySplitting,'no')
            % Enumerate convolution terms by integrating over time.
            for jTime = 1:iTime-1
                tmpDomainSizes = zeros(nDomains,1);                             % Temporary domain sizes.
                tau = timeValues(jTime);                                        % Temporary time.
                tauTimeTransform = relativeTransformFunction(tau,t,2);          % Temporary time transformation.
                
                % Calculate temporary domain size.
                for iDomain = 1:nDomains
                    tmpDomainSizes(iDomain) = domainGrowthFunction(tau,iDomain);
                end
                
                tmpAllBoundaries = [-flipud(cumsum(tmpDomainSizes));0;cumsum(tmpDomainSizes)];      % Temporary boundary locations.
                % Temporary relevant x values.
                relevantX = (validX - domainSizes(1))/domainSizes(2)*domainGrowthFunction(tau,2) + domainGrowthFunction(tau,1);
                relScale = domainGrowthFunction(tau,2)/domainGrowthFunction(t,2);                   % Ratio of temporary domain size to current domain size.
                % Calculate the appropriate weights.
                fluxScale = relativeD(2)*sqrt(relativeTransformFunction(tau,t,2))*domainGrowthFunction(tau,1)/domainGrowthFunction(t,1);
                fluxScale = fluxScale/(relativeD(1)*sqrt(relativeTransformFunction(tau,t,1))*domainGrowthFunction(tau,1)/domainGrowthFunction(t,1) + ...
                    relativeD(2)*sqrt(relativeTransformFunction(tau,t,2))*domainGrowthFunction(tau,2)/domainGrowthFunction(t,2));              
                
                % Enumerate convolution terms.
                convolutionTerms(iTime,indices) = convolutionTerms(iTime,indices) + ...
                    rightFlux(jTime) * dt * ...
                    1/sqrt(pi*Dn*tauTimeTransform) * relScale * fluxScale * ...
                    exp(-(relevantX-tmpAllBoundaries(4)).^2/(4*Dn*tauTimeTransform));
                
                % Include terms from the summation
                for iTerm = 1:nTerms
                    convolutionTerms(iTime,indices) = convolutionTerms(iTime,indices) + ...
                        (-1)^iTerm * rightFlux(jTime) * dt * ...
                        1/sqrt(pi*Dn*tauTimeTransform)* relScale * fluxScale * ( ...
                        exp(-(relevantX-(-2*iTerm*tmpAllBoundaries(5) + ...
                        (-1)^iTerm*tmpAllBoundaries(4))).^2/(4*Dn*tauTimeTransform)) + ...
                        exp(-(relevantX-(2*iTerm*tmpAllBoundaries(5) + ...
                        (-1)^iTerm*tmpAllBoundaries(4))).^2/(4*Dn*tauTimeTransform)));
                end
                
            end
        end
        
        % Update solution profile, total mass, and mass in the left outer
        % domain.
        C(iTime,indices) = C(iTime,indices) + convolutionTerms(iTime,indices);
        totalMass(iStart,iTime) = totalMass(iStart,iTime) + trapz(C(iTime,indices))*dx(4,iTime);
        rightMass(iStart,iTime) = rightMass(iStart,iTime) + trapz(C(iTime,indices))*dx(4,iTime);
        
        %% Calculate splitting probabilities
        
        % Do not calculate solution if only the final splitting
        % probabilities are wanted.
        if strcmpi(onlySplitting,'no')
            % Loop over summation terms.
            for iTerm = 0:nTerms
                % Loop over time.
                for jTime = 1:iTime-1
                    phi = timeValues(jTime);                                            % Temporary time.
                    tmpDomainSizes = zeros(nDomains,1);                                 % Temporary domain sizes.
                    phiTimeTransform_Domain1 = relativeTransformFunction(phi,t,1);      % Temporary time transform for inner domain.
                    phiTimeTransform_Domain2 = relativeTransformFunction(phi,t,2);      % Temporary time transform for outer domain.
                    
                    % Calculate temporary domain sizes.
                    for iDomain = 1:nDomains
                        tmpDomainSizes(iDomain) = domainGrowthFunction(phi,iDomain);
                    end
                    
                    % Evaluate the right and left splitting probabilities.
                    %% WARNING - ONLY CORRECT IF STATED CONDITION HOLDS FOR INTEGRAL TO SIMPLIFY %%
                    rightSplit(iStart,iTime) = rightSplit(iStart,iTime) + 2*sqrt(D(1))/(sqrt(D(1))+sqrt(D(2))) * ...
                        rightFlux(jTime)*(-1)^iTerm*erfc(((2*iTerm+1)*(tmpDomainSizes(1)+tmpDomainSizes(2)) - (-1)^iTerm*tmpDomainSizes(1)) / ...
                        sqrt(4*D(2)*phiTimeTransform_Domain2)) * dt;
                    leftSplit(iStart,iTime) = leftSplit(iStart,iTime) + 2*sqrt(D(1))/(sqrt(D(1))+sqrt(D(2))) * ...
                        leftFlux(jTime)*(-1)^iTerm*erfc(((2*iTerm+1)*(tmpDomainSizes(1)+tmpDomainSizes(2)) - (-1)^iTerm*tmpDomainSizes(1)) / ...
                        sqrt(4*D(2)*phiTimeTransform_Domain2)) * dt;
                    
                end
            end
        end
        
        % If only final splitting probabilities are desired, calculate
        % splitting probability only at final time step.
        if strcmpi(onlySplitting,'yes') && iTime == nMOISavePoints
            % Loop over summation terms.
            for iTerm = 0:nTerms
                % Loop over time.
                for jTime = 1:iTime-1
                    phi = timeValues(jTime);                                            % Temporary time.
                    tmpDomainSizes = zeros(nDomains,1);                                 % Temporary domain sizes.
                    phiTimeTransform_Domain1 = relativeTransformFunction(phi,t,1);      % Temporary time transformation for inner domain.
                    phiTimeTransform_Domain2 = relativeTransformFunction(phi,t,2);      % Temporary time transformation for outer domain.
                    
                    % Calculate current domain sizes.
                    for iDomain = 1:nDomains
                        tmpDomainSizes(iDomain) = domainGrowthFunction(phi,iDomain);
                    end
                    
                    % Evaluate the right and left splitting probabilities.
                    %% WARNING - ONLY CORRECT IF STATED CONDITION HOLDS FOR INTEGRAL TO SIMPLIFY %%
                    rightSplit(iStart,iTime) = rightSplit(iStart,iTime) + 2*sqrt(D(1))/(sqrt(D(1))+sqrt(D(2))) * ...
                        rightFlux(jTime)*(-1)^iTerm*erfc(((2*iTerm+1)*(tmpDomainSizes(1)+tmpDomainSizes(2)) - (-1)^iTerm*tmpDomainSizes(1)) / ...
                        sqrt(4*D(2)*phiTimeTransform_Domain2)) * dt;
                    leftSplit(iStart,iTime) = leftSplit(iStart,iTime) + 2*sqrt(D(1))/(sqrt(D(1))+sqrt(D(2))) * ...
                        leftFlux(jTime)*(-1)^iTerm*erfc(((2*iTerm+1)*(tmpDomainSizes(1)+tmpDomainSizes(2)) - (-1)^iTerm*tmpDomainSizes(1)) / ...
                        sqrt(4*D(2)*phiTimeTransform_Domain2)) * dt;
                    
                end
            end
        end
    end
    
    % Store final solution profile.
    C_finalStore(iStart,:) = C(iTime,:);
end