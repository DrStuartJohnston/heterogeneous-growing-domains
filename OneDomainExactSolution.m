%% Script to enumerate exact single omain solutions for "Exact solutions for diffusive
%% transport on heterogeneous growing domains" by Johnston and Simpson, 2023.

totalMass = zeros(nStartingPoints,nMOISavePoints);                              % Total mass across the domain.
leftSplit = zeros(nStartingPoints,nMOISavePoints);                              % Left splitting probability.
rightSplit = zeros(nStartingPoints,nMOISavePoints);                             % Right splitting probability.
C_finalStore = zeros(nStartingPoints,2*nMOISpacePoints);                        % Store final solution profile.

% Loop over each starting location.
for iStart = 1:nStartingPoints
    C = zeros(nMOISavePoints,2*nMOISpacePoints);                                % All solution profiles for a single starting point.
    rightFlux = zeros(1,nMOISavePoints);                                        % Flux through the right outer boundary.
    leftFlux = zeros(1,nMOISavePoints);                                         % Flux through the left outer boundary.         
    timeValues = linspace(tEnd/nMOISavePoints,tEnd,nMOISavePoints);             % Time points where solution is enumerated.
    domainSizes = 0;                                                            % Current domain size.
    dx = zeros(2,nMOISavePoints);                                               % Current node spacing (space step).
    
    % Loop over time points.
    for iTime = 1:nMOISavePoints

        t = timeValues(iTime);                                                  % Current time.
        
        % Calculate current domain sizes.
        for iDomain = 1:nDomains
            domainSizes(iDomain) = domainGrowthFunction(t,iDomain);
        end
        
        allBoundaries = [-flipud(cumsum(domainSizes));0;cumsum(domainSizes)];   % Current boundary locations.
        timeTransform = transformFunction(t,1);                                 % Perform time transformation.
        currentValidX = zeros(1,2*nMOISpacePoints);                             % Current relevant x values.
        
        %% Solution for -L_1(t) < x < 0 Domain.
        
        xScale = L0(1)/domainSizes(1);                                          % Rescaling of x.
        Dn = D(1);                                                              % Relevant diffusivity.
        validX = linspace(allBoundaries(1),allBoundaries(2),nMOISpacePoints);   % Relevant x values.
        dx(1,iTime) = validX(2)-validX(1);                                      % Current space step.
        indices = 1:nMOISpacePoints;                                            % Indices for storage of solution.
        currentValidX(indices) = validX;                                        % Relevant x values.
        
        % Enumerate solution in the domain.
        C(iTime,indices) = 1/sqrt(4*pi*Dn*timeTransform)*xScale*exp(-(validX*xScale-x0(iStart)).^2/(4*Dn*timeTransform));
        for iTerm = 1:nTerms
            C(iTime,indices) = C(iTime,indices) + 1/sqrt(4*pi*Dn*timeTransform)*xScale*( ...
                (-1)^iTerm*exp(-(validX*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) + ...
                (-1)^iTerm*exp(-(validX*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)));
        end
        
        % Enumerate flux through the outer left boundary.
        leftFlux(iTime) = 2*Dn*(allBoundaries(1)*xScale-x0(iStart))/(sqrt(pi)*(4*Dn*timeTransform)^(3/2))*xScale^2* ...
            exp(-(allBoundaries(1)*xScale-x0(iStart)).^2/(4*Dn*timeTransform));
        for iTerm = 1:nTerms
            leftFlux(iTime) = leftFlux(iTime) + 1/(sqrt(pi)*(4*Dn*timeTransform)^(3/2))*xScale^2*( ...
                (-1)^iTerm*exp(-(allBoundaries(1)*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) * ...
                2*Dn*(allBoundaries(1)*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))) + ...
                (-1)^iTerm*exp(-(allBoundaries(1)*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) * ...
                2*Dn*(allBoundaries(1)*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))));
        end
        
        % Store total mass.
        totalMass(iStart,iTime) = totalMass(iStart,iTime) + trapz(C(iTime,indices))*dx(1,iTime);
        
        %% Solution for 0 < x < L_1(t) Domain.
        
        xScale = L0(1)/domainSizes(1);                                          % Rescaling of x.
        Dn = D(1);                                                              % Relevant diffusivity.
        validX = linspace(allBoundaries(2),allBoundaries(3),nMOISpacePoints);   % Relevant x values.   
        dx(2,iTime) = validX(2)-validX(1);                                      % Current space step.
        indices = nMOISpacePoints+1:2*nMOISpacePoints;                          % Indices for storage of solution.    
        currentValidX(indices) = validX;                                        % Relevant x values.
        
        % Enumerate solution in the domain.
        C(iTime,indices) = 1/sqrt(4*pi*Dn*timeTransform)*xScale*exp(-(validX*xScale-x0(iStart)).^2/(4*Dn*timeTransform));
        for iTerm = 1:nTerms
            C(iTime,indices) = C(iTime,indices) + 1/sqrt(4*pi*Dn*timeTransform)*xScale*( ...
                (-1)^iTerm*exp(-(validX*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) + ...
                (-1)^iTerm*exp(-(validX*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)));
        end
        
        % Enumerate flux through the outer right boundary.
        rightFlux(iTime) = 2*Dn*(allBoundaries(3)*xScale-x0(iStart))/(sqrt(pi)*(4*Dn*timeTransform)^(3/2))*xScale^2* ...
            exp(-(allBoundaries(3)*xScale-x0(iStart)).^2/(4*Dn*timeTransform));
        for iTerm = 1:nTerms
            rightFlux(iTime) = rightFlux(iTime) + 1/(sqrt(pi)*(4*Dn*timeTransform)^(3/2))*xScale^2*( ...
                (-1)^iTerm*exp(-(allBoundaries(3)*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) * ...
                2*Dn*(allBoundaries(3)*xScale-(-2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))) + ...
                (-1)^iTerm*exp(-(allBoundaries(3)*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))).^2/(4*Dn*timeTransform)) * ...
                2*Dn*(allBoundaries(3)*xScale-(2*iTerm*L0(1)+(-1)^iTerm*x0(iStart))));
        end
        
        % Store total mass.
        totalMass(iStart,iTime) = totalMass(iStart,iTime) + trapz(C(iTime,indices))*dx(2,iTime);
        
        % Calculate left and right splitting probabilities.
        for i = 0:nTerms
            leftSplit(iStart,iTime) = leftSplit(iStart,iTime) + (-1)^i*erfc(((2*i+1)*L0(1)+(-1)^i*x0(iStart))./(sqrt(4*D(1)*timeTransform)));
            rightSplit(iStart,iTime) = rightSplit(iStart,iTime) + (-1)^i*erfc(((2*i+1)*L0(1)-(-1)^i*x0(iStart))./(sqrt(4*D(1)*timeTransform)));
        end
    end
    
    % Store final solution profile.
    C_finalStore(iStart,:) = C(iTime,:);
end
