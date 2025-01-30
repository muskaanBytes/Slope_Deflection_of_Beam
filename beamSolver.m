% Prompt the user to specify the length of the beam
while true
    BeamLength = input('Enter the length of the beam: ');
    if BeamLength > 0
        break;
    else
        disp('Error: Beam length must be positive. Please enter a valid length.');
    end
end

% Prompt the user to specify the number of Concentrated Loads
while true
    numConcLoads = input('Enter the number of Concentrated Loads: ');
    if numConcLoads >= 0 
        break;
    else
        disp('Error: Number of concentrated loads must be non-negative. Please enter a valid number.');
    end
end

% Prompt the user to specify the number of Uniformly distributed Loads
while true
    numUdLoads = input('Enter the number of Uniformly Distributed Loads: ');
    if numUdLoads >= 0
        break;
    else
        disp('Error: Number of UDLs must be non-negative. Please enter a valid number.');
    end
end

% Initialize arrays to store load values, distances, and lengths
ConcLoadValues = zeros(1, numConcLoads);
ConcLoadDistance = zeros(1, numConcLoads);
ConcLoadMoment = zeros(1, numConcLoads);
UdLoadValues = zeros(1, numUdLoads);
UdLoadLength = zeros(1, numUdLoads);
UdLoadDistance = zeros(1, numUdLoads);
UdloadMoment = zeros(1, numUdLoads);
MaxmBendingMomentDistance = 0;
MaxmBendingMoment= 0;


% Prompt the user to input values for each concentrated load and its units
for i = 1:numConcLoads
    ConcLoadValues(i) = input(['Enter the value for Concentrated Load ' num2str(i) ' with sign (downward is negative) : ']);
end
while true
    ConcLoadUnits = input('Enter the units for Concentrated Loads (KN or N): ', 's'); % 's' reads input as string
    if strcmpi(ConcLoadUnits, 'KN') || strcmpi(ConcLoadUnits, 'N')
        break;
    else
        disp('Error: Invalid unit. Please enter either "KN" or "N".');
    end
end

% Prompt the user to specify the distance of Concentrated Loads from LHS
for i = 1:numConcLoads
    while true
        ConcLoadDistance(i) = input(['Enter the distance of Concentrated Load ' num2str(i) ' from LHS: ']);
        if ConcLoadDistance(i) >= 0 && ConcLoadDistance(i) <= BeamLength
            break;
        else
            disp('Error: Load distance must be non-negative and not exceed beam length. Please enter a valid distance.');
        end
    end
end

% Prompt the user to input values for each UDL and ensure total length does not exceed beam length
for i = 1:numUdLoads
    while true
        UdLoadValues(i) = input(['Enter the value for Uniformly Distributed Load ' num2str(i) ' with sign (downward is negative) : ']);
        UdLoadDistance(i) = input(['Enter the distance of Uniformly Distributed Load ' num2str(i) ' from LHS: ']);

        % Check for negative distance
        if UdLoadDistance(i) >= 0
            UdLoadLength(i) = input(['Enter the length for Uniformly Distributed Load ' num2str(i) ': ']);
            % Check for negative length
            if UdLoadLength(i) >= 0
                totalUdLoadLength = UdLoadDistance(i) + UdLoadLength(i);
                if totalUdLoadLength <= BeamLength
                    break;
                else
                    disp('Error: Total length of UDLs exceeds beam length. Please enter valid lengths and distances.');
                end
            else
                disp('Error: Length of UDL must be non-negative. Please enter a valid length.');
            end
        else
            disp('Error: Distance from the left-hand side (LHS) must be non-negative. Please enter a valid distance.');
        end
    end
end
while true
    UdLoadUnits = input('Enter the units for Uniformly Distributed Loads (KN/m or N/m): ', 's'); % 's' reads input as string
    if strcmpi(UdLoadUnits, 'KN/m') || strcmpi(UdLoadUnits, 'N/m')
        break;
    else
        disp('Error: Invalid unit. Please enter either "KN/m" or "N/m".');
    end
end

% Calculating Reactions at Supports
% y-axis equation
sumConcLoads = -sum(ConcLoadValues);
sumUdLoads = -sum(UdLoadValues .* UdLoadLength); % Sum of uniformly distributed loads

sumLoads = sumConcLoads + sumUdLoads;

%Moment Calculation
for i = 1:numConcLoads
    ConcLoadMoment(i) = -(ConcLoadValues(i) * ConcLoadDistance(i));
end


for i = 1:numUdLoads
    UdloadMoment(i) = -(UdLoadValues(i) * UdLoadLength(i) * (UdLoadLength(i) / 2 + UdLoadDistance(i)));
end

% Solving equations to find support reactions
syms R1 R2
eq1 = R1 + R2 == sumLoads;
eq2 = R2 * BeamLength == sum(ConcLoadMoment) + sum(UdloadMoment);

% Solve the equations
solution = solve([eq1, eq2], [R1, R2]);
R1 = double(solution.R1);
R2 = double(solution.R2);

disp(['R1 = ', num2str(R1), ' ', ConcLoadUnits]);
disp(['R2 = ', num2str(R2), ' ', ConcLoadUnits]);

% Create a new figure
figure;

% Plotting the shear force diagram
subplot(1, 2, 1); % Create subplot for Shear Force Diagram
x = linspace(0, BeamLength, 1000); % Define positions along the beam
V = zeros(size(x)); % Initialize shear force array

% Loop over each position along the beam
for i = 1:length(x)
    % Initialize shear force at this position without considering R2
    V(i) = R1;

    % Calculate contributions from concentrated loads
    for j = 1:numConcLoads
        if x(i) > ConcLoadDistance(j)
            V(i) = V(i) + ConcLoadValues(j); 
        end
    end

    % Calculate contributions from uniformly distributed loads
    for k = 1:numUdLoads
        if x(i) >= UdLoadDistance(k)
            % Determine the portion of the load acting at this position
            if x(i) - UdLoadDistance(k) <= UdLoadLength(k)
                load_portion = x(i) - UdLoadDistance(k);
            else
                load_portion = UdLoadLength(k);
            end
            % Add the contribution from this load to the shear force
            V(i) = V(i) + UdLoadValues(k) * load_portion;
        end
    end
    
    % Add R2 when reaching the end of the beam
    if x(i) == BeamLength
        V(i) = V(i) + R2;
    end
end

%Finding Maximum Bending Moment Distance 
for i = 2:length(x) % Start from the second index to avoid the first point
    if sign(V(i)) ~= sign(V(i-1)) % Check for sign change
        MaxmBendingMomentDistance = x(i);
        break; % Break the loop when the first sign change is detected
    end
end

% Plotting the shear force diagram
plot(x, V, 'b', 'LineWidth', 1.5);
xlabel('Position along the beam (m)'); % Append units to the axis label
ylabel(['Shear Force (' ConcLoadUnits ')']); % Append units to the axis label
title('Shear Force Diagram');
grid on;

% Plotting the Bending Moment diagram
subplot(1, 2, 2); % Create subplot for Bending Moment Diagram
M = zeros(size(x)); % Initialize bending moment array

% Loop over each position along the beam
for i = 1:length(x)
    % Initialize bending moment at this position without considering R2
    M(i) = R1 * x(i);

    % Calculate contributions from concentrated loads
    for j = 1:numConcLoads
        if x(i) > ConcLoadDistance(j)
            M(i) = M(i) + ConcLoadValues(j) * (x(i) - ConcLoadDistance(j)); 
        end
    end

    % Calculate contributions from uniformly distributed loads
    for k = 1:numUdLoads
        if x(i) >= UdLoadDistance(k)
            % Determine the portion of the load acting at this position
            if x(i) - UdLoadDistance(k) <= UdLoadLength(k)
                load_portion = x(i) - UdLoadDistance(k);
            else
                load_portion = UdLoadLength(k);
            end
            % Add the contribution from this load to the bending moment
            M(i) = M(i) + UdLoadValues(k) * load_portion * (x(i) - UdLoadDistance(k) - 0.5 * load_portion);
        end
   
    end
end
          % Finding Maximum Bending Moment
          for i = 1:length(x)
    if x(i) == MaxmBendingMomentDistance
        MaxmBendingMoment = M(i);
    end
          end
% Finding Point of Contraflexure (if bending moment doesn't change sign)
pointOfContraflexure = NaN; % Initialize with NaN (Not a Number)
for i = 2:(length(x)-1) % Start from the second index and avoid the last index
    if (M(i) > 0 && M(i+1) < 0) || (M(i) < 0 && M(i+1) > 0) && M(i) ~= 0 % Check for transition from positive to negative or vice versa
        pointOfContraflexure = x(i);
        break; % Break the loop when the first transition is detected
    end
end

   
% Display maximum bending moment distance and maximum bending moment 
disp(['Maximum Bending Moment Distance From LHS: ', num2str(MaxmBendingMomentDistance), ' m']);
disp(['Maximum Bending Moment: ', num2str(MaxmBendingMoment), ' ', UdLoadUnits]);
if isnan(pointOfContraflexure)
    disp('No point of contraflexure found.');
else
    disp(['Point of Contraflexure: ', num2str(pointOfContraflexure), ' m']);
end

% Plotting the Bending Moment diagram
plot(x, M, 'r', 'LineWidth', 1.5);
xlabel('Position along the beam (m)'); % Append units to the axis label
ylabel(['Bending Moment (' UdLoadUnits ')']); % Append units to the axis label
title('Bending Moment Diagram');
grid on;