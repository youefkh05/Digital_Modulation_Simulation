clear; clc; close all;

%--------Part 1----------
% ========================
% Simulation Parameters
% ========================
bits_Num = 48;                                  % Number of bits to transmit
mod_types = {'BPSK', 'QPSK', '8PSK', '16-QAM'}; % Cell array of modulation types

% Generate random bits (same for all modulations for fair comparison)
Tx_bits = randi([0 1], 1, bits_Num);

% Loop through all modulation types
for mod_idx = 1:length(mod_types)
    mod_type = mod_types{mod_idx};
    
    fprintf('\n=== Testing %s Modulation ===\n', mod_type);
    
    % ========================
    % 1. Mapping (Modulation)
    % ========================
    [tx_symbols, constellation] = mapper(Tx_bits, mod_type);
    
    % ========================
    % 2. Display Constellation
    % ========================
    drawConstellation(constellation, mod_type);
    title(sprintf('%s Constellation', mod_type));
    
    % ========================
    % 3. Add Channel Noise
    % ========================
    %rx_symbols = awgn(tx_symbols, SNR_dB, 'measured');
    rx_symbols = tx_symbols;
    % ========================
    % 4. Demapping (Demodulation)
    % ========================
    Rx_bits = demapper(rx_symbols, mod_type);
    
    % ========================
    % 5. Display Results
    % ========================
    % Calculate BER
    [BER, bit_errors] = calculateBER(Tx_bits, Rx_bits);
    
    % Display input/output comparison
    fprintf('Original bits:\n');
    disp(reshape(Tx_bits, 16, [])'); % Display in 16-bit groups
    
    fprintf('Received bits:\n');
    disp(reshape(Rx_bits(1:bits_Num), 16, [])'); % Display in 16-bit groups
    
    fprintf('Bit errors: %d\n', bit_errors);
    fprintf('BER: %.2e\n\n', BER);
end


% ========================
% Functions
% ========================

function [Tx_Vector, Table] = mapper(bits, mod_type)
    % MAPPER Digital modulation mapper with explicit symbol table
    % Inputs:
    %   bits     - Binary input array (row vector)
    %   mod_type - 'BPSK', 'QPSK', '8PSK', 'BFSK', '16QAM'
    % Outputs:
    %   Tx_Vector - Complex modulated symbols
    %   Table     - Constellation points (M-ary symbols)

    % Ensure bits are row vector
    bits = bits(:)';
    
    % Define modulation parameters
    switch upper(mod_type)
        case 'BPSK'
            n = 1;  % bits per symbol
            M = 2;   % constellation size
            Table = [-1, 1];  % BPSK symbols (real)
            
        case 'QPSK'
            n = 2;
            M = 4;
            Table = [-1-1j, -1+1j, 1-1j, 1+1j];  % QPSK symbols
            
        case '8PSK'
            n = 3;
            M = 8;
            angles =[0, 1, 3, 2, 7, 6, 4, 5]*pi/4;  % Gray-coded 8PSK
            Table = exp(1j*angles);
            
        case 'BFSK'
            error('BFSK requires time-domain implementation (see alternative)');
            
        case '16-QAM'
            n = 4;
            M = 16;
            % 16-QAM with unit average power (normalized)
            Table = [-3-3j, -3-1j, -3+3j, -3+1j, ...
                     -1-3j, -1-1j, -1+3j, -1+1j, ...
                      3-3j,  3-1j,  3+3j,  3+1j, ...
                      1-3j,  1-1j,  1+3j,  1+1j];
            
        otherwise
            error('Unsupported modulation type: %s', mod_type);
    end
    
    % Pad bits if not multiple of n
    if mod(length(bits), n) ~= 0
        bits = [bits zeros(1, n - mod(length(bits), n))];
    end
    
    % Reshape into n-bit groups
    bit_groups = reshape(bits, n, [])';
    
    % Convert to decimal symbols (0 to M-1)
    Array_symbol = bi2de(bit_groups, 'left-msb') + 1;  % MATLAB uses 1-based indexing
    
    % Map to constellation points
    Tx_Vector = Table(Array_symbol);
end

function drawConstellation(Table, mod_type)
    % DRAWCOnSTELLATION Enhanced constellation visualization
    % Inputs:
    %   Table - Constellation points (complex numbers)
    %   mod_type - Modulation type ('BPSK', 'QPSK', etc.)
    %   show_regions - true to show colored regions, false for boundaries only
    
    if nargin < 3
        show_regions = true; % Default to showing regions
    end
    
    figure;
    hold on;
    
    % Ensure Table is column vector and get points
    Table = Table(:);
    points = [real(Table), imag(Table)];
    
    % Create grid for visualization
    x_range = linspace(min(points(:,1))-1, max(points(:,1))+1, 200);
    y_range = linspace(min(points(:,2))-1, max(points(:,2))+1, 200);
    [x_grid, y_grid] = meshgrid(x_range, y_range);
    grid_points = x_grid(:) + 1j*y_grid(:);
    
    % =============================================
    % 1. Decision Visualization
    % =============================================
    if length(Table) > 2  % Voronoi needs at least 3 points
        [vx, vy] = voronoi(points(:,1), points(:,2));
        plot(vx, vy, 'k-', 'LineWidth', 1.5);
    else
        % For BPSK, draw simple decision boundary
        plot([0 0], ylim, 'k--', 'LineWidth', 1.5);
    end
       
   
    % =============================================
    % 2. Constellation Points
    % =============================================
    scatter(points(:,1), points(:,2), 100, 'filled', 'k');
    
    % =============================================
    % 3. Binary Labels
    % =============================================
    switch upper(mod_type)
        case 'BPSK'
            n = 1; 
        case 'QPSK'
            n = 2; 
        case '8PSK'
            n = 3; 
        case {'16QAM', '16-QAM'}
            n = 4;
        otherwise
            error('Unsupported modulation type');
    end
    
    for i = 1:length(Table)
        bin_str = dec2bin(i-1, n);
        % Position text slightly offset from the point
        text(real(Table(i)) + 0.05, imag(Table(i)) + 0.05, bin_str, ...
            'FontSize', 10, 'Color', 'r');
    end
    
    % =============================================
    % 4. Plot Formatting
    % =============================================
    title(sprintf('%s Constellation', mod_type));
    xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
    grid on;
    axis equal;
    
    % Center axes
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    
    % Set axis limits
    max_val = max([abs(points(:))]) * 1.3;
    xlim([-max_val, max_val]);
    ylim([-max_val, max_val]);
    
    
    hold off;
end

function [received_bits] = demapper(received_symbols, mod_type)
    % DEMAPPER Digital demodulation demapper
    % Inputs:
    %   received_symbols - Complex received symbols
    %   mod_type        - Modulation type ('BPSK', 'QPSK', etc.)
    % Output:
    %   received_bits   - Demodulated bit stream
    
    % Get constellation table from mapper
    [~, Table] = mapper([1], mod_type);
    
    % Determine bits per symbol
    switch upper(mod_type)
        case 'BPSK'
            n = 1;
        case 'QPSK'
            n = 2;
        case '8PSK'
            n = 3;
        case {'16QAM', '16-QAM'}
            n = 4;
        otherwise
            error('Unsupported modulation type');
    end
    
    % Initialize output bits
    received_bits = zeros(1, length(received_symbols)*n);
    
    % Demodulate each symbol
    for i = 1:length(received_symbols)
        % Find nearest constellation point
        [~, idx] = min(abs(received_symbols(i) - Table));
        
        % Convert to binary (0-based index)
        bin_str = dec2bin(idx-1, n);
        
        % Store bits
        received_bits((i-1)*n+1:i*n) = bin_str - '0';
    end
end

function [BER, bit_errors] = calculateBER(original_bits, received_bits)
    % CALCULATEBER Compute Bit Error Rate between original and received bits
    % Inputs:
    %   original_bits - Transmitted bit sequence (1D array)
    %   received_bits - Received/demodulated bit sequence (1D array)
    % Outputs:
    %   BER - Bit Error Rate (ratio of incorrect bits)
    %   bit_errors - Absolute number of bit errors
    
    % Ensure both inputs are row vectors
    original_bits = original_bits(:)';
    received_bits = received_bits(:)';
    
    % Trim received bits if longer (due to padding)
    if length(received_bits) > length(original_bits)
        received_bits = received_bits(1:length(original_bits));
    end
    
    % Calculate errors
    bit_errors = sum(original_bits ~= received_bits);
    BER = bit_errors / length(original_bits);
    
    % Display results
    %fprintf('Bit errors: %d\n', bit_errors);
    %fprintf('BER: %.2e\n', BER);
end