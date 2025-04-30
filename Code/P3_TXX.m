clear; clc; close all;

%--------Part 1----------
% ========================
% Simulation Parameters
% ========================
%bits_Num = 6 * 2^10;                                  % Number of bits to transmit
bits_Num = 48;  
mod_types = {'BPSK', 'QPSK', '8PSK', '16-QAM'}; % Cell array of modulation types
SNR_db_range = -4:1:14;

% Generate random bits (same for all modulations for fair comparison)
Tx_bits = randi([0 1], 1, bits_Num);

% ========================
% Initialize storage matrices
% ========================

% Initialize rx_symbols_all as 2D cell matrix
% Rows: modulation types, Columns: SNR values
rx_symbols_all = cell(length(mod_types), length(SNR_db_range));

% Initialize storage for Energy Bits
Eb_all = cell(1, length(mod_types));

% Initialize storage for Error
BER_all = zeros(length(mod_types), length(SNR_db_range));
error_count_all = zeros(length(mod_types), length(SNR_db_range));

%   Loop through all modulation types
for mod_idx = 1:length(mod_types)
    mod_type = mod_types{mod_idx};
    
    fprintf('\n=== %s Modulation ===\n', mod_type);
    
    % ========================
    % 1. Mapping (Modulation)
    % ========================
    [tx_symbols, constellation,~,Eb] = mapper(Tx_bits, mod_type);
    
    % Store Energy of bit for this modulation type
    Eb_all{mod_idx} = Eb;
    
    % ========================
    % 2. Display Constellation
    % ========================
    drawConstellation(constellation, mod_type, 1);
    title(sprintf('%s Constellation', mod_type));
    
    % ========================
    % 3. Channel Transmission
    % ========================
    % Get noisy symbols for all SNR values
    rx_noisy_symbols = addAWGNChannel(SNR_db_range, tx_symbols, Eb);
    
    % Store in 2D cell matrix
    rx_symbols_all(mod_idx, :) = rx_noisy_symbols;
    
    % ========================
    % 4. Demapping (Demodulation)
    % ========================
    Rx_bits = demapper(rx_noisy_symbols, mod_type);
    
    % ========================
    % 5. Calculate and Store Results
    % ========================
    fprintf('\nSNR Results:\n');
    fprintf('------------\n');
    
    for snr_idx = 1:6:length(SNR_db_range)
        [BER_all(mod_idx, snr_idx), error_count_all(mod_idx, snr_idx)] = ...
            calculateBER(Tx_bits, Rx_bits{snr_idx});
        
        % Display results for each SNR
        fprintf('SNR: %6.1f dB | BER: %8.2e | Errors: %4d/%d\n', ...
            SNR_db_range(snr_idx), ...
            BER_all(mod_idx, snr_idx), ...
            error_count_all(mod_idx, snr_idx), ...
            length(Tx_bits));
    end
    

end

%   Display Noise
drawNoisyConstellations(rx_symbols_all, SNR_db_range, mod_types);

disp(size(BER_all));
disp(size(SNR_db_range));
disp(length(mod_types));


%   Graph BER Vs SNR
plot_BER_vs_SNR(BER_all, SNR_db_range, mod_types);

% ========================
% Functions
% ========================

function [Tx_Vector, Table, Eavg, Eb] = mapper(bits, mod_type)
    % MAPPER Digital modulation mapper with explicit symbol table and energy calculation
    % Inputs:
    %   bits     - Binary input array (row vector)
    %   mod_type - 'BPSK', 'QPSK', '8PSK', 'BFSK', '16-QAM'
    % Outputs:
    %   Tx_Vector - Complex modulated symbols
    %   Table     - Constellation points (M-ary symbols)
    %   Eavg      - Average symbol energy (normalized)
    %   Eb        - Energy per bit

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
    
    % Calculate average symbol energy
    Eavg = mean(abs(Table).^2);
    
    % Calculate average bit energy
    Eb =  Eavg / n;
    
    % Reshape into n-bit groups
    bit_groups = reshape(bits, n, [])';
    
    % Convert to decimal symbols (0 to M-1)
    Array_symbol = bi2de(bit_groups, 'left-msb') + 1;  % MATLAB uses 1-based indexing
    
    % Map to constellation points
    Tx_Vector = Table(Array_symbol);
end

function drawConstellation(Table, mod_type, showdetails)
    % DRAWCOnSTELLATION Enhanced constellation visualization
    % Inputs:
    %   Table - Constellation points (complex numbers)
    %   mod_type - Modulation type ('BPSK', 'QPSK', etc.)
    %   showdetails- true to show colored regions, false for boundaries only
    
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
    if showdetails == 1
        if length(Table) > 2  % Voronoi needs at least 3 points
            [vx, vy] = voronoi(points(:,1), points(:,2));
            plot(vx, vy, 'k-', 'LineWidth', 1.5);
        else
            % For BPSK, draw simple decision boundary
            plot([0 0], ylim, 'k--', 'LineWidth', 1.5);
        end
    end
   
    % =============================================
    % 2. Constellation Points
    % =============================================
    if showdetails == 1
        scatter(points(:,1), points(:,2), 100, 'filled', 'k');
    else
        scatter(points(:,1), points(:,2), 20, 'filled', 'k');
    end
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
   
    if showdetails == 1
        for i = 1:length(Table)
            bin_str = dec2bin(i-1, n);
            % Position text slightly offset from the point
            text(real(Table(i)) + 0.05, imag(Table(i)) + 0.05, bin_str, ...
                'FontSize', 10, 'Color', 'r');
        end
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

function drawNoisyConstellations(rx_symbols_all, SNR_db_range, mod_types)
    % DRAWNOISYCONSTELLATIONS Plot constellations with noisy received points
    % Inputs:
    %   rx_symbols_all - Cell array, rx_symbols_all{mod_idx, snr_idx}
    %   SNR_db_range   - Vector of SNR values (dB)
    %   mod_types      - Cell array of modulation type strings (e.g., {'BPSK', 'QPSK'})
    
    % Validate inputs
    if ~iscell(rx_symbols_all) || ~iscell(mod_types)
        error('rx_symbols_all and mod_types must be cell arrays.');
    end

    num_mods = numel(mod_types);
    num_snr = numel(SNR_db_range);
    
    for mod_idx = 1:num_mods
        mod_type = mod_types{mod_idx};
        
        % Generate constellation table for this modulation
        [~, Table] = mapper([1], mod_type);
        
        for snr_idx = 1:3:num_snr
            rx_symbols = rx_symbols_all{mod_idx, snr_idx};
            snr_db = SNR_db_range(snr_idx);          
            
            % Center axes
            ax = gca;
            ax.XAxisLocation = 'origin';
            ax.YAxisLocation = 'origin';
            
            % Plot decision regions and ideal points
            drawConstellation(Table, mod_type, 0);
            title(sprintf('%s Constellation at SNR = %d dB', mod_type, snr_db));
            xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
            grid on;
            axis equal;
            hold on;
            % Plot noisy received symbols
            scatter(real(rx_symbols), imag(rx_symbols), 10, 'b', 'filled', 'MarkerFaceAlpha', 0.4);
            
            % Set axis limits a bit bigger to fit noisy points
            max_val = 4;
            xlim([-max_val, max_val]);
            ylim([-max_val, max_val]);
            
            hold off;
        end
    end
end

function [received_bits] = demapper(received_symbols, mod_type)
    % DEMAPPER Digital demodulation demapper
    % Inputs:
    %   received_symbols - Complex received symbols (array or cell array)
    %   mod_type        - Modulation type ('BPSK', 'QPSK', etc.)
    % Output:
    %   received_bits   - Demodulated bit stream (array or cell array)
    
    % Check if input is cell array (multiple SNR cases)
    if iscell(received_symbols)
        % Process each SNR case
        received_bits = cell(size(received_symbols));
        for i = 1:numel(received_symbols)
            received_bits{i} = demodulate_symbols(received_symbols{i}, mod_type);
        end
    else
        % Single SNR case
        received_bits = demodulate_symbols(received_symbols, mod_type);
    end
end

function bits = demodulate_symbols(symbols, mod_type)
    % Helper function for actual demodulation
    
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
    bits = zeros(1, length(symbols)*n);
    
    % Demodulate each symbol
    for i = 1:length(symbols)
        % Find nearest constellation point
        [~, idx] = min(abs(symbols(i) - Table));
        
        % Convert to binary (0-based index)
        bin_str = dec2bin(idx-1, n);
        
        % Store bits
        bits((i-1)*n+1:i*n) = bin_str - '0';
    end
end

function noisy_signals = addAWGNChannel(SNR_range_db, clean_signal, Eb)
    % ADDAGWNCHANNEL General AWGN channel noise adder
    % Inputs:
    %   SNR_range_db - Array of SNR values in dB
    %   clean_signal - Input signal (vector or matrix)
    %   Eb - Energy per bit
    % Output:
    %   noisy_signals - Cell array of noisy signals for each SNR
    
    % Initialize output cell array
    noisy_signals = cell(length(SNR_range_db), 1);
    
    % Get size of input signal
    signal_size = size(clean_signal);
    
    % Process each SNR point
    for i = 1:length(SNR_range_db)
        % Convert SNR from dB to linear scale
        SNR_linear = 10^(SNR_range_db(i)/10);
        
        % Calculate noise power (N0)
        N0 = 1 / SNR_linear;
        
        % Generate proper noise
        if isreal(clean_signal)
            % Real noise for real signals
            noise = sqrt(Eb*N0/2) * randn(signal_size);
        else
            % Complex noise for complex signals
            noise = sqrt(Eb*N0/2) * (randn(signal_size) + 1j*randn(signal_size));
        end
        
        % Add noise to the signal
        noisy_signals{i} = clean_signal + noise;
    end
    
    % If only one SNR point was requested, return array instead of cell
    if length(SNR_range_db) == 1
        noisy_signals = noisy_signals{1};
    end
end

function [BER, bit_errors] = calculateBER(original_bits, received_bits)
    % CALCULATEBER Compute Bit Error Rate for single or multiple SNR cases
    % Inputs:
    %   original_bits - Transmitted bit sequence (1D array)
    %   received_bits - Received bits (1D array or cell array for multiple SNR)
    % Outputs:
    %   BER - Bit Error Rate (scalar or array matching received_bits input)
    %   bit_errors - Number of errors (scalar or array)
    
    % Ensure original bits are row vector
    original_bits = original_bits(:)';
    
    % Handle cell array input (multiple SNR cases)
    if iscell(received_bits)
        BER = zeros(size(received_bits));
        bit_errors = zeros(size(received_bits));
        
        for i = 1:numel(received_bits)
            [BER(i), bit_errors(i)] = calculateSingleBER(original_bits, received_bits{i});
        end
    else
        % Single SNR case
        [BER, bit_errors] = calculateSingleBER(original_bits, received_bits);
    end
end

function [BER, bit_errors] = calculateSingleBER(original_bits, received_bits)
    % Helper function for single SNR case BER calculation
    
    % Ensure received bits are row vector
    received_bits = received_bits(:)';
    
    % Trim received bits if longer (due to padding)
    if length(received_bits) > length(original_bits)
        received_bits = received_bits(1:length(original_bits));
    end
    
    % Calculate errors
    bit_errors = sum(original_bits ~= received_bits);
    BER = bit_errors / length(original_bits);
end

function plot_BER_vs_SNR(BER_all, SNR_Range, Mod_Types)
    % This function plots BER vs SNR for multiple modulation types
    % Inputs:
    %   BER_all    : matrix (SNR points × modulation types)
    %   SNR_Range  : vector of SNR values in dB
    %   Mod_Types  : cell array of modulation type names (strings)

    % Transpose BER_all if it has the wrong dimensions
    if size(BER_all, 1) ~= length(SNR_Range)
        BER_all = BER_all.';
    end
    
    % Number of modulation types
    num_mods = length(Mod_Types);

    % Define colors and markers for different mod types
    colors = ['b', 'r', 'g', 'k', 'm', 'c', 'y'];
    markers = ['o', 's', '^', 'd', 'x', '+', '*'];

    % Loop over each modulation type and create a new figure for each
    for idx = 1:num_mods
        % Create a new figure for each modulation type
        figure;
        hold on;
        grid on;

        % Plot simulated BER
        semilogy(SNR_Range, BER_all(:, idx), ...
                 [colors(mod(idx-1,length(colors))+1) '-' markers(mod(idx-1,length(markers))+1)], ...
                 'LineWidth', 1.5, 'MarkerSize', 6);

        EbNo = 10.^(SNR_Range/10); % Convert SNR from dB to linear

        % Plot theoretical or tight upper bound BER
        switch Mod_Types{idx}
            case 'BPSK'
                BER_theory = 0.5 * erfc(sqrt(EbNo));
            case 'QPSK'
                BER_theory = 0.5 * erfc(sqrt(EbNo)); % same as BPSK
            case '8PSK'
                BER_theory = erfc(sin(pi/8) * sqrt(3 * EbNo)) / 3;
            case '16-QAM'
                BER_theory = (3/8)*erfc(sqrt((2/5)*EbNo));
            case '64qam'
                BER_theory = (7/24)*erfc(sqrt((7/21)*EbNo));
            case 'BFSK'
                BER_theory = 0.5*erfc(sqrt(0.5*EbNo));
            otherwise
                warning('No theoretical curve for %s. Skipping.', Mod_Types{idx});
                BER_theory = nan(size(EbNo));
        end

        % If theoretical BER is computed, plot it
        if ~any(isnan(BER_theory))
            semilogy(SNR_Range, BER_theory, ...
                     [colors(mod(idx-1,length(colors))+1) '--'], ...
                     'LineWidth', 1.5);
        end

        % Labels and title
        xlabel('E_b/N_0 (dB)');
        ylabel('Bit Error Rate (BER)');
        title(['BER vs. E_b/N_0 for ' Mod_Types{idx}]);

        % Add a legend
        legend_entries = {['Simulated (' Mod_Types{idx} ')'], ['Theoretical (' Mod_Types{idx} ')']};
        legend(legend_entries, 'Location', 'southwest');

        % Set plot limits
        ylim([1e-5 1]);
        xlim([min(SNR_Range) max(SNR_Range)]);
        
        hold off;
    end
end


function displayBitComparison(Tx_bits, Rx_bits, bit_errors, BER, bits_per_group)
    % DISPLAYBITCOMPARISON Display input/output bit comparison and BER results
    %
    % Inputs:
    %   Tx_bits - Transmitted bit sequence
    %   Rx_bits - Received bit sequence
    %   bit_errors - Number of bit errors
    %   BER - Bit Error Rate
    %   bits_per_group - Number of bits to display per row (default: 16)
    
    if nargin < 5
        bits_per_group = 16; % Default to 16-bit groups
    end
    
    % Ensure inputs are row vectors
    Tx_bits = Tx_bits(:)';
    Rx_bits = Rx_bits(:)';
    
    % Display original bits
    fprintf('Original bits:\n');
    disp(reshape(Tx_bits, bits_per_group, [])');
    
    % Display received bits (trimmed to original length)
    fprintf('\nReceived bits:\n');
    disp(reshape(Rx_bits(1:length(Tx_bits)), bits_per_group, [])');
    
    % Display error statistics
    fprintf('\nError Analysis:\n');
    fprintf('Bit errors: %d\n', bit_errors);
    fprintf('BER: %.2e\n', BER);
    
end