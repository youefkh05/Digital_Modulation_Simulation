%close all;
clear;
N_Bits = 120000;
% Modulation Schemes Tables
M_BPSK_Table = [-1, 1];
M_QPSK_Grey_Table = [-1-1i, -1+1i, 1-1i, 1+1i];
M_8PSK_Table = [
1, 1/sqrt(2)+1i/sqrt(2), -1/sqrt(2)+1i/sqrt(2), 1j,...
1/sqrt(2)-1i/sqrt(2), -1i, -1, -1/sqrt(2)-1i/sqrt(2)
];
M_16QAM_Table = [
-3-3j, -3-1j, -3+3j, -3+1j,...
-1-3j, -1-1j, -1+3j, -1+1j,...
3-3j, 3-1j, 3+3j, 3+1j,...
1-3j, 1-1j, 1+3j, 1+1j
];
M_BFSK_Table = [1, 1j];
M_QPSK_Not_Grey_Table = [-1-1i, -1+1i, 1+1i, 1-1i];
% Random Bit Stream
rawMessage = randi([0 1], [1 N_Bits]);
% Define SNR range
EBoverN0_dB = -4:1:12;
EBoverN0_linear = 10 .^ (EBoverN0_dB / 10);
% Cell Array carrying tables for modulation schemes
modulationSchemes = {
M_BPSK_Table, M_QPSK_Grey_Table, M_8PSK_Table,...
M_16QAM_Table, M_BFSK_Table, M_QPSK_Not_Grey_Table
};
modulationSchemesNames = {'BPSK', 'QPSK (Grey)', '8PSK', '16QAM', 'BFSK', 'QPSK (Not Grey)'};
BER_theoritical = {
0.5*erfc(sqrt(EBoverN0_linear)),...
0.5*erfc(sqrt(EBoverN0_linear)),...
erfc(sin(pi/8) * sqrt(3 * EBoverN0_linear)) / 3,...0.375*erfc(sqrt(0.4*EBoverN0_linear)),...
0.5*erfc(sqrt(0.5*EBoverN0_linear)),...
0.5*erfc(sqrt(EBoverN0_linear))
};
BER_actual = cell(1, length(modulationSchemes));
% The BER simulation for all schemes goes in a loop
for i = 1:length(modulationSchemes)
modulationScheme = modulationSchemes{i}; % Choose Scheme
M_Name = modulationSchemesNames{i}; % Get scheme name
bitsPerSymbol = log2(length(modulationScheme)); % Obtain number of bits per symbol
Eavg = mean(abs(modulationScheme) .^ 2); % Obtain Average Energy of Constellation
N0 = Eavg ./ (bitsPerSymbol * EBoverN0_linear); % Adjust Noise according to SNR and Eavg
% Generate Baseband Equivalent Signal from Raw Data Using Modulation Scheme
transmittedMessage = mapper(rawMessage, modulationScheme);
% Initialize BER vector
BER_actual{i} = zeros(1, length(EBoverN0_dB));
for m = 1:length(EBoverN0_dB)
% Create Channel Noise
NI = sqrt(N0(m) / 2) * randn(1, N_Bits / bitsPerSymbol);
NQ = sqrt(N0(m) / 2) * randn(1, N_Bits / bitsPerSymbol);
noise = NI + 1i * NQ;
receivedMessage = transmittedMessage + noise; % Add Channel Noise
% Decision Block
estimatedMessage = demapper(receivedMessage, modulationScheme);
BER_actual{i}(m) = length(find(rawMessage ~= estimatedMessage));
end
BER_actual{i} = BER_actual{i} / N_Bits;
% Plot Theoritical vs Simulated Overlaid
figure;
semilogy(EBoverN0_dB, BER_actual{i});
hold on;
semilogy(EBoverN0_dB, BER_theoritical{i});
grid on;
xlabel('E_b/N_o (dB)'); ylabel('BER');
title(['BER of ' M_Name]);
legend({['Simulated BER of ' M_Name], ['Theoritical BER of ' M_Name]});
ylim([10e-5 1]);
end
% Plot BPSK, QPSK, 8PSK, and 16-QAM overlaid
figure;
for i = 1:4
semilogy(EBoverN0_dB, BER_actual{i});
hold on;
end
xlabel('E_b/N_o (dB)'); ylabel('BER');
title('Simulated BER of the four modulation schemes');
ylim([10e-5 1]);
grid on; hold off;
legend(modulationSchemesNames);
% Plot Grey vs Non-Grey Encoded QPSK
figure;
semilogy(EBoverN0_dB, BER_actual{2});
hold on;
semilogy(EBoverN0_dB, BER_actual{6});
xlabel('E_b/N_o (dB)'); ylabel('BER');
title('Simulated BER of Grey Encoded vs Non-Grey Encoded QPSK');
ylim([10e-5 1]);
grid on; hold off;
legend({modulationSchemesNames{2}, modulationSchemesNames{6}});

Eb = 1;
Tb = 1;
N_Bits = 100;
N_Samples = 7;
N_realizations = 500;
Fc = 1/Tb;
F1 = 1/Tb;
F2 = 2/Tb;
iteration = 1000;
PLR_S_cor = zeros(1,N_Bits*N_Samples);
t = 0:1/N_Samples:(1-1/N_Samples);
S1 = sqrt(2*Eb/Tb)*cos(2*pi*(F1-Fc)*t)+1i*sqrt(2*Eb/Tb)*sin(2*pi*(F1-Fc)*t);
S2 = sqrt(2*Eb/Tb)*cos(2*pi*(F2-Fc)*t)+1i*sqrt(2*Eb/Tb)*sin(2*pi*(F2-Fc)*t);

for k=1:1:iteration
initial_start = randi([0,N_Samples-1],N_realizations,1);
initial_ensemble = repelem(randi([0,1],N_realizations,N_Bits),1,N_Samples);
PLR_ensemble = zeros(N_realizations,N_Bits*N_Samples);
for i=1:1:N_realizations
for j=1:N_Samples:N_Bits*N_Sample
    if(initial_ensemble(i,j)==0)
    initial_ensemble(i,j:j+N_Samples-1) = S1;
    else
    initial_ensemble(i,j:j+N_Samples-1) = S2;
    end
end
end
for i=1:500
PLR_ensemble(i,:) = circshift(initial_ensemble(i,:),initial_start(i));
end
for t=1:(N_Bits*N_Samples)
PLR_S_cor(t) = PLR_S_cor(t)+sum(conj(PLR_ensemble(:,N_Bits*N_Samples/2))...
.*PLR_ensemble(:,t))/N_realizations;
end
end

Fs = N_Samples;
n = N_Bits*N_Samples;
f=(-n/2:n/2-1)*Fs/n;
figure;
plot(f,abs(fftshift(fft(PLR_S_cor)))/(Fs*iteration),'linewidth',1);
hold on;
plotTheoriticalBFSK(f, Eb, Tb);
title('PSD of BFSK','fontsize',18);
annotation('arrow',[0.44 0.44],[0.11 0.92],'color','r','linewidth',3);
annotation('arrow',[0.595 0.595],[0.11 0.92],'color','r','linewidth',3);
xlabel('f','fontsize',18);
ylabel('S_x(f)','fontsize',18);
legend('actual PSD','Theoretical PSD');
xlim([-2 3]);
ylim([0 2]);

function plotTheoriticalBFSK(f, Eb, Tb)
phi = 8*Eb*cos(pi*Tb*(f-0.5/Tb)).^2 ./ (pi^2 * (4*Tb^2*(f-0.5/Tb).^2 - 1).^2);
phi(find(phi == Inf)) = 0;
PSD = phi + diracDelta(length(f), find(f == 0)) + diracDelta(length(f), find(f== 1));
plot(f, PSD);
end

function delta = diracDelta(vectorSize, deltaPosition, deltaValue)
if (~exist('deltaValue', 'var'))
    deltaValue = 100;
end
delta = [zeros(1, deltaPosition - 1) deltaValue zeros(1, vectorSize -deltaPosition)];
end


function transmittedMessage = mapper(rawMessage, modulationScheme)
N_Bits = length(rawMessage);
bitsPerSymbol = log2(length(modulationScheme));
transmittedMessage = zeros(1, N_Bits / bitsPerSymbol);
for j = 1:N_Bits/bitsPerSymbol
startBit = bitsPerSymbol * (j - 1) + 1;
endBit = bitsPerSymbol * j;
% An efficient binary to decimal conversion
index = bin2dec(string(char(rawMessage(startBit:endBit) + '0'))) + 1;
transmittedMessage(j) = modulationScheme(index);
end
end

function estimatedMessage = demapper(receivedMessage, modulationScheme)
bitsPerSymbol = log2(length(modulationScheme));
N_Bits = bitsPerSymbol * length(receivedMessage);
estimatedMessage = zeros(N_Bits / bitsPerSymbol, bitsPerSymbol);
for k = 1:N_Bits/bitsPerSymbol
% Apply maximum likelihood rule
% Find the symbol with minimum Euclidean distance
[~, index] = min(abs(receivedMessage(k) - modulationScheme));
estimatedBits = myDecimalToBinaryVector(index - 1, bitsPerSymbol);
estimatedMessage(k,:) = estimatedBits;
end
estimatedMessage = reshape(estimatedMessage', [1, N_Bits]);
end

function binaryVector = myDecimalToBinaryVector(decimal, numBits)
% An efficient and optimized decimal to binary conversion
binaryVector = dec2bin(decimal) - '0';
binaryVector = [zeros(1, numBits - length(binaryVector)) binaryVector];
end
