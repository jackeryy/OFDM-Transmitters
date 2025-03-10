% CMPEN 462 Mini Project 1 - Basic OFDM Transmitters
% Jack O'Brien

% Generate message and a bit stream
message = 'WirelessCommunicationSystemsandSecurityJackHenryOBrien';
bitstream = reshape(dec2bin(uint8(message), 8).', 1, []) - '0';
bitstream_vector = bitstream(:);

%OFDM parameters
N_fft = 4096;      %FFT size
subcarrier_bw = 60e3;   % (Hz)
slot_time = 0.25e-3;    % 25 milliseconds
fs = subcarrier_bw * N_fft;       %Sample rate (Hz)

%Determine Total bits for one time slot 64QAM
Tc = 1/ fs;     %Fundamental time unit
T_symbol = N_fft * Tc;      %OFDM symbol duration
symbols = slot_time/T_symbol;      %OFDM symbols per slot

% Bits per symbol
bits = 6; 
bits_needed = ceil(symbols * bits * N_fft);

% CP cycles
cp_normal =  2304; % (From Table 7.12)

% Expand bitstream for modulation
bitstream = repmat(bitstream,1,ceil(bits_needed/length(bitstream)));
bitstream = bitstream(1:bits_needed);


% BPSK and Pi/2 BPSK Mapping
BPSK = ((1 - 2 * bitstream) + 1j * (1 - 2 * bitstream))/sqrt(2);
phaseshift = exp(1j * (pi/2) * mod(0:length(bitstream)-1,2));
pi2_BPSK = (phaseshift .* ((1 - 2 * bitstream) + 1j * (1 - 2 * bitstream)))/sqrt(2); 

% QPSK Maps ( 2 bits per symbol )
QPSK_bits = reshape(bitstream, 2, []).';
QPSK_map = [1 + 1j, -1 + 1j, 1 - 1j, -1 - 1j] / sqrt(2);   

%change 2-bit pairs into decimals
QPSK_indices = QPSK_bits(:,1) * 2 + QPSK_bits(:,2) + 1;   
QPSK = QPSK_map(QPSK_indices);

% 64QAM ( 6 bits per symbol ) - gray mapping
QAM_bits = reshape(bitstream, 6, []).';
I_comp = (1 - 2 * QAM_bits(:,1)) .* (4 - (1 - 2 * QAM_bits(:,3)) .* (2 - (1 - 2 * QAM_bits(:,5))));
Q_comp = (1 - 2 * QAM_bits(:,2)) .* (4 - (1 - 2 * QAM_bits(:,4)) .* (2 - (1 - 2 * QAM_bits(:,6))));
QAM64 = (I_comp + 1j * Q_comp) / sqrt(42);


% Serial to parallel conversion
BPSK_parallel = reshape(BPSK, N_fft, []);
pi2_BPSK_parallel = reshape(pi2_BPSK, N_fft, []);
QPSK_parallel = reshape(QPSK, N_fft, []);
QAM_parallel = reshape(QAM64, N_fft, []);



% Perform 4096-point IFFT
BPSK_ofdm = ifft(BPSK_parallel, N_fft);
pi2_BPSK_ofdm = ifft(pi2_BPSK_parallel, N_fft);
QPSK_ofdm = ifft(QPSK_parallel, N_fft);
QAM_ofdm = ifft(QAM_parallel, N_fft);


% Add cyclic prefix before modulations
BPSK_ofdm_cp = [BPSK_ofdm(end-cp_normal+1:end, :); BPSK_ofdm];
pi2_BPSK_ofdm_cp = [pi2_BPSK_ofdm(end-cp_normal+1:end, :); pi2_BPSK_ofdm];
QPSK_ofdm_cp = [QPSK_ofdm(end-cp_normal+1:end, :); QPSK_ofdm];
QAM_ofdm_cp = [QAM_ofdm(end-cp_normal+1:end, :); QAM_ofdm];


% Convert OFDM symbols to 1D vectors
pi2_BPSK_vector = pi2_BPSK_ofdm_cp(:);
BPSK_vector = BPSK_ofdm_cp(:);
QPSK_vector = QPSK_ofdm_cp(:);
QAM64_vector = QAM_ofdm_cp(:);


% Number of samples for two OFDM symbols 
num_samples = 2 * cp_normal;  

% Time vector for two OFDM symbols in microseconds
t = (0:num_samples-1) * Tc * 1e6; 

% Plot two OFDM symbols for BPSK
figure(1);
plot(t, (BPSK_ofdm_cp(1:num_samples, 1)), 'b');
xlabel('Time (µs)');
ylabel('Amplitude');
title('BPSK - Two OFDM Symbols');
grid on;

% Plot two OFDM symbols for Pi/2 BPSK
figure(2);
plot(t, (pi2_BPSK_ofdm_cp(1:num_samples, 1)), 'r');
xlabel('Time (µs)');
ylabel('Amplitude');
title('Pi/2 BPSK - Two OFDM Symbols');
grid on;

% Plot two OFDM symbols for QPSK
figure(3);
plot(t, (QPSK_ofdm_cp(1:num_samples, 1)), 'g');
xlabel('Time (µs)');
ylabel('Amplitude');
title('QPSK - Two OFDM Symbols');
grid on;

% Plot two OFDM symbols for 64QAM
figure(4);
plot(t,(QAM_ofdm_cp(1:num_samples, 1)), 'm');
xlabel('Time (µs)');
ylabel('Amplitude');
title('64QAM - Two OFDM Symbols');
grid on;
