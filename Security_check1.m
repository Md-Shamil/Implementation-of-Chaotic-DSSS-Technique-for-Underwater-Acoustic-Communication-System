clc;
clear all;
close all;

% Parameters
bitrate = 1;                  % Bitrate for Polar NRZ
bits = [1 0 1 1 0 1 0 1];     % Input binary sequence
n = 100;                      % Number of samples per bit
f_carrier = 5;                % Carrier frequency for BPSK
T = length(bits) / bitrate;   % Total time duration
dt = 1 / (bitrate * n);       % Time resolution
t = 0:dt:T-dt;                % Time vector

% Polar NRZ Encoding
polar_nrz = zeros(1, length(t));
for i = 0:length(bits)-1
    if bits(i+1) == 1
        polar_nrz(i*n+1:(i+1)*n) = 1;   % +1 for bit 1
    else
        polar_nrz(i*n+1:(i+1)*n) = -1;  % -1 for bit 0
    end
end

% BPSK Modulation (without DSSS)
carrier = cos(2 * pi * f_carrier * t);  % Carrier signal
bpsk_signal = polar_nrz .* carrier;     % BPSK modulation (Polar NRZ * Carrier)

% Jamming Signal (white noise added to the signal)
jamming_signal = 0.5 * randn(1, length(t));  % Jamming noise
bpsk_jammed_signal = bpsk_signal + jamming_signal;  % Add jamming noise to BPSK signal

% Plot BPSK Jammed Signal (Without DSSS)
figure;
plot(t, bpsk_jammed_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('BPSK modulated signal with Jamming (Without DSSS)');
axis([0 T -2 2]);

% Chaotic Sequence Generation (using a logistic map for DSSS)
chaotic_sequence = zeros(1, length(t));
x0 = 0.7;  % Initial condition for logistic map
r = 3.999; % Control parameter for chaotic behavior

for i = 2:length(chaotic_sequence)
    x0 = r * x0 * (1 - x0);  % Logistic map
    chaotic_sequence(i) = x0;  % Store chaotic values
end

% Normalize chaotic sequence to be in the range [-1, 1]
chaotic_sequence = 2 * (chaotic_sequence - 0.5);

% Chaotic DSSS (BPSK * Chaotic Sequence)
chaotic_spread_signal = bpsk_signal .* chaotic_sequence;

% Add Jamming Signal to the Chaotic Spreaded Signal
chaotic_jammed_signal = chaotic_spread_signal + jamming_signal;

% Plot Chaotic Spreaded Signal with Jamming
figure;
plot(t, chaotic_jammed_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic DSSS signal with Jamming');
axis([0 T -2 2]);

%% Receiver (With Chaotic DSSS)

% Stage 1: Despreading (multiply by the same chaotic sequence)
despread_signal = chaotic_jammed_signal .* chaotic_sequence;

% Stage 2: Coherent Demodulation (DSSS)
demodulated_signal = despread_signal .* carrier;

% Stage 3: Low-Pass Filter (DSSS)
% Modify the low-pass filter design to make it stronger (e.g., lower cutoff frequency)
[b,a] = butter(7, 0.15, 'low');  % 7th order Butterworth LPF, lower cutoff
filtered_signal = filtfilt(b, a, demodulated_signal);  % Apply zero-phase filtering

% Plot Filtered Signal (With DSSS)
figure;
plot(t, filtered_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Demodulated Signal with DSSS (Filtered)');
axis([0 T -2 2]);

%% Bit Detection

% Integrate and Dump (DSSS)
integrated_signal_dsss = zeros(1, length(bits));
for i = 1:length(bits)
    bit_segment = filtered_signal((i-1)*n+1:i*n);
    integrated_signal_dsss(i) = sum(bit_segment) * dt;
end

% Bit Detection
detected_bits_dsss = integrated_signal_dsss > 0;

% Plot Detected Bits (DSSS)
figure;
stairs(0:length(detected_bits_dsss)-1, detected_bits_dsss, 'LineWidth', 1.5);
xlabel('Bit Index');
ylabel('Detected Bits');
title('Detected Bits with DSSS');
axis([-0.5 length(bits)-0.5 -0.5 1.5]);

% Display Results
disp('Original bits:');
disp(bits);
disp('Detected bits (With DSSS):');
disp(detected_bits_dsss);

% Check if recovered bits match original bits
if isequal(bits, detected_bits_dsss)
    disp('DSSS: Recovered bits match the original bits');
else
    disp('DSSS: Error in bit recovery');
end
