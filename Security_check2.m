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

% Plot Polar NRZ signal (Binary data sequence)
figure;
plot(t, polar_nrz, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Binary Data Sequence');
axis([0 T -2 2]);

% BPSK Modulation
carrier = cos(2 * pi * f_carrier * t);  % Carrier signal
bpsk_signal = polar_nrz .* carrier;     % BPSK modulation (Polar NRZ * Carrier)

% Plot BPSK Modulated Signal
%figure;
%plot(t, bpsk_signal, 'LineWidth', 1.5);
%xlabel('Time');
%ylabel('Amplitude');
%title('BPSK Modulated Signal');
%axis([0 T -2 2]);

% Chaotic Sequence Generation for Transmitter
chaotic_sequence = zeros(1, length(t));
x0 = 0.7;  % Initial condition for logistic map
r = 3.999; % Control parameter for chaotic behavior

for i = 2:length(chaotic_sequence)
    x0 = r * x0 * (1 - x0);  % Logistic map
    chaotic_sequence(i) = x0;  % Store chaotic values
end

% Normalize chaotic sequence to be in the range [-1, 1]
chaotic_sequence = 2 * (chaotic_sequence - 0.5);

% Chaotic Spreading (BPSK * Chaotic Sequence)
chaotic_spread_signal = bpsk_signal .* chaotic_sequence;

% Plot Chaotic Spread Signal
%figure;
%plot(t, chaotic_spread_signal, 'LineWidth', 1.5);
%xlabel('Time');
%ylabel('Amplitude');
%title('Chaotic Spread Signal');
%axis([0 T -2 2]);

%% Receiver

% Stage 1: Despreading (multiply by the same chaotic sequence)
despread_signal = chaotic_spread_signal .* chaotic_sequence;

% Stage 2: Coherent Demodulation
demodulated_signal = despread_signal .* carrier;

% Plot Demodulated Signal (continuous)
figure;
subplot(2,1,1)
plot(t, chaotic_spread_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic Spread Signal');
axis([0 T -2 2]);
subplot(2,1,2);
plot(t, demodulated_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Demodulated Signal (Continuous)');
axis([0 T -2 2]);

% Stage 3: Low-Pass Filter (LPF)
[b,a] = butter(5, 0.2, 'low');  % Design a 5th order Butterworth LPF
filtered_signal = filtfilt(b, a, demodulated_signal);  % Apply zero-phase filtering

% Plot Filtered Signal
figure;

subplot(3,1,1);
plot(t, filtered_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Filtered Signal');
axis([0 T -2 2]);

% Stage 4: Integrate and Dump
integrated_signal = zeros(1, length(bits));
for i = 1:length(bits)
    bit_segment = filtered_signal((i-1)*n+1:i*n);
    integrated_signal(i) = sum(bit_segment) * dt;
end

% Plot Integrated Signal

subplot(3,1,2);
stem(1:length(bits), integrated_signal, 'LineWidth', 1.5);
xlabel('Bit Index');
ylabel('Amplitude');
title('Integrated Signal');

% Stage 5: Bit Detection
detected_bits = integrated_signal > 0;

% Plot Detected Bits

subplot(3,1,3);
stairs(0:length(detected_bits)-1, detected_bits, 'LineWidth', 1.5);
xlabel('Bit Index');
ylabel('Detected Bits');
title('Detected Bits');
axis([-0.5 length(bits)-0.5 -0.5 1.5]);

% Display results
disp('Original bits:');
disp(bits);
disp('Recovered bits:');
disp(detected_bits);

% Check if recovered bits match original bits
if isequal(bits, detected_bits)
    disp('Success: Recovered bits match the original bits!');
else
    disp('Error: Recovered bits do not match the original bits.');
    disp('Mismatched positions:');
    disp(find(bits ~= detected_bits));
end

%% Eavesdropper Simulation

% Generate chaotic sequence for the eavesdropper (randomly)
eavesdropper_chaotic_sequence = 2 * (rand(1, length(t)) - 0.5);

% Eavesdropper attempts to despread and demodulate the signal
eavesdropper_despread_signal = chaotic_spread_signal .* eavesdropper_chaotic_sequence;  % Despreading
eavesdropper_demodulated_signal = eavesdropper_despread_signal .* carrier;  % Coherent Demodulation

% Stage 3: Low-Pass Filter for Eavesdropper
eavesdropper_filtered_signal = filtfilt(b, a, eavesdropper_demodulated_signal);  % Apply zero-phase filtering

% Plot Eavesdropper Filtered Signal
figure;
subplot(3,1,1);
plot(t, eavesdropper_filtered_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Eavesdropper Filtered Signal');
axis([0 T -2 2]);

% Stage 4: Integrate and Dump for Eavesdropper
eavesdropper_integrated_signal = zeros(1, length(bits));
for i = 1:length(bits)
    bit_segment = eavesdropper_filtered_signal((i-1)*n+1:i*n);
    eavesdropper_integrated_signal(i) = sum(bit_segment) * dt;
end

% Plot Eavesdropper Integrated Signal

subplot(3,1,2);
stem(1:length(bits), eavesdropper_integrated_signal, 'LineWidth', 1.5);
xlabel('Bit Index');
ylabel('Amplitude');
title('Eavesdropper Integrated Signal');

% Stage 5: Bit Detection for Eavesdropper
eavesdropper_detected_bits = eavesdropper_integrated_signal > 0;

% Plot Eavesdropper Detected Bits

figure;
stairs(0:length(eavesdropper_detected_bits)-1, eavesdropper_detected_bits, 'LineWidth', 1.5);
xlabel('Bit Index');
ylabel('Detected Bits (Eavesdropper)');
title('Eavesdropper Detected Bits');
axis([-0.5 length(bits)-0.5 -0.5 1.5]);

% Display results for eavesdropper
disp('Eavesdropper detected bits:');
disp(eavesdropper_detected_bits);

% Check if eavesdropper's recovered bits match original bits
if isequal(bits, eavesdropper_detected_bits)
    disp('Eavesdropper Success: Recovered bits match the original bits!');
else
    disp('Eavesdropper Error: Recovered bits do not match the original bits.');
    disp('Mismatched positions:');
    disp(find(bits ~= eavesdropper_detected_bits));
end
