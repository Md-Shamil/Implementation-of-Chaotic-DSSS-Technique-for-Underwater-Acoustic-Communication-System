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

% Carrier Signal
carrier = cos(2 * pi * f_carrier * t);  % Carrier signal
bpsk_signal = polar_nrz .* carrier;     % BPSK modulation (Polar NRZ * Carrier)

% Plot Polar NRZ and BPSK Signal with PSDs
figure;
subplot(3,2,1);
plot(t, polar_nrz, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Polar NRZ Signal');
axis([0 T -2 2]);

subplot(3,2,2);
pwelch(polar_nrz, [], [], [], 1/dt);
title('PSD of Polar NRZ Signal');

subplot(3,2,3);
plot(t, bpsk_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('BPSK Signal');
axis([0 T -2 2]);

subplot(3,2,4);
pwelch(bpsk_signal, [], [], [], 1/dt);
title('PSD of BPSK Signal');

% Chaotic Sequence Generation
x0 = 0.7;  % Initial condition for logistic map
r = 3.999; % Control parameter for chaotic behavior
chaotic_sequence = zeros(1, length(t));
for i = 2:length(chaotic_sequence)
    x0 = r * x0 * (1 - x0);
    chaotic_sequence(i) = x0;
end
chaotic_sequence = 2 * (chaotic_sequence - 0.5);

% Plot Chaotic Sequence with PSD
figure;
subplot(3,2,1);
plot(t, chaotic_sequence, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic Sequence');
axis([0 T -2 2]);

subplot(3,2,2);
pwelch(chaotic_sequence, [], [], [], 1/dt);
title('PSD of Chaotic Sequence');

% Apply Environmental Effects to Chaotic Sequence

% 1. Noise Effect
noise_level = 0.2;
chaotic_sequence_noise = chaotic_sequence + noise_level * randn(size(chaotic_sequence));

% Plot Chaotic Sequence with Noise and PSD
subplot(3,2,3);
plot(t, chaotic_sequence_noise, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic Sequence with Noise');
axis([0 T -2 2]);

subplot(3,2,4);
pwelch(chaotic_sequence_noise, [], [], [], 1/dt);
title('PSD of Chaotic Sequence with Noise');

% 2. Multipath Effect
delay_samples = round(0.001 / dt); % 0.001-second delay
multipath_attenuation = 0.5;
multipath_signal = [zeros(1, delay_samples), chaotic_sequence(1:end-delay_samples)];
chaotic_sequence_multipath = chaotic_sequence + multipath_signal * multipath_attenuation;

% Plot Chaotic Sequence with Multipath and PSD
subplot(3,2,5);
plot(t, chaotic_sequence_multipath, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic Sequence with Multipath');
axis([0 T -2 2]);

subplot(3,2,6);
pwelch(chaotic_sequence_multipath, [], [], [], 1/dt);
title('PSD of Chaotic Sequence with Multipath');

% Doppler Shift
doppler_shift = 0.5;
doppler_carrier = cos(2 * pi * (f_carrier + doppler_shift) * t);
chaotic_sequence_doppler = chaotic_sequence .* doppler_carrier;

% Attenuation
distance_factor = linspace(1, 0.5, length(t));
chaotic_sequence_attenuation = chaotic_sequence .* distance_factor;

% Final Chaotic Sequence with Combined Effects
chaotic_sequence_final = chaotic_sequence_noise .* chaotic_sequence_multipath ...
                         .* chaotic_sequence_doppler .* chaotic_sequence_attenuation;

% Plot Final Chaotic Sequence with Combined Effects and PSD
figure;
subplot(2,1,1);
plot(t, chaotic_sequence_final, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Final Chaotic Sequence with Combined Effects');
axis([0 T -2 2]);

subplot(2,1,2);
pwelch(chaotic_sequence_final, [], [], [], 1/dt);
title('PSD of Final Chaotic Sequence');

% Chaotic Spreading
chaotic_spread_signal = bpsk_signal .* chaotic_sequence_final;

% Plot Chaotic Spread Signal and PSD
figure;
subplot(2,1,1);
plot(t, chaotic_spread_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic Spread Signal');
axis([0 T -2 2]);

subplot(2,1,2);
pwelch(chaotic_spread_signal, [], [], [], 1/dt);
title('PSD of Chaotic Spread Signal');

% Receiver Processing

% Stage 1: Despreading
despread_signal = chaotic_spread_signal .* chaotic_sequence_final;

% Stage 2: Coherent Demodulation
demodulated_signal = despread_signal .* carrier;

% Stage 3: Low-Pass Filter (LPF)
[b, a] = butter(5, 0.2, 'low');
filtered_signal = filtfilt(b, a, demodulated_signal);

% Stage 4: Integrate and Dump
integrated_signal = zeros(1, length(bits));
for i = 1:length(bits)
    bit_segment = filtered_signal((i-1)*n+1:i*n);
    integrated_signal(i) = sum(bit_segment) * dt;
end

% Bit Detection
detected_bits = integrated_signal > 0;

% Display results
disp('Original bits:');
disp(bits);
disp('Recovered bits with Environmental Effects on Chaotic Sequence:');
disp(detected_bits);

% Verification of recovery
if isequal(bits, detected_bits)
    disp('Success: Recovered bits match the original bits!');
else
    disp('Error: Recovered bits do not match the original bits.');
    disp('Mismatched positions:');
    disp(find(bits ~= detected_bits));
end
