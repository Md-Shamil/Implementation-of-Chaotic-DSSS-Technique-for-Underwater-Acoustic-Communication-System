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

% Plot Polar NRZ and BPSK Signal
figure;
subplot(2,1,1);
plot(t, polar_nrz, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Polar NRZ Signal');
axis([0 T -2 2]);

subplot(2,1,2);
plot(t, bpsk_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('BPSK Signal');
axis([0 T -2 2]);

%% Chaotic Sequence Generation
x0 = 0.7;  % Initial condition for logistic map
r = 3.999; % Control parameter for chaotic behavior
chaotic_sequence = zeros(1, length(t));
for i = 2:length(chaotic_sequence)
    x0 = r * x0 * (1 - x0);
    chaotic_sequence(i) = x0;
end
chaotic_sequence = 2 * (chaotic_sequence - 0.5);

% Plot Chaotic Sequence
figure;
plot(t, chaotic_sequence, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic Sequence');
axis([0 T -2 2]);

%% Apply Environmental Effects to Chaotic Sequence
doppler_shift=50
% 1. Noise Effect
noise_level = 0.2 + abs(doppler_shift) * 0.05;  % Increase noise as Doppler shift increases
chaotic_sequence_noise = chaotic_sequence + noise_level * randn(size(chaotic_sequence));

% Plot Chaotic Sequence with Noise
figure;
plot(t, chaotic_sequence_noise, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic Sequence with Noise');
axis([0 T -2 2]);

% 2. Multipath Effect
delay_samples = round(0.001 / dt); % 0.001-second delay
multipath_attenuation = 0.5;
multipath_signal = [zeros(1, delay_samples), chaotic_sequence(1:end-delay_samples)];
chaotic_sequence_multipath = chaotic_sequence + multipath_signal * multipath_attenuation;

% Plot Chaotic Sequence with Multipath
figure;
plot(t, chaotic_sequence_multipath, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic Sequence with Multipath');
axis([0 T -2 2]);

% Doppler Shift - Apply to the carrier frequency
doppler_shift = 50; % The Doppler shift frequency (This value will cause incorrect detection beyond a threshold)
f_carrier_doppler_shifted = f_carrier + doppler_shift; % Carrier frequency with Doppler shift
doppler_carrier = cos(2 * pi * f_carrier_doppler_shifted * t);  % Adjusted carrier frequency with Doppler shift

% Modify the chaotic sequence with Doppler-shifted carrier
chaotic_sequence_doppler = chaotic_sequence .* doppler_carrier;  % Multiplying chaotic sequence with Doppler-shifted carrier

% Plot Chaotic Sequence with Doppler Shift
figure;
plot(t, chaotic_sequence_doppler, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title(['Chaotic Sequence with Doppler Shift (Doppler Shift = ', num2str(doppler_shift), ')']);
axis([0 T -2 2]);


% 4. Attenuation
distance_factor = linspace(1, 0.5, length(t));
chaotic_sequence_attenuation = chaotic_sequence .* distance_factor;

% Plot Chaotic Sequence with Attenuation
figure;
plot(t, chaotic_sequence_attenuation, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic Sequence with Attenuation');
axis([0 T -2 2]);
permittivity=1.05;
pressure=0.05;
% Select one or combine the effects for final chaotic sequence
chaotic_sequence_final = chaotic_sequence_noise .* chaotic_sequence_multipath ...
                         .* chaotic_sequence_doppler .* chaotic_sequence_attenuation.*1.05.*14.7;

% Plot Final Chaotic Sequence after Combined Effects
figure;
plot(t, chaotic_sequence_final, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Final Chaotic Sequence with Combined Effects');
axis([0 T -20 20]);

%% Chaotic Spreading
chaotic_spread_signal = bpsk_signal .* chaotic_sequence_final;

% Plot Chaotic Spread Signal
figure;
plot(t, chaotic_spread_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic Spread Signal');
axis([0 T -20 20]);

%% Receiver Processing

% Stage 1: Despreading
despread_signal = chaotic_spread_signal .* chaotic_sequence_final;

% Plot Despread Signal
figure;
plot(t, despread_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Despread Signal');
axis([0 T -100 100]);

% Stage 2: Coherent Demodulation
demodulated_signal = despread_signal .* carrier;

% Plot Demodulated Signal (Continuous)
figure;
plot(t, demodulated_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Demodulated Signal (Continuous)');
axis([0 T -50 50]);

% Stage 3: Low-Pass Filter (LPF)
[b, a] = butter(5, 0.2, 'low');
filtered_signal = filtfilt(b, a, demodulated_signal);

% Plot Filtered Signal
figure;
plot(t, filtered_signal, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Filtered Signal');
axis([0 T -100 100]);

% Stage 4: Integrate and Dump
integrated_signal = zeros(1, length(bits));
for i = 1:length(bits)
    bit_segment = filtered_signal((i-1)*n+1:i*n);
    integrated_signal(i) = sum(bit_segment) * dt;
end

% Plot Integrated Signal
figure;
stem(1:length(bits), integrated_signal, 'LineWidth', 1.5);
xlabel('Bit Index');
ylabel('Amplitude');
title('Integrated Signal');

% Stage 5: Bit Detection
detected_bits = integrated_signal > 0;

% Plot Detected Bits
figure;
stairs(0:length(detected_bits)-1, detected_bits, 'LineWidth', 1.5);
xlabel('Bit Index');
ylabel('Detected Bits');
title('Detected Bits');
axis([-0.5 length(bits)-0.5 -0.5 1.5]);

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
