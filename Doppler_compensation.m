%% 
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

% Chaotic Sequence Generation
x0 = 0.7;  % Initial condition for logistic map
r = 3.999; % Control parameter for chaotic behavior
chaotic_sequence = zeros(1, length(t));
for i = 2:length(chaotic_sequence)
    x0 = r * x0 * (1 - x0);
    chaotic_sequence(i) = x0;
end
chaotic_sequence = 2 * (chaotic_sequence - 0.5);

% Doppler Effect - Define a frequency range
f_doppler_min = 1000;            % Minimum Doppler frequency shift (transmit)
f_doppler_max = 100000;           % Maximum Doppler frequency shift (receive)

% Generate Doppler-shifted waveforms for the transmitter and receiver side
doppler_waveform_transmit = cos(2 * pi * (f_carrier + f_doppler_min) * t);
doppler_waveform_receive = cos(2 * pi * (f_carrier + f_doppler_max) * t);

% Add the Doppler waveform to the chaotic sequence for combined effect
chaotic_sequence_with_doppler = chaotic_sequence + doppler_waveform_receive;

% Chaotic Spreading
chaotic_spread_signal = bpsk_signal .* chaotic_sequence_with_doppler;

% Receiver Processing
% Despreading
despread_signal = chaotic_spread_signal .* chaotic_sequence_with_doppler;

% Coherent Demodulation
demodulated_signal = despread_signal .* carrier;

% Low-Pass Filter (LPF)
[b, a] = butter(5, 0.2, 'low');
filtered_signal = filtfilt(b, a, demodulated_signal);

% Integrate and Dump
integrated_signal = zeros(1, length(bits));
for i = 1:length(bits)
    bit_segment = filtered_signal((i-1)*n+1:i*n);
    integrated_signal(i) = sum(bit_segment) * dt;
end

% Bit Detection
detected_bits = integrated_signal > 0;

% Plotting
figure;

% Plot Original Chaotic Sequence
subplot(5,1,1);
plot(t, chaotic_sequence, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Original Chaotic Sequence');
axis([0 T -2 2]);

% Plot Doppler Shifted Waveform for Transmitter
subplot(5,1,2);
plot(t, doppler_waveform_transmit, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Doppler Shifted Waveform (Transmitter)');
axis([0 T -2 2]);

% Plot Doppler Shifted Waveform for Receiver
subplot(5,1,3);
plot(t, doppler_waveform_receive, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title(['Doppler Shifted Waveform ' ...
    '(Receiver)']);
axis([0 T -2 2]);

% Plot Chaotic Sequence with Combined Doppler Effect
subplot(5,1,4);
plot(t, chaotic_sequence_with_doppler, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Chaotic Sequence with Combined Doppler Effect (Transmitter + Receiver)');
axis([0 T -4 4]);

% Plot Final Outcome: Detected Bits
subplot(5,1,5);
stairs(0:length(detected_bits)-1, detected_bits, 'LineWidth', 1.5);
xlabel('Bit Index');
ylabel('Detected Bits');
title('Detected Bits After Receiver Processing');
axis([-0.5 length(bits)-0.5 -0.5 1.5]);

% Display results
disp('Original bits:');
disp(bits);
disp('Recovered bits with Doppler Effect and Chaotic Sequence:');
disp(detected_bits);

% Verification of recovery
if isequal(bits, detected_bits)
    disp('Success: Recovered bits match the original bits!');
else
    disp('Error: Recovered bits do not match the original bits.');
    disp('Mismatched positions:');
    disp(find(bits ~= detected_bits));
end
