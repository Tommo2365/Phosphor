clear all;
close all;

% Load Data
filename = ['G:\.shortcut-targets-by-id\1B8SEBzrQvBxWdQRqdwpsGrI0tFR0vHlg\' ...
    'Tom\100 Phosphor Thermometry\102 Data\102 Luminescence Decay Time Thermometry\101 Work Area\' ...
    '101 Proof of concept\101 Decay Time Cal Data\120' ...
    'C.csv'];

n_rows = 300;
raw_data = readmatrix(filename, 'NumHeaderLines', n_rows); % ignore first n row

% Given Data
x = raw_data(:,1);
y = raw_data(:,2);

% Detect edges based on predefined indices
edge_indices = [100, 500, 900, 1300, 1700];

% Split waves into segments
waveform_length = 40; % chop data at specific indices

% time dependent x-axis
x_time_inc = abs(raw_data(2,1) - raw_data(1,1)); % Time incriments
time_vector = linspace(0, x_time_inc, (waveform_length));

% Segment Data into Waveforms
app_waveforms = [];
for i = 1:length(edge_indices)
    if edge_indices(i) + waveform_length - 1 <= length(y)
        single_waveform = y(edge_indices(i):edge_indices(i)+ waveform_length-1);
        app_waveforms(:,i) = single_waveform;
    end
end

% Average and normalize the waveforms
average_waveform = mean(app_waveforms, 2);

% Apply exponential least squares fit
f = fit(time_vector', average_waveform, 'exp1');

% Plot the data and the fit
figure;
plot(time_vector, average_waveform, 'o'); % Plot the data as circles
hold on;
plot(f, time_vector, average_waveform); % Plot the fit
xlabel('Time');
ylabel('Average Waveform');
%legend('Data', 'Exponential Fit');
title('Least Squares Exponential Fit');
hold on;


% Corrected Successive Integration (CSI) method - y(t) = a*int(yt)dt + b*t + c
integrated_waveform = cumsum(average_waveform) * x_time_inc; % computes the cumulative sum of the elements of a vector. For discrete data, this is analogous to summing up small "areas" (or values) as you move along the vector. When you have data sampled at regular intervals, the cumulative sum approximates the integral of the data.

% y(t) = average_waveform 
% int(y(t)) = intergrated_average_waveform
% t = time_vector

% Linear regression of the sum
design_matrix = [integrated_waveform, time_vector', ones(length(time_vector), 1)]; % This gives the least squares solution for the coefficients of the linear model


coeff = design_matrix \ average_waveform; % Solve for coefficients

coeff_a = coeff(1); % coefficient of the inergrated waveform, int(y(t))
coeff_b = coeff(2); % coefficient for the time vector and is used to determine the decay constant.
coeff_c = coeff(3); % the constant term or y-intercept, c

tau = -1/coeff_b; % Extract decay constant
fprintf('the decay constant, b = Tau = %f\n',  tau)

plot(time_vector, average_waveform);
title(['Decay constant (lambda) = ', num2str(tau)]);
legend('Data', 'Exponential Fit', 'CSI Method')

% Fourier Transform method
%F = fft(average_waveform); % Compute the Fourier transform of the waveform
%omega = 2 * pi * (0:length(time_vector)-1) / (time_vector(end) - time_vector(1)); % Frequency vector

% Extract decay constant from Fourier transform
%lambda_FT = real(F(2:end)) ./ imag(F(2:end)); % Exclude the DC component (index 1)

% The decay constant is typically extracted from the first non-zero frequency component
%lambda_from_FT = lambda_FT(1);
%plot(time_vector, average_waveform);
%title(['Decay constant (lambda) from Fourier Transform = ', num2str(lambda_from_FT)]);

% Trapezoidal Integration
integrated_waveform_trapz = cumtrapz(time_vector, average_waveform) * x_time_inc;

% Linear regression using trapezoidal integration
X_trapz = [integrated_waveform_trapz, time_vector', ones(length(time_vector), 1)];
coeff_trapz = X_trapz \ average_waveform;

% Extract the initial approximate decay constant
lambda_initial = coeff_trapz(2);

% Apply the correction formula
lambda_corrected = 1/ (log(lambda_initial^2 + 1) / (lambda_initial^2 - 1));

% Display the corrected decay constant
disp(['Corrected Decay Constant (lambda) = ', num2str(lambda_corrected)]);
