%% Load DICOM data with file named "filename"
filename = '';
tic
data = dicomread('filename');
toc
readFileTime = toc


%% Getting the useful images from the whole image sequence.

% Total number of saved dicom image
total_frame = size(data,4);

% Extract background information
full_bg = rgb2gray(data(:,:,:,1));

% Boundary of useful signal
signal_y = 60:560; % range of y index includes useful signal
signal_x = 1:880; % range of x index includes useful signal

full_bg(y,x) = 0; %% background of image inclusing all the text
cropped_bg = full_bg(signal_y,signal_x);
bg_noise = mean(mean(cropped_bg));

% I use bg_noise = 40 initially
bg_noise = 40;

% Number of frames processed by this program
% use processed_frame = total_frame for full processing.
processed_frame = total_frame;

% Useful signal boundary size.
size_y = length(signal_y);
size_x = length(signal_x);

% initial_signal = zeros(size(full_bg,1),size(full_bg,2),processed_frame);
% cropped_signal = zeros(size_y,size_x,processed_frame);

curr_index = 1;
for j = 1:processed_frame
    if (curr_index > total_frame)
        processed_frame = j-1;
        break;
    end
    % signal_initial is D mode signal with background information
    % initial RGB image in DICOM file changed into greyscale
    signal_initial(:,:,j) = rgb2gray(data(:,:,:,curr_index));
    % curr_index is the index number of frame with no overlapping signal
    curr_index = floor(curr_index + 221.5);
    % signal_cropped is D mode signal without background information    
    signal_cropped(:,:,j) = signal_initial(signal_y,signal_x,j) - (cropped_bg);
end
signal_cropped(signal_cropped < bg_noise) = 0;
signal_cropped(signal_cropped >= bg_noise) = 1;

%% Acquire the envolope of ultrasound signal from DICOM image. 
% The envolope is the magnitude of the signal and velocity of trachea
% vibration, in cm/s.

signal_combined = signal_cropped(:,:,1);
for j = 2:processed_frame
    signal_combined = horzcat(signal_combined,signal_cropped(:,:,j));
end
se = strel('disk',6);
signal_closed = imclose(signal_combined,se);

% envolope of the signal
signal_amplitude = sum(signal_closed,1)';
signal_amplitude = signal_amplitude/size_y * 30;

%% Using the peaks of signal amplitude to count inspiration cycle
% The assumeption is that inspiration peaks is higher than expiration peaks
%
% NEED MODIFICATION SINCE THE ASSUMPTION IS NOT ALWAYS CORRECT.

[PKS,LOCS,W,P] = findpeaks(signal_amplitude, 'MinPeakHeight', ...
    0.5*max(signal_amplitude), 'MinPeakProminence',0.2*max(signal_amplitude));

% inspiration_dicom is the inspiration count found through DICOM data
inspiration_dicom = length(PKS);

%% Load and process spirometer data.
spirometer_filename = '';
spirometer = importdata(spirometer_filename);
% Load spirometer volume data
volume = spirometer.data(:,1);
% Since there is shift in the raw data, detrend is performed to remove it.
volume_detrend = detrend(volume);
% Change the baseline of detrended volume to 0 since detrend() changed it.
volume_detrend = volume_detrend - max(volume_detrend);
% Inverse the sign of volume so it's positive when inspiration ends so that
% volumes_inverse are all non-negative numbers. (For easier correlation
% later with integration of ultrasound signal).
volume_inverse = - volume_detrend; 

% Load spirometer flow data;
flow = spirometer.data(:,2);
% Inverse the sign of flow so positive is inspiration.
flow_inverse = -(flow); 
% Absolute value of low, since we need compare it with ultrasound data.
% Only magnitude matters when compare the signal amplitude later.
flow_abs = abs(flow);

%% Find peaks in spirometer volume and flow.
% peaks in volume_detrend means the end/begning of one
% inspiration-expiration cycle. 
[PKS_v,LOCS_v,W_v,P_v] = findpeaks(volume_detrend, ...
    'MinPeakProminence',0.1*max(volume_inverse));

% inspiration_volume = number of counts of inspiration/expiration cycles
% using volume data.
inspiration_volume = length(PKS_v);

% find the peaks in flow and use the inspiration peaks to count the cycles
% of inspiration/expieration
[PKS_f,LOCS_f,W_f,P_f] = findpeaks(flow_inverse, 'MinPeakHeight',...
    0.5*max(flow_inverse),'MinPeakProminence',0.5*max(flow_inverse));
inspiration_flow = length(PKS_f);

%% Resample data in ultrasound signal to match the points of spirometer 

size_spirometer = size(volume,1);
size_dicom = size(signal_amplitude,1);

index_dicom = linspace(1,size_dicom,size_spirometer);
signal_resampled = zeros(size_spirometer,1);

% resample the signal, use 1st order interpolation
for i = 1:size_spirometer
    signal_resampled(i) = signal_amplitude(floor(index_dicom(i)));
    if index_dicom(i) > floor(index_dicom(i))
        signal_resampled(i) = (index_dicom(i)-floor(index_dicom(i))) * ...
           (signal_resampled(i+1)-signal_resampled(i))+signal_resampled(i);
    end
end
    

%% Calculate Correlation between ultrasound signal and spirometer flow
[Pearson_RHO_Flow,Pearson_PVAL_Flow] = ...
    corr(signal_resampled,flow_abs,'Type','Pearson')
[Spearman_RHO_Flow,Spearman_PVAL_Flow] = ...
    corr(signal_resampled,flow_abs,'Type','Spearman')

signal_inspiration = signal_resampled(flow_match>=0);
signal_expiration = signal_resampled(flow_match<=0);
flow_inspiration = flow_abs(flow_match>=0);
flow_expiration = flow_abs(flow_match<=0);

[Pearson_RHO_Flow_in,Pearson_PVAL_Flow_in] = ...
    corr(signal_inspiration,flow_inspiration,'Type','Pearson')
[Spearman_RHO_Flow_in,Spearman_PVAL_Flow_in] = ...
    corr(signal_inspiration,flow_inspiration,'Type','Spearman')

[Pearson_RHO_Flow_ex,Pearson_PVAL_Flow_ex] = ...
    corr(signal_expiration,flow_expiration,'Type','Pearson')
[Spearman_RHO_Flow_ex,Spearman_PVAL_Flow_ex] = ...
    corr(signal_expiration,flow_expiration,'Type','Spearman')

%% Calculate correlation between ultrasound signal and spirometer volume
matched_auc_in = zeros(size(signal_resampled));
matched_auc_ex = zeros(size(signal_resampled));
index = 1;
for i = 1 : length(flow_match)
    if (flow_match(i) < 0)
        matched_auc_in(i) = 0;
        index = i;
    else
        matched_auc_in(i) = 0.01*sum(signal_resampled(index:i));
    end
end
index = length(flow_match);
for i = length(flow_match):-1:1
    if (flow_match(i) < 0)
        matched_auc_ex(i) = 0.01*sum(signal_resampled(i:index));
    else
        matched_auc_ex(i)= 0;
        index = i;
    end
end

%% Display correlation result in a table

T = table({'Inspiration_Flow'; 'Expiration_Flow';'Overall_FLow';...
    'Inspiration_Volume'; 'Expiration_Volume'}, ...
    [Pearson_RHO_Flow_in; Pearson_RHO_Flow_ex; Pearson_RHO_Flow;...
    Pearson_RHO_Volume_in;Pearson_RHO_Volume_ex], ...
    [Pearson_PVAL_Flow_in; Pearson_PVAL_Flow_ex; Pearson_PVAL_Flow; ...
    Pearson_PVAL_Volume_in; Pearson_PVAL_Volume_ex], ...
    [Spearman_RHO_Flow_in; Spearman_RHO_Flow_ex; Spearman_RHO_Flow;...
    Spearman_RHO_Volume_in;Spearman_RHO_Volume_ex], ...
    [Spearman_PVAL_Flow_in; Spearman_PVAL_Flow_ex; Spearman_PVAL_Flow; ...
    Spearman_PVAL_Volume_in; Spearman_PVAL_Volume_ex]);

T.Properties.VariableNames = {'Type' 'Pearson_Corr' 'P_p_value' ...
    'Spearman_Corr' 'S_p_value'};
      
disp(T)
      
      
      