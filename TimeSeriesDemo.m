% %imports data from the sound file 18apr25_2.wav as variable mic_dat
% %note fields of mic_dat are 'data' corresponding to the waveform and 'fs'
% %corresponding to the sampling frequency
% mic_dat = importdata('Data\18apr25_2.wav');
% % create variable 'micSig' for the waveform
% micSig = mic_dat.data;
% % create variable 'fs' for the sampling rate
% fs = mic_dat.fs;
% 
% % create time base (T)
% T = [1:size(micSig,1)]/fs;
% % create figure
% figure;
% % plot micSig as a function of time
% plot (T,micSig(:,1));
% xlabel('Time (s)');
% ylabel('Amplitude (V)');
% title('Microphone Signal');
% 
% %rescale over a 100 ms window between 5 and 5.1 seconds
% axis([5, 5.1, -.2 .2]);
% 
% %plot individual data points
% hold on
% plot(T,micSig(:,1),'rx')
% hold off
% 
% %QUESTIONS:
% % 1) Draw the waveform of the microphone signal
% % 2) What is the frequency of the waveform?
% % 3) How many data points are in each waveform?
% 
% %%
% %periodogram - amplitude vs frequency (i.e. frequency domain) MATLAB
% %DEFAULTS
% 
% figure
% [pxx, f] = periodogram(micSig(:,1),[],[],fs);
% plot(f, pxx);
% xlabel('Frequency (Hz)')
% ylabel('Power (V^2)')
% title('Power Spectrum (Periodogram) of Mic Signal');
% % conversion to decibles 10 log (x)
% plot(f, 10*log10(pxx))
% xlabel('Frequency (Hz)')
% ylabel('dB Power (V^2)')
% title('Power Spectrum (Periodogram) of Mic Signal');
% %QUESTIONS:
% % 1) What is the highest frequency signal that can be analyzed if the
% % sampling rate is 22050 Hz?
% % 2) What is the lowest frequency signal that can be measured for this
% % recording?
% 
% %rescale periodogram
% axis([0 1000 -120 -30]);
% 
% %QUESTIONS:
% % 1) What is the dominant frequency on the periodogram?
% % 2) What biological/mechanical process does this represent?
% 
% %periodogram - amplitude vs frequency (i.e. frequency domain) Windowed to
% %our needs
% % Want to calculate power @ each frequency in range 60 Hz - 400 Hz in 2 Hz
% % intervals
% hold all
% [pxx2, f2] = periodogram(micSig(:,1),[],[60:2:400],fs);
% plot(f2,10* log10(pxx2),'r-');
% 
% 
% %% Filtering (part 1 - RC circuits)
% %default filter, 3dB cut-off @ 300 Hz
% LP1=designfilt('lowpassiir', 'FilterOrder', 1, 'HalfPowerFrequency', 300, 'SampleRate', 22050);
% %uncomment the line below to design filter LP1
% %LP1 = designfilt;
% LoPassSig = filter(LP1,micSig(:,1));
% figure;
% plot(T,micSig(:,1))
% hold on
% plot(T,LoPassSig,'r')
% xlabel ('Time (s)')
% ylabel ('Amplitude (V)')
% title ('effect of a LP filter cut-off @ 300 Hz')
% legend([{'Original'};{'300 Hz Low-Pass'}]);
% axis([5, 5.1, -.2 .2]);
% audiowrite('Data\LoPass.wav',LoPassSig,fs);
% % QUESTIONS
% % 1) What components of the original signal remain? What components are
% % attenuated?
% % 2) Why is there a right-shift in the signal?
% HP1=designfilt('highpassiir', 'FilterOrder', 1, 'HalfPowerFrequency', 300, 'SampleRate', 22050);
% HiPassSig = filter(HP1,micSig(:,1));
% plot(T,HiPassSig,'g');
% audiowrite('Data\HiPass.wav',HiPassSig,fs);
% legend([{'Original'};{'300 Hz Low-Pass'};{'300 Hz High-Pass'}]);
% % QUESTIONS
% % 1) What components of the original signal remain? What components are
% % attenuated?
% % 2) If you wanted to attenuate the 'Treble' component which filter would
% % you use?
% % 3) If you wanted to attenuate the 'Bass' component which filter would you
% % use?
% % 4) If you were to physically build this filter (i.e. resistors,
% % capacitiors), what values of resistors and capacitors would you use?
% 
% % Frequency domain analysis
% figure;
% plot(f, 10*log10(pxx))
% hold all
% [LP_pxx, f] = periodogram(LoPassSig,[],[],fs);
% plot(f,10*log10(LP_pxx));
% axis([0 1000 -120 -30]);
% 
% %A more extreme Low-Pass Filter. Cut-off = 50 Hz
% LP2=designfilt('lowpassiir', 'FilterOrder', 1, 'HalfPowerFrequency', 50, 'SampleRate', 22050);
% LoPassSig2 = filter(LP2,micSig(:,1));
% audiowrite('Data\LoPass2.wav',LoPassSig2,fs);
% [LP_pxx2, f] = periodogram(LoPassSig2,[],[],fs);
% plot(f,10*log10(LP_pxx2));
% 
% 
% figure;
% plot(f, 10*log10(pxx))
% hold all
% [HP_pxx, f] = periodogram(HiPassSig,[],[],fs);
% plot(f,10*log10(HP_pxx));
% axis([0 1000 -120 -30]);
% 
% 
% 
% %% Filtering (part 2)
% % Running average is a type of Low pass filter. The length of the 'run' is
% % inversly proptional to the cut-off frequency
% %
% % Signal = X1, X2, X3, X4, X5 ...
% % Window size = n points, filtered value F1 = 1/n*(X1+X2+X3...Xn).
% % average of the first n points
% % known as a finite impulse response filter (finite because there is a
% % defined window size)
% 
% %building a 50 Hz running average window
% %window size = sampling rate * inverse of cutoff freq.
% windowSize = round(fs*1/50);
% LoPassSig3 = filter(1/windowSize*ones(1,windowSize),1,micSig(:,1));
% figure;
% plot(T,micSig(:,1))
% hold all
% plot(T,LoPassSig2)
% plot(T,LoPassSig3)
% axis([5, 5.1, -.2 .2]);
% xlabel('Time (s)')
% ylabel('Amplitude (V)')
% legend([{'Original'};{'50 Hz Low-Pass Butterworth (RC)'};{'50 Hz Low-Pass FIR rect window'}]);
% 
% %Questions
% % 1) Edit code for a window with of 500 Hz. How many samples are averaged
% % for each point of this filter?
% % 2) How would you make a FIR high-pass filter?

%% Importing Image data - It's really a 2D Array

neuronPhase = imread('Data\8035-03d_0024-site2-phase.png');
neuronPhase(1:10,1:20,1)
neuronPhase(1:10,1:20,2)
neuronPhaseG=rgb2gray(neuronPhase);
figure; imshow(neuronPhaseG);
title('original image');

% Histogram of image (x axis is intensity, y axis is counts [frequency])
figure;
h = histogram(neuronPhaseG);
title('histogram of intensities, original image');
xlabel('intensity')
ylabel('Count (# pixels');

% Questions
% 1) In a 16 bit image, what is the minimum pixel value?, maximum pixel
% value? number of distinct pixel values?
% 2) How would you increase the brightness of the image?
% 3) How would you increase the contrast of the image?

%multiply each value by 2
figure; 
imshow(neuronPhaseG*2)
title('Image x2 intensity')
figure; histogram(neuronPhaseG*2);
title('Image x2 intensity')
xlabel('Intensity');
ylabel('Count (# pixels');

%add 20k to each value
figure;
imshow(neuronPhaseG+20000)
title('Image + 20K intensity')
figure; histogram(neuronPhaseG+20000);
title('histogram (image +20K)')
xlabel('Intensity');
ylabel('Count (# pixels');

%counversion to Double
neuronPhaseGD = double(neuronPhaseG);

%Max pixel value
maxV = max(max(neuronPhaseGD))
%Min pixel value
minV = min(min(neuronPhaseGD))



% 99%tile
maxV = quantile(neuronPhaseGD,.99,'all');
% 1%tile
minV = quantile(neuronPhaseGD,.01,'all');


imRange = maxV-minV;

figure;
neuronPhaseNormalized = uint16(2.^16/imRange*(neuronPhaseGD-minV));

imshow(neuronPhaseNormalized);
title('Hi-Contrast Image (normalized)')
figure; histogram(neuronPhaseNormalized);
title('histogram (normalized)')
xlabel('Intensity');
ylabel('Count (# pixels');

%% Spatial filtering

% create a 'small' 10 x 10 pixel image with values from 0 - 10
DemoMat =randi(10,10)

% create spatial filter
% filter subtracts the left value and adds the right value.
h = [-1 0 1];

% Filtered matrix
FiltDemo = imfilter(DemoMat,h)

% Questions
% 1) Can you come up with an 'averaging filter', what effect would this
% have on the image?
% 2) The filter above is a difference filter. The reported value is the
% difference between the adjacent right and left pixels. What features of
% an image would this accentuate? How might this be useful?


%fsize filter size
fsize = 250;
hGaussian = fspecial('gaussian', fsize, 20);

neuronPhaseNorm_GaussFilt = imfilter(neuronPhaseNormalized,hGaussian);


%Show gaussian filtered image
figure; 
imshow(neuronPhaseNorm_GaussFilt);
title('Gaussian Filtered Image')

%Questions 
% 1) what does the following image show? 
figure;
imshow(uint16(neuronPhaseNormalized-neuronPhaseNorm_GaussFilt))

% 2) how could a gaussian filter be used to correct a gradient in an image? 


% median filtering to 'denoise' an image

figure;
imshow(medfilt2(neuronPhaseNormalized, [15 15]));
title('Median Filter 15 x 15 neighborhood')

neuronPhaseMedFilt = medfilt2(neuronPhaseNormalized, [15 15]);


%% Edge Detection in an image
% use a high-pass spatial filter (Canny) to find the max local deriatives
% of an image:


neuronEdge = edge(neuronPhaseMedFilt,'Canny',[.1 .3]);
figure;
imshow(neuronEdge);
title('edge detection');























