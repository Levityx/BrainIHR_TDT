function splcal_trevors(filename)
% Calibrates main phys lab sound system to SPL
% 
% Noise bursts played to both ears, but only single microphone used so you
% need to swap between ears.
%
% Filename is where calibration will be stored - remember to give different
% names for left and right
%
% Original file name: splcal

mainphys_param.mic_sensitivity = -38.9; %dB re 1V per Pa (94 dB)
%mainphys_param.mic_sensitivity = -60.9; %1/8" dB re 1V per Pa (94 dB)
mainphys_param.amp_gain = 50; %dB
mainphys_param.adjust = 0; %dB
mainphys_param.probe_file = 'Path/to/short_B&K_probe.spl';
mainphys_param.noise_attn = 0; %dB

corr_param = mainphys_param;

if nargin < 1
    filename = 'Path/to/'
end

system('Path/to/set_headamp_00.exe');

% Set up the zbus control. 
tdt.zBus=actxcontrol('ZBUS.x',[1 1 1 1]);
tdt.zBus.ConnectZBUS('GB');

% Set up  processors.
tdt.RX = setUpRX8_cal('Path/to/noise_calibration.rcx');
for i=1:2  tdt.PA(i)= connectPA5(i); end % changed into a cell.

no_repeats = 10;
noise_amp = 1; opch = 3; 
record_duration = 1;
stimulus_duration = 0.8;

% play a short noise burst to prime buffers and correct any DAC/ADC drift
% over time or since rebooting PC (no idea if this is necessary!)
play_noisecal(noise_amp,opch,tdt, 0.2, 0.1, corr_param.noise_attn);

for rep = 1:no_repeats
    fprintf('Repeat %3.0f ',rep);

    [noise_data{rep} samplerate] ...
    = play_noisecal(noise_amp,opch,tdt, ... 
              record_duration, stimulus_duration, corr_param.noise_attn);          

    [corr_response(:,rep) freqs correction]...
    = correct_frequency_response(noise_data{rep}, ...
                                    record_duration, mainphys_param );

    figure(1);
    avg_corr_rep=mean(corr_response,2);
    semilogx (freqs, avg_corr_rep);
    axis([50 50000 60 140]);
    drawnow
end

% Alban: to improve GUI
[pname,fname, ext] = fileparts(filename);
if strcmp(ext, '.mat')
    filename = fullfile(pname, fname);
end

if isempty(strfind(filename, 'freq_resp'))
    fname_resp = [filename '_freq_resp.mat'];
else
    fname_resp = [filename '.mat'];
end
fprintf('\nSaving frequency response %s\n', fname_resp)
save(fname_resp, 'avg_corr_rep', 'freqs');
filename = strrep(filename,'freq_resp','data');
if isempty(strfind(filename, '_data'))
    fname_data = [filename '_data.mat' ];
else
    fname_data = [filename '.mat'];
end

fprintf('Saving raw data %s\n', fname_data)
save(fname_data, 'noise_data', 'samplerate', 'correction', 'mainphys_param');




