classdef SweepTDT < handle
    
   % Class to contain a signal to play from the GUI guideSetupSystem3. 
   % This class is only used as input to some properties of the class
   % AlbanTDT
   %
   % Methods are built so that when this object is played through System3,
   % it will be played at desiredLevel (property of AlbanTDT object).
   %
   % Most properties are private (for writing only) to ensure that the
   % whole thing behaves as expected. Run SweepTDT.explainAttenuation() for
   % explanations.
   
   properties
       % Filename of the wav file if sound is loaded with audioread
       origin@char
       % Text being said or description of the waveform
       text@char
   end
   
   properties (SetAccess = private)
       % Object of class AlbanTDT used to create a signal that will behave
       % according to how the TDT system is setup (this is controlled by 
       % the AlbanTDT object).
       tdt@AlbanTDT
       % Original frequency of the signal (differs from sound_sampling_freq if resampled)
       originalfreq@double
       % The digital signal that we want to play at the desired level, if
       % the system was flat (but it's not, so we distort the signal into
       % signal_left and signal_right so that we counter-distort the signal)
       signal@double
       % Signal distorted to compensate the system's unflatness for left ear
       signal_left@double
       % Signal distorted to compensate the system's unflatness for right ear
       signal_right@double
       % Attenuation for left PA5 (together with left RX8, should be 30)
       attenPA5_left@double
       % Attenuation for right PA5 (together with right RX8, should be 30)
       attenPA5_right@double
       % Attenuation on left signal within RX8 (together with left PA5, should be 30)
       attenRX8_left@double
       % Attenuation on right signal within RX8 (together with right PA5, should be 30)
       attenRX8_right@double
       % Maximal tone output of the signal (normally 120dB)
       maxToneOutput
       % Headamp attenuation value (normally 20dB)
       headampAtten
       % This rms value is calculated after removing silences from the
       % raw signal (they might be long for a speech waveform)
       signalRmsNoSil@double
       % Digitally, the sweep should have an rms of 10/sqrt(2), which is
       % equivalent to a sound played at 120dB (assuming pure tone and flat
       % system). This will be reduced to desiredLevel by the RX8, PA5 and
       % headamp devices.
       desiredRms = 10/sqrt(2)
       % This is the maximum (in absolute value) the signals played can
       % have without clipping by the RX8. To ensure this, we check the
       % maximal value of the signal_left and _right, and divide it by 10 as
       % many times necessary. We reduce
       maxValSignal
       % Duration of the stimulus (in ms, set when loading the signal)
       stimulus_duration@double
       % Delay in the RX8 system (in ms, set when initialising the SweepTDT object)
       stimulus_delay@double
       % Duration of the recording by RZ2 (in ms, set when loading the signal)
       record_duration@double
       % Sound sampling at which the sounds will be played (in Hz, set when initialising)
       sound_sampling_fq
       % Sound level (in dBSPL) at which the sounds should be played by
       % TDT. Since this is used to decide the calibrations,
       % setSignals_leftAndRight is used whenever this property is changed
       % with the property setDesiredLevel
       desiredLevel@double     % 70dB normally
   end
   
   methods 
       
       %%%%%%%%%%%%%%%%%   INITIALISATION   %%%%%%%%%%%%%%%%%
       
       function sweep = SweepTDT(tdt)
           sweep.tdt                = tdt;
           sweep.maxToneOutput      = tdt.maxToneOutput;
           sweep.maxValSignal       = tdt.maxValSignal;
           sweep.desiredLevel       = tdt.desiredLevel;
           sweep.headampAtten       = tdt.headampAtten;
           sweep.stimulus_delay     = tdt.stimulus_delay;
           sweep.sound_sampling_fq  = tdt.sound_sampling_fq;
           if ~isequal(sweep.headampAtten, 20)
               error('IHR convention: headampAtten should be 20!');
           end
       end
       
       %%%%%%%%%%%%%%%% DO IT ALL %%%%%%%%%%%%%%%%%%
       
       function setItAll(sweep, wavOrSignal)
           % Sets the signal given (as path, cell or array) and sets the
           % left and right parameters. If you want to use a specific sound
           % level, use setDesiredLevel BEFORE calling setItAll, for
           % efficiency.
           sweep.setSignal(wavOrSignal)
           sweep.setSignals_leftAndRight();
       end
       
       %%%%%%%%%%%%%%    SETTING THE INITIAL SIGNAL   %%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%        AT THE RIGHT LEVEL      %%%%%%%%%%%%%%%%
       
       function setSignal(sweep, wavOrSignal)
           % Set the signal after obtaining its rms without long silences
           
           % Unwrap the input
           [sig, Fs] = unwrapSignal(sweep,wavOrSignal);
           
           % Resample if necessary
           sweep.originalfreq = Fs;
           if Fs ~= sweep.sound_sampling_fq
               sig = resample(sig, sweep.sound_sampling_fq, Fs);
               Fs = sweep.sound_sampling_fq;
           end
           
           % Some options
           option.doPlots = 0;
           
           % Calculate the RMS without the long silences (if any)
           sweep.signalRmsNoSil = sweep.getRmsWithoutSilences(sig,Fs,option);
           
           % Must be column otherwise RX8 won't understand (tdt.continuousAcquire)!!
           if size(sig,2)==1
               sig = sig';
            end
           
           % Save in sweep.signal after renormalising to desiredRms. Should
           % be slightly smaller than 7, due to the silences
           sweep.signal = sig * sweep.desiredRms / sweep.signalRmsNoSil;
           
           % Set additional info:stimulus_duration, record_duration (in ms)
           sweep.stimulus_duration = 1000 * length(sig) / sweep.sound_sampling_fq;
           sweep.record_duration = sweep.stimulus_duration + sweep.stimulus_delay;
       end
       
       function [sig, Fs] = unwrapSignal(sweep,wavOrSignal)
           % Wrapper to obtain the digital signal from path, cell
           
           switch class(wavOrSignal)
              
               case 'char'   % Read file
                   if ~exist(wavOrSignal, 'file')
                       error('The string is assumed to be a path to a wav file');
                   end
                   [sig, Fs] = audioread(wavOrSignal);
                   sweep.origin = wavOrSignal;
                   
               case 'cell'   % Take inputs
                   switch length(wavOrSignal)
                       case 1,  sig = wavOrSignal{1}; Fs = sweep.tdt.sound_sampling_fq; 
                       case 2, [sig, Fs] = wavOrSignal{:};
                       otherwise, error('If wavOrSignal is a cell, should have 1 or 2 elements.');
                   end
                  
               case 'double'  % Signal given as array
                    sig = wavOrSignal;
                    Fs = sweep.tdt.sound_sampling_fq; % 48818Hz noramally

               otherwise
                   error('wavOrSignal''s class not recognised.');
           end
       end
       
       function  signalRms = getRmsWithoutSilences(sweep,signal,Fs,option)
           % Calculate the RMS of the signal without its long silences
           
           % Get indices of beginning and end of silences
           [beginSil,endSil] = sweep.getSilenceBeginEndIndices(signal,Fs,option);
           
           % Get rid of long silences
           [signalNoSil] = sweep.cutOffLongSilences(signal,beginSil,endSil,option);
           
           % Calculate rms
           signalRms = norm(signalNoSil,2)/sqrt(length(signalNoSil));
       end
       
       function [beginSil,endSil] = getSilenceBeginEndIndices(~,signal,Fs,option)
           % Get a list of indices, when silences begin and end
           
           % Threshold: 1/15th of the max
           m = max(signal);
           th = signal>(m/15);

           % A 'long' silence lasts for at least 0.05s 
           durationLongSilence = 0.05; % in s
           
           % Number of bins to consider a long silence
           timebinLongSilence = 2*floor(floor(durationLongSilence*Fs)/2); 
           
           % Convolve and shift to recenter
           con = circshift(conv(ones(timebinLongSilence,1),double(th)), [-timebinLongSilence/2,0]);
           
           % 'Silent' bins
           sil = (m*con/max(con)<0.005);
           
           if size(sil,2)==1
               sil = sil';
           end
                      
           % Indices when silences start and end
           chang    = find([1, diff(sil)]);
           beginSil = chang(sil(chang)==1);
           endSil   = chang(sil(chang)==0);
           if endSil(end)<beginSil(end)
               % to simplify the next operation
               endSil(end+1) = length(signal) - 10; 
           end
           if endSil(1) == 1
               endSil = endSil(2:end);
           end
           
           if length(beginSil)~=(length(endSil))
               error('Not the same length, problem');
           end
       end
       
       
       function [signalNoSil] = cutOffLongSilences(~,signal,beginSil,endSil,option)
           % Makes the signal silence-free
           
           % Replace silences
           for kk=1:(length(beginSil))
               cBeg = beginSil(kk);
               endSil(kk) = min(length(signal),endSil(kk)); % Not thoroughly tested, seems ok
               if cBeg>length(signal)
                   continue;
               end
               signal((cBeg+1):(cBeg+10)) = linspace(signal(cBeg), signal(endSil(kk)),10);
               signal((cBeg+11):(endSil(kk)-1)) = NaN;
           end
           signalNoSil=signal(~isnan(signal));
           
           
           if option.doPlots
               for jj=beginSil % where long silence starts
                   plot(jj,0.25, '.b');
               end
               for jj=endSil % where long silence starts
                   plot(jj,0.225, '.r');
               end
               plot(0.35+signalNoSil/(10*max(signalNoSil)))
           end
           
       end
       
       %%%%%%%%%%%%%%%%%  SET LEFT AND RIGHT SIGNALS AND  %%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%          ATTENUATIONS            %%%%%%%%%%%%
       
       function setSignals_leftAndRight(sweep)
           % Sets left and right signals
           if isempty(sweep.signal)
               error('sweep.signal should already be set.');
           end
           len = length(sweep.signal);
           NFFT = 2^nextpow2(len); % Next power of 2 from length of y
           Y = fft(sweep.signal,NFFT);
           sweep.setSignal_left_right(Y,len,'left');
           sweep.setSignal_left_right(Y,len,'right');
       end
       
       function setSignal_left_right(sweep,Y,len,side)
           % Set signal_left (or _right) using the calibration. Also sets
           %  sweep.(['attenRX8_' side])  and   sweep.(['attenPA5_' side])
           switch side
               case {'left', 'right'}
               otherwise, error('%s should be left or right');
           end
           if isempty(sweep.tdt.(['calib_polyfit_' side '_log10']))
           %if isempty(sweep.tdt.(['calib_polyfit_' side]))
               if isempty(sweep.tdt.(['calib_' side]))
                   warning('You need to load the %s calibrations before calculating the sweep correction.\n', side);
                   return;
               end
               %warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
               %cFit = polyfit(sweep.tdt.calib_freqs, sweep.tdt.(['calib_' side]),3);
               %sweep.tdt.(['calib_polyfit_' side])  = cFit;
               %warning('on','MATLAB:polyfit:RepeatedPointsOrRescale')
               %else
               %    cFit       = sweep.tdt.(['calib_polyfit_' side]);
               
                sweep.tdt.(['calib_polyfit_' side '_log10']) = ...
                    sweep.tdt.fitCalibration(sweep.tdt.calib_freqs, sweep.tdt.calib_right);
           end
           
           cFit_log10 = sweep.tdt.(['calib_polyfit_' side '_log10']); 
           NFFT = length(Y);   
           Fs = sweep.sound_sampling_fq; % sampling freq of signal
           f = Fs/2*linspace(0,1,NFFT/2+1);
                      
        
           
           if 1==0 % TESTS only!!!
               f_log = logspace(1,log10(Fs/2),NFFT/2+1);
               % freqToUse = sweep.tdt.calib_freqs
               % calibToUse = sweep.tdt.(['calib_' side]) ;
               f_minToUse = 100;   %Hz, decided by Trevor on Juy 12th 2016
               f_maxToUse = 14000; %Hz CHANGED AFTERWARDS
               indToUse   = ...
                   sweep.tdt.calib_freqs < Fs/2        & ...
                   sweep.tdt.calib_freqs > f_minToUse  & ...
                   sweep.tdt.calib_freqs < f_maxToUse;
               freqToUse  = sweep.tdt.calib_freqs        ( indToUse );
               calibToUse = sweep.tdt.(['calib_' side])  ( indToUse );
               % Log-spread them
               %freqToUse  = logspace(log10(freqToUse(1)), log10(freqToUse(end)), length(freqToUse))';
               indToUse = floor([1:9 logspace(1.05, log10(length(freqToUse)),60)])';
               freqToUse  = freqToUse(indToUse);
               calibToUse = calibToUse(indToUse);
               
               % Make it flat under f_min and above f_max (Trevor's choice)
               freq_before_fMin = 1:10:f_minToUse-1;
               freq_after_fMax = logspace(log10(f_maxToUse+1),log10(Fs/2), 30);
               freqToUse  = [freq_before_fMin   freqToUse'   freq_after_fMax];
               
               calib_before_fMin = linspace(120, calibToUse(1), length(freq_before_fMin)); %     calibToUse(1)*ones(size(freq_before_fMin)) ;
               calib_after_fMax  = calibToUse(end)*ones(size(freq_after_fMax));
               calibToUse = [calib_before_fMin    calibToUse'    calib_after_fMax];
               
               % DVPT: to delete
               for lev = 9 %3:9;
                   % cFit_log10 = polyfit(log10(sweep.tdt.calib_freqs(2:end)), sweep.tdt.calib_left(2:end),lev);
                   cFit_log10 = polyfit(log10(freqToUse(1:end)), calibToUse(1:end),lev);
                   
                   % Approximate calibration at the required frequencies, using the polynomial fit
                   % calculated by tdt.load_calibration_single (or above) when loading a calibration
                   %curvef = polyval(cFit,f_log);
                   curvef_log10 = polyval(cFit_log10,log10(f_log));
                   
                   figure;
                   subplot(1,2,1)
                   %  plot(f_log, curvef, freqToUse, calibToUse );
                   subplot(1,2,2)
                   plot(log10(f_log), curvef_log10, log10(freqToUse),calibToUse);
                   ylim([70 140])
                   title(lev)
               end
               
               % DVPT: with splines
               cFit_spline_log10 = spline(log10(freqToUse), calibToUse, log10(f));
               figure; semilogx(freqToUse, calibToUse, (f), cFit_spline_log10);
               title('spline')
               figure; plot(freqToUse, calibToUse, (f), cFit_spline_log10);
               
               
               % DVPT: with interp1
               cFit_interp1_log10 = interp1(log10(freqToUse),calibToUse, log10(f),'spline');
               figure; semilogx(freqToUse,calibToUse, (f), cFit_interp1_log10);
               title('interp1')
               
     
               
               %DVPT: with smoothing splines : DON'T HAVE THE TOOLBOX
               %p=0.1:0.1:1;
               %xxi=log10(f);
               %yy = zeros(length(p),length(xxi));
               %for j=1:length(p)
               %  yy(j,:) = csaps(log10(sweep.tdt.calib_freqs), sweep.tdt.(['calib_' side]),p(j),xxi);
               %end
               
               %curvef = polyval(cFit,f);
               % curvef = polyval(cFit_log10,f);
               %curvef_log10 = polyval(cFit_log10,log10(f));
               
               %cFit_spline_log10
               %curvef_log10 = polyval(cFit_log10,log10(f_log))
               
               %cFit_spline_log10 = spline(log10(freqToUse), calibToUse, log10(f));
               %   figure; semilogx(freqToUse, calibToUse, (f), cFit_spline_log10);
               
               
           end
           
            switch class(cFit_log10)
                case 'double', curvef = polyval(cFit_log10,log10(f)); % We used polyfit
                case 'function_handle', curvef = cFit_log10(f);       % We used splines     
                otherwise, error('Class %s not recognised',class(cFit_log10));
            end
        
          %curvef = curvef_log10(2:end);  % first is Inf given how we calculated
          curvef(1) = sweep.maxToneOutput;
          curvef(2) = sweep.maxToneOutput; % second cn also be terribly bad
          
          if any(curvef < 90 | curvef > 150)
              warning('Some values of curvef are problematic');
          end
          
           % figure; plot(curvef, (freqs), calib); xlabel('Frequency (Hz)'); ylabel('Sound level (dB SPL)');
           % legend('Recorded','Fitted Curve') ; ylim([70 125]); set(gcf, 'PaperPosition',[0 0 6 4]); set(gcf, 'PaperSize',[6 4]);
           %  a = 'path\to\folder'; saveas(gcf, [a 'pdf']); saveas(gcf, [a 'fig']);
           
           % In other files, applying the calibration should look like this:
           %[~, i] = arrayfun(@(x) min(abs( x - freqs )), f);
           %scalars = 10 .^  ((calib(i) - maxToneOutput) / 20);
           
           % To correct the calibration, we want the dB to add to signal
           vecCalib = 10 .^ ((sweep.maxToneOutput - curvef) / 20);
           if size(curvef,2)==1
               vecCalib = vecCalib';
           end
           
           % Apply the weights and inverts the Fourier representation
           % Y starts at 2 because we removed the first value of curve_f
           newSignal = ifft([Y(1:NFFT/2).*vecCalib(1:end-1), Y(NFFT/2+1:NFFT).*flip(vecCalib(2:end))],'symmetric');
           
           % Remove the padding
           sig = newSignal(1:len); 

           % The attenuation of the system RX8+PA5 (30dB, normally)
           attenuationRX8PA5 = sweep.maxToneOutput - sweep.headampAtten - sweep.desiredLevel;

           % Calculate appropriate attenuations in RX8 and PA5 for this side
           maxVal = max(abs(sig));           
           if maxVal > sweep.maxValSignal
               % Calculate the minimal attenuation to put in RX8 such that
               % the values it outputs are in [-10,10]. This maximises the
               % dynamic range used by the system.
               attenRX8 = ceil(20*log10(maxVal/sweep.maxValSignal)*10)/10; 
               % (*10)/10 because PA5 precise up to 0.1
               % ceil because RX8 should rather attenuate a bit more than
               % not enough (clipping)
               attenPA5 = attenuationRX8PA5 - attenRX8;
           else
               % RX8 doesn't need to attenuate
               attenRX8 = 0;
               attenPA5 = attenuationRX8PA5;
           end
      
           % Final check
           if attenPA5<0 % ERROR 
               if attenPA5>-20 % Attempt to still play it (only for development)
                   warning (['attenPA5 is negative (%.2f), we can''t play the sound at the desired level. Played at %.2fdB instead of %.2fdB.' ...
                       ' For now we make it 0 (we assume with a good calibration this won''t happen)'], ...
                       attenPA5,...
                       sweep.maxToneOutput - sweep.headampAtten - attenPA5 - attenRX8,...
                       sweep.desiredLevel);
                   attenPA5 = 0;  
               else % Too loud, might damage the hardware
                   error(['The PA5 attenuation is set to a very negative value. ' ...
                       'This is dangerous!!! The calibration is probably wrong, ' ...
                       'hence we''d need to amplify the sounds.\n']);
               end
           end
           
           % Save the side's signal, RX8 atten, PA5 attenuation
           if size(sig,2)==1 % Must be column otherwise RX8 can't understand (tdt.continuousAcquire)!!
               sig = sig';
           end
           sweep.(['signal_' side]) = sig;
           sweep.(['attenRX8_' side]) = attenRX8;
           sweep.(['attenPA5_' side]) = attenPA5;
           
       end
       
       %%%%%%%%%%%%% TO SET LEVEL %%%%%%%%%%%%%%%%%
       
       function setDesiredLevel(sweep, desiredLevel)
           % This function sets the desiredLevel and calls setSignals_leftAndRight
           % if a signal was already in sweep.signal.
           sweep.desiredLevel = desiredLevel;
           if ~isempty(sweep.signal)
               sweep.setSignals_leftAndRight();
           end
           
       end
   end
  
   
   %%%%%%%%%%%%%%%%%    Explain the system and its constraints   %%%%%%%%%%%%%%%%
   
   methods (Static)
       function explainAttenuation()
           % Display text to explain the system and its constraints
           s={};
           s{end+1}='                                                                                                                                              ';                            
           s{end+1}=' This aims at explaining the construction of the SweepTDT method under System3 constraints.                                                   ';
           s{end+1}=' We work under the following constraints:                                                                                                     ';
           s{end+1}='    - The desired sound level is given (desiredLvl=70dB normally),                                                                            ';
           s{end+1}='    - The headamp attenuator is meant to make sure sounds are not too loud (IHR''s convention: headampAtten=20dB),                            ';
           s{end+1}='    - The RX8 can only output values between -10 and 10 (voltage, analog signal),                                                             ';
           s{end+1}='    - Except for the relationship rms=10/sqrt(2) <-> level=120dB from Matlab, we cannot translate a given sound level into an rms value.      ';
           s{end+1}='                                                                                                                                              ';                            
           s{end+1}='  So the system RX8+PA5 is meant to use a maximal dynamic range (analog values as high as possible)                                           ';
           s{end+1}='  while delivering a signal at an appropriate level (Atten=120-20-70=30 normally)).                                                           ';
           s{end+1}='                                                                                                                                              ';                            
           s{end+1}='                                         SCHEMATIC PRESENTATION OF ATTENATION PROTOCOLE                                                       ';
           s{end+1}=' ____________________________________________________________________________________________________________________________________________ ';
           s{end+1}='|                                                                                                                                            |';
           s{end+1}='|                          lvl=120dB                       lvl=lvl-headampAtten-desiredLvl   lvl=lvl-headampAtten      lvl==desiredLvl       |';
           s{end+1}='|                       (rms=10/sqrt(2))                                                                                                     |';
           s{end+1}='|                                       _______________________________                                                                      |';
           s{end+1}='|                                      |                               |                                                                     |';
           s{end+1}='|  RawSignal --> Matlab  ----------------> RX8 ----------------> PA5 ----------------> headamp ----------------> speakers ----------> animal |';
           s{end+1}='|   (.wav)                 (digital)             (analog)               (analog)                   (analog)                (sound)           |';
           s{end+1}='|____________________________________________________________________________________________________________________________________________|';
           s{end+1}='                                                                                                                                              ';
           s{end+1}=' Matlab also uses the system''s calibration to distort the signal such that the response is flat                                              ';
           s{end+1}=' (when playing pure tones are 120dB, this calibration) ensures that it''s indeed as loud as expected).                                        ';
           s{end+1}=' Look at applyCalibration.m for details on how this is done.                                                                                  ';                            
           fprintf(1, '%s\n',s{:});
       end
   end
end