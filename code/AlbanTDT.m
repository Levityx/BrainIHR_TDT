classdef AlbanTDT < handle
    % Class to run an experiment with TDT hardware.
    % Entirely managed by guideSetupSystem3.m
    %
    % Author: Alban Levy, April 2016
    
    properties
        circuitPathReco@char
        circuitPathPlay@char
        listWavFiles % @cell||struct
        
        % Folder where we will record the current data, plots...
        local@char
        
        no_repeats = 1
        noise_amp = 1
        
        % internal gain in circuit
        recordinggain = 1
        
        % opch=3 is the expected value: previously used to decide between 1. left ear 2. right ear 3. both
        opch = 3
        
        % Rate at which the sound will be played by RX8, in Hz
        sound_sampling_fq %= 48828
        
        % Number of channels to record (128 by default)
        nb_channels@double
        
        % Array of electrodes to record from (delete the others)
        arrayChanRec@double
        
        % Size of the RX8 to load the wave files (implies a limitation on the sound's length)
        sizeBuf =  1280000; % 10000128: too big
    
        % Connect to devices and zBus control, load circuits:
        RZ
        RX
        zBus
        PA_left
        PA_right
        
        % stimulus delay before recording the neural responses (in ms)
        stimulus_delay = 50; %ms
        
        % Default sampling frequency we want
        RZ_SamplingFreq % = 24414 Hz
        
        % Filter coefficients used on the recording ([b,a] = butter...)
        channels_filt_b@double
        
        % Filter coefficients used on the recording ([b,a] = butter...)
        channels_filt_a@double
        
        % Figure to plot the data being acquired
        figure
        
        % Running index (updated while playing different .wav fils)
        runningWavIndex@double
        
        % Boolean property to stop the loop from the GUI
        stopTheLoop@double
        pauseTheLoop@double
        
        % Sound level parameters
        maxValSignal = 10 % signal given to System 3 in [-10 10]
        maxToneOutput = 120 %dB
        desiredLevel = 70 %dB
        headampAttenExpected = 20 %dB
        headampAtten
     
        % Properties to be defined at each experiment
       
        % Directory in which are the folders, scripts, ...
        masterDirectory 
        % ID used for puretonograms and experiment
        ID@char
        % path to the folder 'data'; the IP of RS4 changes everytime
        RS4IP@char
        % Depth of the electrode
        depth@double
        % Minimum time at which we calculate the activity to plot (in ms)
        minTimeActivity@double
        % Max time at which we calculate the activity to plot (in ms)
        maxTimeActivity@double
        
        
        % Calibration values (left) to use when we give sounds to the system 3
        calib_left@double
        % Calibration values (right) to use when we give sounds to the system 3
        calib_right@double
        % Calibration values (frequencies) to use when we give sounds to the system 3
        calib_freqs@double
        % Keep path of the left calibration file
        calib_left_file@char
        % Keep path of the right calibration file
        calib_right_file@char
        
        % Polynomial fits of the calibration: polyfit(freqs,calib_L,3);
        % calculated by load_calibration_single when loading a calibration
        calib_polyfit_left  
        % Same polynomial fit, but applied on the logarithme scale (log10
        % was applied on the X-axis: polyfit(log10(freqs),calib_L,3))
        calib_polyfit_left_log10
        % Polynomial fits of the calibration: polyfit(freqs,calib_R,3);
        % calculated by load_calibration_single when loading a calibration
        calib_polyfit_right
        % Same as left case
        calib_polyfit_right_log10
        
        % Sweep object that will contain the processed signal to play
        sweep@SweepTDT
        
    end
    
    methods
        function tdt = AlbanTDT()
            tdt.stopTheLoop = 0;
            tdt.pauseTheLoop = 0;
            tdt.minTimeActivity = 0;
            tdt.maxTimeActivity = Inf;
            tdt.desiredLevel = 70; % dB
            tdt.nb_channels = 128;
            tdt.arrayChanRec = 1:tdt.nb_channels;
            
            tdt.RZ_SamplingFreq = 24414;   % record activity, in Hz
            tdt.sound_sampling_fq = 48828; % play sounds, inHz
            HPFreq = 300; 
            [tdt.channels_filt_b, tdt.channels_filt_a] = butter(4,HPFreq/(tdt.RZ_SamplingFreq/2), 'high'); % third order, cut-off 500Hz, high-pass
        end
        
        function tdt = connect(tdt, varargin)
            % Connect to devices and zBus control, load circuits:
            % - Recording channels: RZ2
            % - Playing sound: RX8
            % - Send triggers A and B: zBus
            % - Attenuate sound for left or roght ear: PA_left, PA_right
            for kk=varargin
                arg = kk{1};
                switch arg
                    case 'RX8', tdt.RX = Circuit_Loader_device('GB', 1, tdt.circuitPathPlay, 'RX8');
                    case 'RZ2', tdt.RZ = Circuit_Loader_device('GB', 1, tdt.circuitPathReco, 'RZ2');
                        tdt.RZ.SetTagVal('sizeBuf',tdt.sizeBuf); % Manually put to a bit more, because otherwise it would get approximated
                        if tdt.sizeBuf~=tdt.RZ.GetTagVal('sizeBuf')
                            error('The buffer size is not as planned: put it to %d manually and reconnect the device. Also try rebooting the hardware', tdt.sizeBuf);
                        end
                    case 'ZBUS',  tdt.zBus = actxcontrol('ZBUS.x',[1 1 1 1]);
                        tdt.zBus.ConnectZBUS('GB');
                    case 'PA_left',   tdt.PA_left  = connectPA5(1);
                    case 'PA_right',  tdt.PA_right = connectPA5(2);
                    otherwise, error('The device %s isn''t known in method connect', arg);
                end
            end
        end
        
        function checkAndSetHeadAmp(tdt)
            % Checks whether tdt.headampAtten was set at 
            % tdt.headampAttenExpected (20dB). If not, call
            %       >> system('Path/to/set_headamp.exe');
            % The user is expected to manually enter 20 in both, press enter,
            % close the small window.
            if isempty(tdt.headampAtten)||tdt.headampAtten~=tdt.headampAttenExpected
                fprintf('Manually enter %d and %d in the window, ''Change'' it and close the window\n',...
                    tdt.headampAttenExpected,tdt.headampAttenExpected);
                h = errordlg(['After changing the headamp to ' num2str(tdt.headampAttenExpected) ...
                    ', press enter'], 'Change headamp');
                error('Provide path to set_headamp');
                system('Path/to/set_headamp.exe');
                close(h);
                tdt.headampAtten = tdt.headampAttenExpected;
                
            end
        end
        
        function [rec, e] = continuousAcquire(tdt, sweep)
            % Main method of the class AlbanTDT. Takes as input an object
            % of class SweepTDT, will load sweep.signal_left and
            % sweep.signal_right into the RX8, set the attenuation at the
            % appropriate levels, start playing to sound to both ear,
            % record on the RZ2 - and on the RS4 if connected to the
            % network.
            %
            % For an explanation on the  attenuation protocole, run
            %     >> SweepTDT.explainAttenuation 
            
            if ~isa(sweep,'SweepTDT')
                error('This function assumes an input of class SweepTDT');
            else % Check a few things
                if any(arrayfun(@isempty,[sweep.attenRX8_right sweep.attenPA5_left sweep.attenRX8_left sweep.attenRX8_right]) )
                    fprintf('ContinuousAcquire: Some fields of sweep are empty. There might be a problem.\n');
                    disp(sweep);
                end
            end
            
            % Set headamp to -20dB
            tdt.checkAndSetHeadAmp();
            
            % Main function to reset the circuits (zBusB trigger) and start
            % running them (zBusA trigger)
            
            % In batch mode, the main worker may be connected to the hardware. We
            % need this worker to connect instead
            if  tdt.RZ.GetTagVal('sizeBuf')==0
                fprintf('continuousAcquire: Connecting to the hardware because it''s a new batch...\n');
                tdt.connect('RX8', 'RZ2', 'ZBUS');
            end
            
            % if figure missing, we start a new one
            if isempty(tdt.figure)||~isvalid(tdt.figure),  
                tdt.figure = guideChannelsPlot();
            end
            
            % Frintfing options
            guidata_tdt = guidata(tdt.figure);
            fprintopt = get(guidata_tdt.checkboxFprintfInfo,'Value');
            
            % We reinitialise tdt.stopTheLoop
            if tdt.stopTheLoop 
                tdt.stopTheLoop=0;
            end
            
            %tdt.RZ.SetTagVal('sizeBuf',tdt.sizeBuf); % 10000128: too big, hence 1280000
            npts              = tdt.RZ.GetTagVal('sizeBuf'); % number of 32bit words to write the interleaved nb_channels channels
            bufpts            = npts / 2;
            rec_freq          = tdt.RZ.GetSFreq;
            chan              = tdt.nb_channels;
            stim_duration     = sweep.stimulus_duration;
            rec_duration      = sweep.record_duration;
            signal_left       = sweep.signal_left;
            signal_right      = sweep.signal_right;
            stimdelay         = sweep.stimulus_delay;
            attenPA5_left     = sweep.attenPA5_left;
            attenPA5_right    = sweep.attenPA5_right;
            attenRX8_left     = sweep.attenRX8_left;
            attenRX8_right    = sweep.attenRX8_right;
         
            try
                e.rec_dur         = tdt.RZ.SetTagVal('rec_dur',         rec_duration);  % in ms
                e.stim_dur        = tdt.RX.SetTagVal('stim_dur',        stim_duration); % in ms
                e.stim_delay      = tdt.RX.SetTagVal('stim_delay',      stimdelay);     % in ms
                e.attenRX8_left   = tdt.RX.SetTagVal('bw_atten_L',      attenRX8_left);
                e.attenRX8_right  = tdt.RX.SetTagVal('bw_atten_R',      attenRX8_right);
                e.signal_Channel  = tdt.RX.SetTagVal('signal_Channel',  tdt.opch); % kept from original circuit in case we want to modify later
                e.attenPA5_left   = tdt.PA_left.SetAtten(attenPA5_left);
                e.attenPA5_right  = tdt.PA_right.SetAtten(attenPA5_right);
                e.datain_left     = tdt.RX.WriteTagV('datain_left',  0, signal_left); % Sample selected in RPvdsEx/Implement/Device Setup.
                e.datain_right    = tdt.RX.WriteTagV('datain_right', 0, signal_right);
                e.getAttenScalarL = tdt.RX.GetTagVal('checkScalarL'); % to check bw_atten_L acted as should
                e.getAttenScalarR = tdt.RX.GetTagVal('checkScalarR');
            catch ME
                disp(ME);
                error('Make sure all the material is on and connected, and try reloading the circuits');
            end
            
            % Number of time the buffer will be (at least partly) filled up
            totalNumberSamples = (rec_duration/1000)*rec_freq*chan;
            nbHalfBuf = ceil(totalNumberSamples/bufpts); % mult. by 2 if we want the number of half-buffers to use
            
            % Number of 32bits float to load from buffer in the last cycle
            sizeLastCycle = ceil((totalNumberSamples - (nbHalfBuf-1)*bufpts)/chan)*chan; % need the 'ceil' trick because freq is not integer
            
            % Allocate memory in the PC RAM to speed up the downloading of half RAMs
            rec = cell(nbHalfBuf, 1);
            
            %recArray = zeros(1,bufpts*nbHalfBuf); % Allocating a big array now to make sure there is enough contiguous memory
            % If not, keep using rec (non-contiguous memory is fine)
            for bb=1:nbHalfBuf
                rec{bb}=zeros(1,bufpts);
            end
            
            % Check the RZ2 is ready, otherwise wait a bit for it
            ADactive = tdt.RZ.GetTagVal('ADactive');
            tAD = tic;
            hasPrinted = 0;
            while ADactive
                tocAD = toc(tAD);
                if tocAD > 10
                    error('Have been waiting for more than 10s, someting seems fishy.');
                end
                if ~hasPrinted
                    fprintf('ContinuousAcquire: Waiting for ADactive to be 0... (error after 10s)\n');
                    hasPrinted = 1;
                end
                ADactive = tdt.RZ.GetTagVal('ADactive');
            end
            
            
            % Resetting RX and RZ using zBus trigB
            tdt.zBus.zBusTrigB(0, 0, 20); % always return 0
            pause(0.02); % pause to get the reset effective
            
            % Check the left PA5 is set (can take a little time)
            timePA = tic;
            while abs(tdt.PA_left.GetAtten-attenPA5_left)>0.2 % should be the PA5 left
                toctimePA = toc(timePA);
                if toctimePA > 10,  error('Shouldn''t take so long, something is wrong'); end
            end
            
            % zBus trigger A to start the operations (play sound, record channels)
            tdt.zBus.zBusTrigA(0, 0, 20);
            ti_startPlaying = tic;
            ti_halfBufferFull=tic;
            
            % Continuous acquisition in the PC's RAM
            for bufCyc=1:nbHalfBuf
                
                to_startPlaying = toc(ti_startPlaying);
                
                curindex = tdt.RZ.GetTagVal('index');
                if fprintopt
                    disp(['  Current index: ' num2str(curindex)]);
                end
                
                % Determines which half of the buffer is being filled (1 or 2)
                which_half = 1 + mod(bufCyc-1,2);
                
                % Wait until current half of Buffer is full
                switch which_half
                    case 1
                        
                        firstInd = 0;
                        % Wait until first half of buffer is full, then index goes above bufpts
                        cond = @(curindex,bufpts)(mod(curindex,npts)<bufpts);
                        
                    case 2
                        
                        firstInd = bufpts;
                        % Wait until second half of buffer is full, then
                        % index goes back to 0 (modulo npts)
                        cond = @(curindex,bufpts)(mod(curindex,npts)>bufpts);
                end
                
                while cond(curindex,bufpts) && (to_startPlaying<rec_duration/1000)
                    curindex = tdt.RZ.GetTagVal('index');
                    to_startPlaying = toc(ti_startPlaying);
                end
                        
                to_halfBufferFull = toc(ti_halfBufferFull); % this buffer half (or less for last cycle) just got full
                ti_halfBufferFull = tic; % The (next) other half just started to be filled in now (except for the last cycle)
                
                % The last cycle may not use the whole half buffer!
                if  bufCyc<nbHalfBuf
                    nbPtsRead = bufpts;
                else
                    nbPtsRead = sizeLastCycle;
                end
                
                % Read first segment (ReadTagV) while second being written
                ti_readHalfBuffer=tic;
                for kread = 1:5
                    % For unknown reasons, this might fail sometimes, giving
                    % back a NaN. While this is possibly a very big problem, we
                    % workaround for now by using a loop.
                    % If you read this message, it means I failed. Humanity is
                    % lost. Be careful. Good luck.
                    % Update: Apparently it was solved. Leaving the msg in
                    % case
                    a = tdt.RZ.ReadTagV('dataout', firstInd, nbPtsRead); % Read from the buffer
                    if ~any(isnan(a))
                        break;
                    end
                    fprintf(' !!!  A NaN was obtained from RZ. Be careful...\n');
                end
                rec{bufCyc,1} = a; 
                to_readHalfBuffer = toc(ti_readHalfBuffer);
                
                if fprintopt
                    fprintf('  (%d,%d) firstInd=%d, nbPtsRead=%d,            recIsNaN=%d\n',...
                        bufCyc, which_half, firstInd, nbPtsRead, any(isnan(rec{bufCyc,1})));
                    
                    % Checking things are fine: reading time should be lower than filling
                    fprintf('  (%d,%d) Time for half buffer to fill: %f\n',bufCyc,which_half,to_halfBufferFull);
                    %0.78s % went up to 1.579 when decreasing the sample rate of PZ2 device from
                    % 50 kHz to 25kHz (which is the sample rate required, right?)
                    fprintf('  (%d,%d) Time to read half buffer to fill: %f\n',bufCyc,which_half,to_readHalfBuffer);
                    %1.27s % Went down to 1.045 when decreasing sample rate to 25kHz!!!
                end
                % Check to see if the data transfer rate is fast enough
                check = to_halfBufferFull < to_readHalfBuffer;
                if check
                    disp('!!!  Transfer rate is too slow (unless this was the last one)');
                end
                
            end
                
        end
        
        function [fig, pureTones, activ] = playPureTones(tdt,varargin)
            % Will play an sequence of pure tones, at various frequencies
            % and sound levels. The input is expected in the form of of a
            % series of 'Property'Value' pairs:
            %   'freq' with an array of frequencies  (in Hz)
            %   'loud' with an array of sound levels (in dB)
            %   'nbRep' with an integer, number of time to play each
            %   'saveMat' and 'saveFig' with a boolean, to save the data
            %      generated and/or the figure showing the receptive fields
            %
            % The pure tones are played in a pseudo-random sequence
            %    >> RandStream('mt19937ar', 'Seed', 200);
            %
            % A SweepTDT object is made by the method playShortBursts each
            % time, then sent to continuousAcquire.
            %
            % The activity is measured as the L1 norm of each channel
            % between tdt.minTimeActivity and tdt.maxTimeActivity (in ms),
            % after applying the filter
            %   >> filterTDT = @(x)filter(tdt.channels_filt_b, tdt.channels_filt_a,x);
            % on each channel.
            tdt.checkAndSetHeadAmp;
            
            %guidata_tdt = guidata(tdt.figure);
            %set(guidata_tdt.radioPlay, 'Value', 1);
            
            % Play a series of pure tones,
            nbRep = 2; % Default : 2 repetitions
            saveMat = 1;
            saveFig = 1;
            for ind = 1:2:(nargin-1)
                arg = varargin{ind};
                switch arg
                    case 'freq', freqInit = varargin{ind+1};
                    case 'loud', loudInit = varargin{ind+1};
                    case 'nbRep', nbRep = varargin{ind+1};
                    case 'saveMat', saveMat = varargin{ind+1};
                    case 'saveFig', saveFig = varargin{ind+1};
                end
            end
            
            % Use a pseudo-random sequence
            stream = RandStream('mt19937ar', 'Seed', 200);
            permLoud = randperm(stream, length(loudInit));
            permFreq = randperm(stream, length(freqInit));
            loud = loudInit(permLoud);
            freq = freqInit(permFreq);            
            % contain the spike times
            activ = cell(length(loud),length(freq));
            % contains the average activity
            mean_activ_1 = zeros(length(loud),length(freq), tdt.nb_channels);
            
            cTitre = sprintf('Receptive Fields: freqs=%sHz, loud=%sdB',...
                sprintf('%.1f ', freqInit(:)),...
                sprintf('%.1f ', loudInit(:)));
            [fig, ~, nl, nc] = tdt.newFigure(cTitre);            
            
            % Array of channels handles, to set their CData easily using
            %       >> findall(h{chanNum}.Children, 'Type', 'Image')
            h = arrayfun(@(chanNum)subplot(nl, nc,chanNum), 1:tdt.nb_channels, 'UniformOutput', false);
            
            % Create the figure with all subplots (takes time)
            fprintf(1,'Creating the window (takes about 10s)');
            
            % Check whether the fields are given in d.
            % If this is initialised, the positionning of the receptive
            % fields is defined by the position of the channels in the
            % electrodes
            d = findall(0,'Tag','figureGUISetupSystem3');
            if ~isempty(d) && isfield(d.UserData,'electrodesConfigs')
                configs = d.UserData.electrodesConfigs;
                cInd = 1;
                Xpos = ones(128,1);
                Ypos = ones(128,1);
                Wpos = ones(128,1); % width
                Hpos = ones(128,1); % height
                numElec = length(configs);
                
                % Checking all configs to optimise style of plot (should be automatise later)
                allElecFiles = arrayfun(@(num)fullfile(d.UserData.electrodesConfigFolder, ...
                    d.UserData.electrodesConfigList{configs(num)}), 1:numElec, 'Unif', false);
                allMat   = cellfun(@(cFile)dlmread(cFile, ' '), allElecFiles, 'Unif', false);
                allSizes =  cellfun(@(mat)[size(mat,1) size(mat,2)],  allMat, 'Unif', false);
                if numElec==1 && allSizes{1}(1) == 16  && allSizes{1}(2) == 1
                    dispo = 'dispo16x1';
                elseif numElec==2 && allSizes{1}(1) == 16  && allSizes{1}(2) == 1 &&...
                        allSizes{2}(1) == 16  && allSizes{2}(2) == 1
                    dispo = 'dispo16x1_16x1';
                else
                    dispo = 'dispo_other';
                end
                
                for confInd = 1:numElec
                    %cFile = fullfile(d.UserData.electrodesConfigFolder, ...
                    %    d.UserData.electrodesConfigList{configs(confInd)});
                    %mat = dlmread(cFile, ' ');
                    %[cSizeX, cSizeY] = size(mat);
                    cFile = allElecFiles{confInd};
                    mat = allMat{confInd};
                    cSizeX = allSizes{confInd}(1);
                    cSizeY = allSizes{confInd}(2);
                    
                    % Manually change things for 16 lines: not readable
                    % enough
                    infoTitle = ''; %#ok
                    if cSizeX == 16
                        mat = reshape(mat, cSizeX/2, cSizeY*2);
                        [cSizeX, cSizeY] = size(mat);
                        infoTitle = 'Column split in 2'; %#ok
                    end
                    
                    % Find the positions
                    try
                        for jj=1:numel(mat)
                            switch dispo
                                case 'dispo16x1'
                                    % Fine for 1 nexus A1x16 (after reshaping):
                                    [cY, cX] = find(mat==jj);
                                    Xpos(cInd) =  0.05  +        (cX-1)/((cSizeY+0.48)*numElec) + (confInd-1)/numElec;
                                    Ypos(cInd) = -0.075 + (cSizeX-cY+1)/((cSizeX+0.50)*numElec);
                                    Wpos(cInd) = 1/(numElec*cSizeY*1.25);
                                    Hpos(cInd) = 1/(numElec*cSizeX*1.10);
                                case 'dispo16x1_16x1'
                                    % Fine for 2 nexus A1x16 (after reshaping):
                                    [cY, cX] = find(mat==jj);
                                    Xpos(cInd) =  0.05  +        (cX-1)/((cSizeY+0.5)*numElec) + (confInd-1)/numElec;
                                    Ypos(cInd) = -0.045 + (cSizeX-cY+1)/((cSizeX+0.5)*1);
                                    Wpos(cInd) = 1/(numElec*cSizeY*1.25);
                                    Hpos(cInd) = 1/(numElec*cSizeX*0.70);
                                case 'dispo_other'
                                    % Fine for 2 nexus A8x8:
                                    [cY, cX] = find(mat==jj);
                                    Xpos(cInd) = -0.025 + (cX)/((cSizeX+1.5)*numElec) + (confInd-1)/numElec;
                                    Ypos(cInd) = -0.075 + (cSizeY-cY+1)/(cSizeY+0.5);
                                    Wpos(cInd) = 1/(numElec*cSizeX*1.25);
                                    Hpos(cInd) = 1/(numElec*cSizeY*0.70);
                            end
                            cInd = cInd+1;
                        end
                    catch ME
                        disp(ME)
                        error('Problem with the configuration %s, check it\n',cFile);
                    end
                end
                totalNbChannels = cInd - 1;
                Xpos = Xpos(1:totalNbChannels);
                Ypos = Ypos(1:totalNbChannels);
                Wpos = Wpos(1:totalNbChannels);
                Hpos = Hpos(1:totalNbChannels);
            end
            
            % Declare the images, invisible for now
            for ll=1:nl
                fprintf(1,'.');
                for cc=1:nc
                    % for chanNum=1:tdt.nb_channels % Need to make the children 'Image' instead of 'Line'
                    chanNum =(cc-1)*nl + ll;
                    %chanNum = (ll-1)*nl + cc;
                    imagesc( mean_activ_1(:,:,chanNum), 'Parent', h{chanNum});
                    if exist('configs', 'var')
                        pos = [Xpos(chanNum)  Ypos(chanNum)  Wpos(chanNum)  Hpos(chanNum)];
                    else
                        pos = [(ll-0.5)/(nl+1) (cc-0.5)/(nc+1) 1/(nl+1) 1/(nc+1)];
                    end
                    set(h{chanNum},'Position', pos);
                    text(0.75,1,num2str(chanNum), 'Color', [1 1 1],'Parent', h{chanNum});
                end
            end
            fprintf(1, ' done!\n');
            set(fig, 'Visible', 'on','Units','normalized',  'Position', [0 1/32 0.999 29.5/32]);
            set(0, 'CurrentFigure', fig);
            figure(tdt.figure); %#ok
            % Get bottom left child
            [~,minPosChan]  = min(cellfun(@(x)x.Position(1), h));
            if length(minPosChan)>1
                % If more than one, minimize Y
                [~,minPosChan] = min(cellfun(@(x)x.Position(2), h(minPosChan)));
            end
            childBottomLeft = h{minPosChan};
            text( nc/2-2,   -0.3, ['X-axis: Frequencies = ' sprintf('%.1f ', freqInit(:)) 'Hz'], 'Parent',childBottomLeft, 'Units','normalized');
            text(-2,        -1.8, ['Y-axis: Loudness  = '   sprintf('%.1f ', loudInit(:)) ' dB SPL'], 'Parent',childBottomLeft, 'Rotation', 90 );
            if exist('numElec', 'var')
                allXs = cellfun(@(x)x.Position(1), h);
                allYs = cellfun(@(x)x.Position(2), h);
                topletfchild  = (max(allYs)-allYs<0.02 ) & (allXs-min(allXs)<0.01);
                toprightchild = (max(allYs)-allYs<0.02 ) & (max(allXs)-allXs<0.01);
                if max(sum(topletfchild),sum(toprightchild)) > 1
                    warning('There is a problem in plotting the electrode''s values. Skipped.');
                else
                    postopleftX  = h{topletfchild}.Position(1);
                    postoprightX = h{toprightchild}.Position(1);
                    % Chosen manually..
                    for ll=1:numElec
                        text( (((numElec-ll+1)*postopleftX+ (ll-1)*postoprightX)/1 +(postopleftX+postoprightX)/1.5)*cSizeX,...
                            1.35,  sprintf('Electrode %d ',ll), 'Parent', h{topletfchild}, 'Units','normalized');
                    end
                end
            end
            
            % Set gui on top
            figure(tdt.figure);  %#ok
            
            % True for surf, false for imagesc (for refreshing the plot)
            existZData = isprop(h{1},'ZData');
            
            % Set filter properties
            if isempty(tdt.channels_filt_b)
                error('tdt.channels_filt_b should have been calculated when loading the calibration');
            end
            b = tdt.channels_filt_b;
            a = tdt.channels_filt_a;
            filterTDT = @(x)filter(b,a,x);
            
            % Start the loop
            for rep=1:nbRep % repetition
                if tdt.stopTheLoop, break; end
                s = RandStream('mt19937ar', 'Seed', rep);
                randMat = reshape(randperm(s, length(loud)*length(freq)),length(loud),length(freq));
                for ii = 1:numel(randMat)
                    [l_ind, f_ind] = find(randMat == ii);
                    if tdt.stopTheLoop, break; end
                    
                    f = freq(f_ind);
                    l = loud(l_ind);
                    durationBurst = 100; %ms
                    [rec,tdt.sweep] = tdt.playShortBursts('type', 'pureTone', 'freq', f, 'loud', l, ... %#ok
                        'burstLength', durationBurst, 'silenceLength', 200, 'repetition', 1);
                    recChannel = arrangeArray(tdt,rec);
                    tdt.sweep.text = sprintf('Freq = %.2fHz at %.2fdB SPL', f,l);
                    %activ{l_ind,f_ind} = recChannel; % gets wayyy too big
                    
                    % Apply high-pass fileter to get rid of noises
                    recChannel = cellfun(@(x)filterTDT(x),recChannel,'UniformOutput', false);
                    
                    % Setting for this run
                    tdt.minTimeActivity = tdt.stimulus_delay;
                    tdt.maxTimeActivity = tdt.stimulus_delay + durationBurst;
                    
                    % indicesAtWhichWeCalculateActivity
                    indCalcActMin = floor(max(tdt.minTimeActivity/1000*tdt.RZ_SamplingFreq, 1));
                    indCalcActMax = floor(min(tdt.maxTimeActivity/1000*tdt.RZ_SamplingFreq, length(recChannel{1})));
                    arrayWithOnes = zeros(size(recChannel{1}));
                    arrayWithOnes(indCalcActMin:indCalcActMax) = 1;
                    
                    % The activity is defined as this. Feel free to change
                    activity = @(x)norm(x.*arrayWithOnes,1);
                    
                    % Accumulate the activity in matrix mean_activ_1
                    mean_activ_1(permLoud(l_ind),permFreq(f_ind),:) = ...
                        mean_activ_1(permLoud(l_ind),permFreq(f_ind),:) + ...
                        reshape(cellfun(@(x)activity(x),recChannel),...
                        1,1,tdt.nb_channels); %cellfun(@mean,recChannel).^2;
                    
                    % Refresh the figure
                    guideChannelsPlot('pushbutton1_Callback', tdt.figure,[],...
                        guidata(tdt.figure), tdt, recChannel, tdt.sweep);
                    %  set(h, 'CData', mean_activ_1(:,:,chanNum));
                    
                    % Refresh all plots
                    for chanNum=1:tdt.nb_channels
                        cImage = findall(h{chanNum}.Children, 'Type', 'Image');
                        if length(cImage)>1
                            warning('If there''s more than one image, might need debugging...');
                        end
                        set( cImage, 'CData', mean_activ_1(:,:,chanNum));
                    end
                    
                    % Set labels
                    arrayfun(@(ax)set(ax,'YTick', 1:length(loud)),fig.Children);
                    arrayfun(@(ax)set(ax,'YTickLabel', []),fig.Children);  % loudInit)
                    arrayfun(@(ax)set(ax,'XTick', 1:length(freq)),fig.Children);
                    arrayfun(@(ax)set(ax,'XTickLabel',[]),fig.Children); % freqInit)
                    
                    % If we plot a surface instead of imagesc, refresh Z values
                    if existZData
                        set(h, 'ZData', mean_activ_1(:,:,chanNum));
                    end
                    
                end
            end
            
            % Reset tdt.stopTheLoop
            tdt.stopTheLoop = 0;
            
            % Add info into the structure pureTones
            pureTones.mean_activ = mean_activ_1;
            pureTones.freq = freqInit;
            pureTones.loud = loudInit;
            pureTones.describ = ['mean_activ_1(j,k,l): mean response of channel l to a pure '...
                'tone of frequency freq(k) played at loud(j) dbSPL (no renormalisation)'];
            if exist(tdt.RS4IP,'dir')
                dirRS4IP = dir(tdt.RS4IP);
                pureTones.RS4 =  dirRS4IP(end).name;
            end
            
            % Saving otptions, in .tmp folder if necessary
            if saveMat || saveFig
                cFol = tdt.local;
                if isempty(cFol)
                    tmpFol = fullfile(tdt.masterDirectory, '.tmp');
                    if ~exist(tmpFol, 'dir')
                        mkdir(tmpFol);
                    end
                    cFol = tmpFol;
                end
                %if ~isempty(tdt.ID)
                
                fprintf('Saving the data... ');
                if saveMat
                    warning('Saving the figure (may take 30s)\n');
                    save(fullfile(cFol, 'pureTones.mat'),'pureTones');
                end
                if saveFig
                   % rf = strrep(sprintf('receptiveFields%s', datestr(now)), ' ', '');
                    rf = strrep(sprintf('receptiveFields%s', ''), ' ', '');
                    saveas(fig, fullfile(cFol, [rf '.fig']));
                    saveas(fig, fullfile(cFol, [rf '.pdf']));
                    saveas(fig, fullfile(cFol, [rf '.jpg']));
                end
                fprintf(' done!\n');
            end
            %elseif saveMat||saveFig
            %    error('You asked to save the receptive fields, but tdt.ID is empty. This data is not saved.');
            %end
        end
        
        function [rec,sweep, e] = playShortBursts(tdt,varargin)
            tdt.checkAndSetHeadAmp();
            % Plays bursts of noise and silences for a given duration or until
            % stopped.
            %
            % Default values (duration empty means the user has to manually stop)
            burstLength = 100; %ms
            silenceLength = 100; %ms
            repetition = []; %ms
            duration = []; %ms
            rec = []; % response
            e = 1;
            
            for ind = 1:2:(nargin-1)
                arg = varargin{ind};
                switch arg
                    case 'burstLength', burstLength = varargin{ind+1};
                    case 'silenceLength', silenceLength = varargin{ind+1};
                    case 'duration',    duration = varargin{ind+1};
                    case 'type', typeNoise = varargin{ind+1};
                    case 'freq', freq = varargin{ind+1};
                    case 'loud', loud = varargin{ind+1};
                    case 'repetition', repetition = varargin{ind+1};
                end
            end
            
            if isempty(duration)&&isempty(repetition) % Continue playing until we manually stop
                continue_playing = 1;
                repSweep = 1; % default number of repetition of sweeps
                %f = figure; %#ok
                loudness = tdt.desiredLevel; % Default signal level
                % If we use guideSetupSystem3, access to the buttons values
                %fig = findall(0, 'Tag','ExploringNoiseBursts');
                tagged = findall(0, 'Tag','loudNoiseBurst');
                numberBurstToPlay = findall(0, 'Tag', 'numberBurstToPlay');
                numberBurstPlayed = findall(0, 'Tag', 'numberBurstPlayed');
                isemptyNumberBurstToPlay = isempty(numberBurstToPlay)||strcmp(get(numberBurstToPlay,'String'),'NaN');
                
                % Prevent user from changing this value
                if ~isempty(isemptyNumberBurstToPlay)
                    set(numberBurstToPlay ,'Enable', 'Off'); 
                end
                
                isemptyNumberBurstPlayed = isempty(numberBurstPlayed);
                while continue_playing %
                    if ~isemptyNumberBurstToPlay
                        burstToPlay = str2double(get(numberBurstToPlay, 'String'));
                        if burstToPlay < 1 % if we're done, exit the while
                            break;
                        end
                    end
                    if ~isempty(tagged)
                        loudness      = eval(get(findall(0, 'Tag','loudNoiseBurst', 'Type', 'uicontrol'), 'String'));
                        burstLength   = eval(get(findall(0, 'Tag','timeBurst',      'Type', 'uicontrol'), 'String'));
                        silenceLength = eval(get(findall(0, 'Tag','timeSil',        'Type', 'uicontrol'), 'String'));
                        freq          =     (get(findall(0, 'Tag','freqOrNaN',      'Type', 'uicontrol'), 'String'));
                        % randomly choose one of the given values
                        loudness      = loudness(floor(rand*length(loudness))+1);
                        burstLength   = burstLength(floor(rand*length(burstLength))+1);
                        silenceLength = silenceLength(floor(rand*length(silenceLength))+1);
                        
                    end
                    if isempty(freq)
                        freq = NaN;
                    else
                        freq = eval(freq);
                        freq = freq(floor(rand*length(freq))+1);
                    end
                    if isnan(freq)
                        tdt.sweep = tdt.make_noise_burst(burstLength, silenceLength, loudness);
                    else
                        tdt.sweep = tdt.make_pure_tone(burstLength,silenceLength,repSweep,freq,loudness);
                    end
                    [rec,e] = tdt.continuousAcquire(tdt.sweep);
                    recChannel = arrangeArray(tdt,rec);
                    guideChannelsPlot('pushbutton1_Callback', tdt.figure,[],...
                        guidata(tdt.figure), tdt, recChannel, tdt.sweep);
                    if ~isemptyNumberBurstToPlay
                        set(numberBurstToPlay, 'String',num2str(burstToPlay-1));
                        
                    end
                    if ~isemptyNumberBurstPlayed
                        burstPlayed = str2double(get(numberBurstPlayed, 'String'));
                        set(numberBurstPlayed, 'String',num2str(burstPlayed+1));
                    end
                    continue_playing = ~tdt.stopTheLoop;
                    while tdt.pauseTheLoop
                        drawnow;
                    end
                end
                tdt.stopTheLoop = 0;
                if ~isemptyNumberBurstToPlay
                    set(numberBurstToPlay ,'Enable', 'On'); 
                end
            else % Play a fixed number of bursts/pure tones
                switch typeNoise
                    case 'whiteNoise', sweep = tdt.make_noise_burst(burstLength, silenceLength, duration);
                    case 'pureTone',   sweep = tdt.make_pure_tone(burstLength,silenceLength,repetition,freq,loud);
                    otherwise, error('The type of sound (%s) is not implemented', typeNoise);
                end
                tdt.sweep = sweep;
                [rec, e] = tdt.continuousAcquire(tdt.sweep); 
            end
            if any(structfun(@not,e))
                disp('There seems there was a problem: e=');
                disp(e);
            end
        end
        
        function   sweep = make_noise_burst(tdt, burstLength, silenceLength,loudness, varargin)
            % Make a sweep or more of white noise
            if ~isempty(varargin)
                duration = varargin{1}/1000;
                numBurstSil = ceil(duration/((burstLength + silenceLength)/1000));
            else
                % duration = (burstLength + silenceLength)/1000;
                numBurstSil = 1;
            end
            
            % make simple noise; gate the first and last 100 elements
            burst = randn(floor(burstLength/1000*tdt.sound_sampling_fq),1);
            burst(1:200) = burst(1:200).*exp(linspace(-10,0,200)');
            burst(end-199:end) = burst(end-199:end).*exp(flip(linspace(-10,0,200))');
            % Add a period of silence
            noise = [burst; zeros(floor(silenceLength/1000*tdt.sound_sampling_fq), 1)]*7/max(abs(burst));
            noise = repmat(noise,numBurstSil,1);
            
            sweep = SweepTDT(tdt);
            % Set loudness before setItAll to avoid recalculating sweep.setSignals_leftAndRight
            sweep.setDesiredLevel(loudness); % desiredLevel is private for robustness
            sweep.setItAll(noise);  % loads signal, renormalise, calibrate
            sweep.text = sprintf('White noise at %0.2fdB', loudness);

        end
        
        
        function  sweep = make_pure_tone(tdt,burstLength,silenceLength,repetition,freq,loud)
            % freq in Hz, loud in decibels
            
            % Make a pure tone signal
            t = linspace(0,burstLength/1000,floor(tdt.sound_sampling_fq*burstLength/1000));
            
            % Renormalisation of the sound using the calibration, based on
            % Path/to/makeRFstimConfig.m
            if isempty(tdt.calib_left)
                error('Please load the calibrations in tdt.calib_left and _right');
            end
            
            % Will be renormalised in the SweepTDT object
            pureTone = sin(2*pi*t*freq); 
            
            % Add silence (we take care of gating problem (clicking sound during discontinuities)
            if silenceLength>0
                sil = zeros(1,floor(tdt.sound_sampling_fq*silenceLength/1000));
                % Gating manually the tone if ;onf enough (should away be
                % the case given the high sampling rate)
                if length(burstLength) > 400
                    pureTone(1:200) = pureTone(1:200).*exp(linspace(-10,0,200)');
                    pureTone(end-199:end) = pureTone(end-199:end).*exp(flip(linspace(-10,0,200))');
                end
                signal_single_sweep = [pureTone sil];
            else
                signal_single_sweep = pureTone;
            end
            sig = repmat(signal_single_sweep, 1, repetition);
            
            % Attenuation based on the calibrations
            % Attenuation of a pure tone can be done using directly the
            % calibrations or their fit. For consistency, we use the fit.
            %    [~, i] = min(abs(tdt.calib_freqs - freq)); % index of closest frequency. All in Herz
            %    attenL = tdt.calib_left(i)  - loud - tdt.headampAtten;
            %    attenR = tdt.calib_right(i) - loud - tdt.headampAtten;
    
            sweep = SweepTDT(tdt);
            sweep.setDesiredLevel(loud); % setDesiredLevel is private. Set before giving signal)
            sweep.setItAll(sig);
            sweep.text = sprintf('Pure tone: dur=%.2fs, freq=%.2f, loud=%.2fdB.',...
                burstLength/1000,freq,loud); 
        end
        
        function [f,h, nl, nc] = newFigure(tdt,~) % %#ok
            % Open a new figure (slow, hence done before running System 3)
            % nb_channel may have been changed by guideSetupSystem3.setElectrodesGivenPopups
            nChan = tdt.nb_channels;
            h = zeros(nChan,1);
            switch nChan
                case 128, nl=8; nc=16;
                case 64,  nl=8; nc=8;
                case 32,  nl=4; nc=8;
                case 16,  nl=4; nc=4;
                case 8,   nl=4; nc=2;
                otherwise, error('Number of channels (%d) is probably wrong');
            end
            
            f = figure('Visible','off'); %#ok
            % title(cTitre);
            for kk=1:nChan
                subplot(nl,nc,kk);
                h(kk)=plot(0,0); % h(kk)=plot(recChannel{kk}) % kept for later
            end
            
        end
        
        %%%%%%%%%%%%%% ID 
        
        function ID = make_ID(~)
            % Make an ID used both for puretonogram and experiment. Fixed after
            % the electrode has been positionned
            ID = datestr(now);
            ID = strrep(ID, ' ', '_');
            ID = strrep(ID, ':', '_');
            ID = strrep(ID, '-', '_');
        end
        
        %%%%%%%%%%%%%%%% SAVE SOFT VERSION OF TDT
        
        function saveLightTdt(tdt,filename, varargin)
            % Save the current tdt variable at filename, without a few things
            % (figure, loaded calibration...) given as fields in varargin
            % ex:  tdt.saveLightTdt(fullfile(tdt.local,'tdt.mat'), 'figure', 'calib_left', 'calib_right')
            % If tdt too big, use -v7.3 can work but too slow. Instead, we
            % got rid of various elements of the object before saving.
            for kk=varargin
                prop = kk{1};
                if isprop(tdt, prop)
                    tdtHold.(prop) = tdt.(prop);
                    tdt.(prop) = [];
                    % Normally they are not handles, hence it does't break anything
                else 
                    fprintf('tdt.saveLightTdt: Input %s is not a property of tdt.\n',prop);
                end
            end
            save(filename,'tdt' );%;, '-v7.3'); % use another version if breaks 
            for kk=varargin
                prop = kk{1};
                if isprop(tdt, prop)
                    tdt.(prop) =  tdtHold.(prop);
                end
            end
        end
        
        %%%%%%%%%%%%%%% CHANNELS MANAGEMENT %%%%%%%%%%%%%%%
        
        function deleteFromNonDesiredChannels(tdt, pathToCurrentRecordings, fprintoptions)
        % Only keeps the channels given in tdt.arrayChanRec from RS4
        cFol = fullfile(tdt.RS4IP, pathToCurrentRecordings);
        if nargin<3
            fprintoptions = 0;
        end
        if ~exist(cFol,'dir')
            fprintf('deleteFromNonDesiredChannels: The given folder (%s) doesn''t exist.\n',pathToCurrentRecordings);
            return;
        end
        arrayToDelete = setdiff(1:128, tdt.arrayChanRec);
        dirCFol = dir(cFol);
        for kk=arrayToDelete
            % Get the file containing ['ch ' kk]
            ind = cellfun(@(s)(strfind(s,['ch' num2str(kk) '.sev'])), {dirCFol.name},'Unif', false);
            cFile = fullfile(cFol,dirCFol(~cellfun(@isempty,ind)).name);
            delete(cFile);
            if fprintoptions
                fprintf('Deleting %s.\n',cFile);
            end
        end
        end
        
        %%%%%%%%%%%%%%%% CALIBRATION %%%%%%%%%%%%%%%%%%%%%%
     
        function load_calibration(tdt, calibLeft, calibRight)
            % Call load_calibration_single on left and right ear. Check
            % they are done at the same frequencies.
            load_calibration_single(tdt, calibLeft, 'left')
            load_calibration_single(tdt, calibRight, 'right')
             if ~isequal( calibL.freqs,  calibR.freqs)
                error('Calibration frequencies are supposed to be equal');
            else
                fprintf('Calibrations succesfully loaded\n');
            end
        end
        
        function load_calibration_single(tdt, calibPath, side)
            % Load calibFilename into tdt.calib_left, _right, _freqs,
            % calculate  tdt.polyfit to be used later with polyval
            calib = load(calibPath);
            switch side
                case 'left'
                    tdt.calib_left_file = calibPath;
                    tdt.calib_freqs = calib.freqs;
                    tdt.calib_left  = calib.avg_corr_rep;
                    %tdt.calib_polyfit_left        = polyfit(tdt.calib_freqs, tdt.calib_left,3);
                    tdt.calib_polyfit_left_log10 =  tdt.fitCalibration(tdt.calib_freqs, tdt.calib_left);
                case 'right'
                    
                    tdt.calib_right_file = calibPath;
                    tdt.calib_freqs = calib.freqs;
                    % Should contain fields freqs and avg_corr_rep
                    tdt.calib_right = calib.avg_corr_rep ;
                    tdt.calib_freqs = calib.freqs;
                    %tdt.calib_polyfit_right = polyfit(tdt.calib_freqs, tdt.calib_right,3);
                    tdt.calib_polyfit_right_log10 =  tdt.fitCalibration(tdt.calib_freqs, tdt.calib_right);
                otherwise
                    error('You choose the wrong side.')
            end
           
        end
        
        
        function cFit_log10 = fitCalibration(tdt, freqs, calib)
            % Uses polyfit to create a polynomial fitting on the
            % log-scale calibration
            
            Fs = tdt.sound_sampling_fq;
            % freqToUse = sweep.tdt.calib_freqs
            % calibToUse = sweep.tdt.(['calib_' side]) ;
            f_minToUse = 100;   %Hz, decided by Trevor on Juy 12th 2016
            f_maxToUse = 12000; %Hz. 14000 was too much, 10000 too low
            indToUse   = ...
                freqs < Fs/2        & ...
                freqs > f_minToUse  & ...
                freqs < f_maxToUse;
            freqToUse  = freqs( indToUse );
            calibToUse = calib( indToUse );
            
            % Log-spread them
            %freqToUse  = logspace(log10(freqToUse(1)), log10(freqToUse(end)), length(freqToUse))';
            indToUse = floor([1:9 logspace(1.05, log10(length(freqToUse)),60)])';
            freqToUse  = freqToUse(indToUse);
            calibToUse = calibToUse(indToUse);
            
            % Make it flat under f_min and above f_max (Trevor's choice)
            freq_before_fMin = 1:10:f_minToUse-1;
            freq_after_fMax = logspace(log10(f_maxToUse+1),log10(Fs/2), 30);
            freqToUse  = [freq_before_fMin   freqToUse'   freq_after_fMax];
            
            %calib_before_fMin = linspace(120, calibToUse(1), length(freq_before_fMin)); %     
            calib_before_fMin = calibToUse(1) * ones(size(freq_before_fMin)) ;
            calib_after_fMax  = calibToUse(end)*ones(size(freq_after_fMax));
            calibToUse = [calib_before_fMin    calibToUse'    calib_after_fMax];
            
            % Trevor's call
            typeFit = 'spline';
            switch typeFit
                case 'polynomial'
                    % Fitting a polynomial of degree 'deg'
                    %warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
                    deg = 7;
                    cFit_log10 = polyfit(log10(freqToUse), calibToUse, deg);
                    %warning('on','MATLAB:polyfit:RepeatedPointsOrRescale')
                    
                case 'spline'
                    % Fitting splines 
                    cFit_log10 = @(f)spline(log10(freqToUse), calibToUse, log10(f));
                    %figure; semilogx(freqToUse, calibToUse, (f), cFit_spline_log10);
                    
                otherwise
                    error('Type %s not implemented', typeFit);
            end
        end
        
        
        
        
        %%%%%%%%%%%%%%%%  CHECKING %%%%%%%%%%%%%%%%%%%%%%%%
    
        
        function tdt = checking_length_wav_files(tdt)
            % Checks the wav files so that they all fit within the buffer defined in
            % the circuit. Outputs 1 if they are all fine, 0 otherwise.
            
            % Check the buffer can load all the wav file. For some reason, doesn't work
            % with RX but does with RZ... Maybe because it's fixed, but dynamically
            % changeable in RZ?
            bufferSize = tdt.RX.GetTagVal('sizeBuf');
            if bufferSize==0
                % disp('We can''t read buffer size because it''s a static parameter. Assume 4e7 (check manually).');
                bufferSize = 1e7; % In RZ2_Manual, p.1-12, Memory = 64MB SDRAM per DSP, I understand this is the max
                % But if I try more than 1e7, the circuit won't load
                % (Cirecuit_Loader_device.m error when loading)
            end
            
            arrayLengths = zeros(size(tdt.listWavFiles));
            for kk=1:length(tdt.listWavFiles)
                switch class(tdt.listWavFiles)
                    case 'cell',
                        wavPath = tdt.listWavFiles{kk};
                        %info = audioinfo(wavPath); % audioinfo in later versions, mmfileinfo in old
                        info = mmfileinfo(wavPath);% for some reason, audioinfo not working...
                        Duration  =  info.Duration;
                    case 'struct',
                        Duration = tdt.listWavFiles(kk).duration;
                    otherwise,     error('Class of listWavFiles not recognised');
                end
                
                % length of sound after resampling
                arrayLengths(kk) = Duration * tdt.sound_sampling_fq;
            end
            
            % Check it's all good
            isbiggerThanBuffer = arrayLengths > bufferSize;
            out = any(isbiggerThanBuffer);
            if out
                fprintf('The following files won''t fit into memory (bufferSiz=%d normally).\n', bufferSize);
                fprintf('Cut them, increase the RAM memory, or write dynamically, as in Continuous_Play.m\n');
                %fprintf('%s\n', obj.listWavFiles.path{isbiggerThanBuffer});
                fprintf('%d ', find(isbiggerThanBuffer));
                error('\nChange the wav file before continuing.');
            end
            
        end
        
        function check_numberCores(~)
            if feature('numCores')<2
                error('According to tests run manually, we need at least two CPUs to go as fast as the TDT buffer loading');
            end
        end
        
        function check_sizeBufferMultipleOfNbChans(obj)
            if floor(obj.sizeBuf/(obj.nb_channels*2))~=obj.sizeBuf/(obj.nb_channels*2)
                error('We require that half of the buffer be a multiple of the number of channels, for simplicity.');
            end
        end
        
        function check_RZSamplingFreq(tdt)
            if tdt.RZ_SamplingFreq ~= floor(tdt.RZ.GetSFreq);
                % could be 97656.25 or 48825 % Chosen to be 24414 for RZ
                error('RZ has a sampling freqeuency of %d, instead of %d',floor(tdt.RZ.GetSFreq),obj.RZ_SamplingFreq);
            end
        end
        
    end
    
end