function varargout = guideSetupSystem3(varargin)
%GUIDESETUPSYSTEM3 M-file for guideSetupSystem3.fig
%      GUIDESETUPSYSTEM3, by itself, creates a new GUIDESETUPSYSTEM3 or raises the existing
%      singleton*.
%
%      H = GUIDESETUPSYSTEM3 returns the handle to a new GUIDESETUPSYSTEM3 or the handle to
%      the existing singleton*.
%
%      GUIDESETUPSYSTEM3('Property','Value',...) creates a new GUIDESETUPSYSTEM3 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to guideSetupSystem3_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUIDESETUPSYSTEM3('CALLBACK') and GUIDESETUPSYSTEM3('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUIDESETUPSYSTEM3.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guideSetupSystem3

% Last Modified by GUIDE v2.5 12-Jul-2016 18:00:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guideSetupSystem3_OpeningFcn, ...
                   'gui_OutputFcn',  @guideSetupSystem3_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% Alban's dev comment, please don't look: 
% optional THINGS TO FINISH 
%   - solve the datenum problem: a device is not well setup, probably RZ2,
%   leading to a wrong date of file creaton. This shouldn't be a problem as
%   we pick the last created file, but still...
% - The edit box focus: a trick was used to make the comments disappear
% after clicking 'enter', but the focus is then automatically set on the
% edit box.
% - Automatising the headamp setting: Would prefer to have it automatically
% set to 20 when we play a sound (except calibration)
% - Generate pdf after the experiment is finished, using all the comments and
% plots
% - Packaging: put all in a Matlab package to import. For some reason, it
% broke...
% - Move all the file from RS4 after an experiment
% - get rid of the ugly bottom left bit in guideChannelsPlot figure
% - check levels with devices
% - Fit in log?
% Add logo? (tried, was ugly)
% Add comment about saving plot of fitting calibration, and alternative
% when tdt.local hasn't been created yet (saved in tdt.masterdirectory/.tmp,
%  then moved the .pdf and .fig to tdt.local when ID is made)
% Write comments even out of experiment? (pure tones,...)


% --- Executes just before guideSetupSystem3 is made visible.
function guideSetupSystem3_OpeningFcn(hObject, event, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for guideSetupSystem3
handles.output = hObject;

%[hObject, event, handles] =
initialiseTdt(hObject, event, handles, varargin);
set(handles.radioLazy, 'Value' ,1);


% Update handles structure
guidata(hObject, handles);


function initialiseTdt(hObject, event, handles, varargin)

% Code to execute when the GUI is open
pathGuideSetup = get(handles.figureGUISetupSystem3, 'Filename');
cLoc = fileparts(pathGuideSetup); %
addpath(fullfile(cLoc)); 

v = regexp(version, '(\w*\)', 'match');
if ~strcmp(v, '(R2015a)')
    % NOT USED AT THE MOMENT
%     warning(['On the IHR''s physiology lab computer, we require 2015a to avoid hardware\n' ...);
%     'bug and use the //processing toolbox (Windows 32-bit not supported in later releases)']);
end


% Object of class AlbanTDT used through the whole process. Run 
% 'doc AlbanTdt' to know more
tdt = AlbanTDT();
tdt.masterDirectory = cLoc;

tdt.circuitPathReco = fullfile(cLoc,'RPvdsEx','ContinuousRecord_128channels.rcx'); % Without CoreSweepControl (OpenEx thing)
tdt.circuitPathPlay = fullfile(cLoc,'RPvdsEx','ContinuousPlay_zBusA.rcx'); % two signals, one for each ear

tdt.connect('RX8', 'RZ2', 'ZBUS', 'PA_left', 'PA_right');

% Check a few things
tdt.check_numberCores;
tdt.check_sizeBufferMultipleOfNbChans;
tdt.check_RZSamplingFreq;

% Name of all variables to load manually
varToLoad = {'calib_left_file', 'calib_right_file', 'depth', 'RS4IP', 'ID'}; 

% Save in UserData so that all functions and objects can access it
handles.figureGUISetupSystem3.UserData.tdt =  tdt;
handles.figureGUISetupSystem3.UserData.variablesToLoad =  varToLoad;
handles.figureGUISetupSystem3.UserData.savePlotCalibsFit = 0;

% Add logo to the GUI; ugly when tried
% logofile = fullfile(fileparts(get(handles.figureGUISetupSystem3, 'Filename')),'Logos', 'MRC_IHR_Nottingham_transp_small2.jpg');
% rgbImage = imread(logofile);
% imshow(rgbImage,'Parent', handles.axesLogoIhr);
% set(handles.axesLogoIhr,'visible','off') 
% set(handles.axesLogoIhr,'xtick', [])
% set(handles.axesLogoIhr,'ytick', []);

% To make the GUI disabled while working with TDT (prperties 'CreateFcn' 
% and 'CloseRequestFcn'of all dialog boxes)
freezeFigureGUISetupSystem3   = 'set(findobj(findall(0,''Tag'',''figureGUISetupSystem3''),''Enable'',''On''),  ''Enable'', ''Off'')';
unfreezeFigureGUISetupSystem3 = 'set(findobj(findall(0,''Tag'',''figureGUISetupSystem3''),''Enable'',''Off''), ''Enable'', ''On''); set(findall(0,''Tag'',''radioLazy''),''Value'',1);';
handles.figureGUISetupSystem3.UserData.freezeFigureGUISetupSystem3    = freezeFigureGUISetupSystem3;
handles.figureGUISetupSystem3.UserData.unfreezeFigureGUISetupSystem3  = unfreezeFigureGUISetupSystem3;
handles.figureGUISetupSystem3.UserData.CreateFcn = freezeFigureGUISetupSystem3;
handles.figureGUISetupSystem3.UserData.CloseRequestFcn =  ['closereq; ' unfreezeFigureGUISetupSystem3];

% Setup the configurations: 2 neuroNexusA8x8 by default (128 channels)
Selectconfigurationforeachelectrode(hObject, event, handles);



% --- Outputs from this function are returned to the command line.
function varargout = guideSetupSystem3_OutputFcn(hObject, event, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushLoadCalibLeft.
function pushLoadCalibLeft_Callback(hObject, event, handles)
% hObject    handle to pushLoadCalibLeft (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, filterIndex] = goGetCalibFileLoad(hObject, event, handles);



% --- Executes on button press in pushLoadCalibRight.
function pushLoadCalibRight_Callback(hObject, event, handles)
% hObject    handle to pushLoadCalibRight (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, filterIndex] = goGetCalibFileLoad(hObject, event, handles);



function [filename, pathname, filterIndex] = goGetCalibFileLoad(hObject, event, handles, cPath)
% filename optional argument if we already know the file (loaded after calibrating)
tdt = handles.figureGUISetupSystem3.UserData.tdt;
if ~exist('cPath', 'var') ||isempty(cPath)
    if isstruct(hObject)
        tag = hObject.Tag;
    else
        tag = get(hObject, 'Tag');
    end
    switch tag
        case 'pushLoadCalibLeft',  side = 'left';  calib = tdt.calib_left;
        case 'pushLoadCalibRight', side = 'right'; calib = tdt.calib_right;
        otherwise, error('Something unexpected: Tag not recognised');
    end
    if ~isempty(calib)
        answer = questdlg(sprintf('A %s calibration is already loaded; do you want to change it?',side),...
            sprintf('Change %s calibration?',side), 'Yes', 'No', 'Yes');
        if strcmp(answer, 'No')
            return;
        else
            tdt.(['calib_' side]) = [];
            set(handles.(['checkboxCalib' upper(side(1)) side(2:end)]), 'Value', 0);
        end
    end
    
    
    if isempty(tdt.local)
        cInitPath = tdt.masterDirectory;
    else
        cInitPath = tdt.local;
    end
    [filename, pathname, filterIndex] = uigetfile(cInitPath );
    if ~filterIndex
        return;
    end
    % Check it has the appropriate shape
    cPath = fullfile(pathname,filename);
end

whoscalib = whos('-file', cPath);
if ~all(ismember({whoscalib.name}, {'freqs','avg_corr_rep'}))
    warning('The file %s doesn''t contain freqs and avg_corr_rep variables, \nhence we assume there''s a problem.',cPath);
    return;
end

% Trick to use the same notations for two different uses
if isstruct(hObject)
    tag = hObject.Tag;
else
    tag = get(hObject, 'Tag');
end
switch tag
    case 'pushLoadCalibRight'
        tdt.calib_right_file = cPath;
        side = 'right';
        tdt.load_calibration_single(cPath, side);
        set(handles.checkboxCalibRight,'Value',1)
    case 'pushLoadCalibLeft'
        tdt.calib_left_file = cPath;
        side = 'left';
        tdt.load_calibration_single(cPath, side);
        set(handles.checkboxCalibLeft,'Value',1)
    otherwise
        error('Wrong hObject...');
end

% Make a plot of the fitting, save it in local folder
fig = figure; %('Visible', 'off');
Fs = tdt.sound_sampling_fq; % sampling freq of signal
f = Fs/2*linspace(0,1,5000);
cFit_log10 = tdt.(['calib_polyfit_' side '_log10']);
switch class(cFit_log10)
    case 'double', curvef = polyval(cFit_log10,log10(f)); % We used polyfit
        cfit = 'Polynomial';
    case 'function_handle', curvef = cFit_log10(f);       % We used splines
        cfit = 'Splines';
    otherwise, error('Class %s not recognised',class(cFit_log10));
end
semilogx(tdt.calib_freqs,tdt.(['calib_' side]), f, curvef);
legend('Calibration', 'Fit');
legend('Location', 'SouthWest');
xlabel('Frequencies (Hz)')
ylabel('Sound Level (dB)');
title(sprintf('Calibration %s ear, with %s fit', side, cfit));
set(fig, 'PaperPosition',[0 0 14 7]);
set(fig, 'PaperSize',[14 7]);
cFol = tdt.local;
if isempty(tdt.local)
    % Save in a temporary folder. When local is selected, empty the tmp
    % folder in it
    cFol = fullfile(tdt.masterDirectory, '.tmp');
    try
        if ~exist(cFol, 'dir')
            mkdir(cFol)
        end
    catch
        warning('Couldn''t create tmp in tdt.masterDirectory');
    end
end
% Save .fig and .pdf (only documents moved from .tmp normally)
cFile = fullfile(cFol, ['calib_' side '_' cfit ]);
if exist(cFol,'dir')
    saveas(fig, [cFile '.png']); % to add to the report (simpler than .pdf)
    saveas(fig, [cFile '.pdf']);
    saveas(fig, [cFile '.fig']);
else
    warning('Couldn''t save %s', cFile);
end
close(fig)


% --- Executes on button press in pushRS4.
function pushRS4_Callback(hObject, event, handles)
% hObject    handle to pushRS4 (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
persistent RS4IP_persi
if isempty(RS4IP_persi)
    defaultText = 'XXX.XXX.X.XX';
else 
    defaultText = RS4IP_persi;
end
RS4IP = handles.figureGUISetupSystem3.UserData.tdt.RS4IP;

if ~isempty(RS4IP)
    answer = questdlg(sprintf('The current RS4 IP in memory is %s, do you want to change it?',RS4IP),...
         'Change RS4 IP?', 'Yes', 'No', 'Yes');
    if strcmp(answer, 'No')
        return;
    else 
        handles.figureGUISetupSystem3.UserData.tdt.RS4IP = '';
        set(handles.checkboxRS4ID, 'Value', 0);
    end
end
     
in = inputdlg('Enter RS4 IP''s adress (Panel ''Status'' in RS4, has to be connected to network)',...
    'RS4 IP Adress', 1, {defaultText});
if isempty(in)
    return
end
RS4IP = ['\\' in{1} '\data'];
if exist(RS4IP, 'dir')
    handles.figureGUISetupSystem3.UserData.tdt.RS4IP = RS4IP;
    set(handles.checkboxRS4ID, 'Value', 1);
else
    set(handles.checkboxRS4ID, 'Value', 0);
    error(['The dynamic IP of RS4 is not ''%s'' (or RS4 not connected). In tab ''Status'' of the '...
        'RS4, read th IP, correct it in Matlab code'], in{1});
end
RS4IP_persi = in{1};




% --- Executes on button press in checkboxCalibLeft.
function checkboxCalibLeft_Callback(hObject, event, handles)
% hObject    handle to checkboxCalibLeft (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCalibLeft


% --- Executes on button press in checkboxCalibRight.
function checkboxCalibRight_Callback(hObject, event, handles)
% hObject    handle to checkboxCalibRight (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCalibRight


% --- Executes on button press in checkboxRS4ID.
function checkboxRS4ID_Callback(hObject, event, handles)
% hObject    handle to checkboxRS4ID (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxRS4ID



function editDepth_Callback(hObject, event, handles)
% hObject    handle to editDepth (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDepth as text
%        str2double(get(hObject,'String')) returns contents of editDepth as a double
str = get(handles.editDepth, 'String');
cs = str2double(str);
if isnan(cs)
    warning('The value set should be a double');
    set(handles.editDepth, 'String', '');
    set(handles.checkboxDepth, 'Value', 0);
else 
    handles.figureGUISetupSystem3.UserData.tdt.depth = cs;
    set(handles.checkboxDepth, 'Value', 1);
end


% --- Executes during object creation, after setting all properties.
function editDepth_CreateFcn(hObject, event, handles)
% hObject    handle to editDepth (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in checkboxDepth.
function checkboxDepth_Callback(hObject, event, handles)
% hObject    handle to checkboxDepth (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDepth


% --- Executes on key press with focus on checkboxCalibLeft and none of its controls.
function checkboxCalibLeft_KeyPressFcn(hObject, event, handles)
% hObject    handle to checkboxCalibLeft (see GCBO)
% event  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushMakeID.
function pushMakeID_Callback(hObject, event, handles)
% hObject    handle to pushMakeID (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkboxMakeID, 'Value')
   answer = questdlg('An ID is already set. Are you sure you want to create a new one?',...
       'Check Change ID', 'Yes', 'No','No');
   switch answer
       case 'No', 
           return;
   end
end

tdt = handles.figureGUISetupSystem3.UserData.tdt;
defaultDirectory = tdt.make_ID();
fprintf('Default directory would be: %s\n', defaultDirectory);

answer = questdlg(['Do you want the default ID (' defaultDirectory ') or reopen an existing one?'],...
    'Check Make ID', 'Default', 'Existing','Existing');
switch answer
    case 'Default'
        cdir = fullfile(tdt.masterDirectory, defaultDirectory);
    case 'Existing'
        cdir = uigetdir(fullfile(tdt.masterDirectory, defaultDirectory),'Choose the ID for the experiment (default ID is current time)');
        if isnumeric(cdir)&&cdir==0
            return;
        end
end


[~, tdt.ID] = fileparts(cdir);
tdt.local = cdir; % Folder for the experiment and puretone

if ~exist(tdt.local, 'dir')
    mkdir(tdt.local);
end

% Move all that was waiting in the tmp folder
tmpfolder =  fullfile(tdt.masterDirectory, '.tmp');
if exist(tmpfolder, 'dir')
    % '*' was breaking... don't put anything important here
    listFormat = {'png', 'pdf', 'fig', 'mat'};
    for kk=1:length(listFormat)
        fmt = listFormat{kk};
        try
            copyfile([tmpfolder filesep '*' fmt], tdt.local);
        catch
        end
    end        
    %rmdir(tmpfolder, 's');
end

set(handles.checkboxMakeID, 'Value', 1);


% --- Executes on button press in checkboxMakeID.
function checkboxMakeID_Callback(hObject, event, handles)
% hObject    handle to checkboxMakeID (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxMakeID


% --- Executes on button press in pushDisp.
function pushDisp_Callback(hObject, event, handles)
% hObject    handle to pushDisp (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Displaying the variables to load for System 3 experiments:');
tdt = handles.figureGUISetupSystem3.UserData.tdt;
for kk=handles.figureGUISetupSystem3.UserData.variablesToLoad
    var=kk{1};
    if isempty(tdt.(var))
        fprintf('!!!\ttdt.%s: empty\n', var);
    else
        switch class(tdt.(var))
            case 'char',    fprintf('\ttdt.%s = %s\n',var, tdt.(var));
            case 'double',  fprintf('\ttdt.%s = %.2f\n',var, tdt.(var));
            otherwise,      fprintf('\ttdt.%s = ',var); disp(tdt.(var));                 
        end
    end
end
if isfield(handles.figureGUISetupSystem3.UserData,'electrodesConfigs')
    electrodesConfigList = handles.figureGUISetupSystem3.UserData.electrodesConfigList;
    pathToAllConfig = cellfun(@(file)fullfile(tdt.masterDirectory, 'electrodeConfigurations', file),...
        electrodesConfigList, 'unif', false);
    electrodesConfigs = handles.figureGUISetupSystem3.UserData.electrodesConfigs;
    numChannelPerElec = cellfun(@(conf)numel(dlmread(conf, ' ')),pathToAllConfig);
    for kk = 1:length(electrodesConfigs)
        fprintf('\tElectrode %d (%d channels) is configured with config %s\n',kk,numChannelPerElec(electrodesConfigs(kk)),...
            electrodesConfigList{electrodesConfigs(kk)} );
    end
end

function goforit = checkVariablesAreLoaded(hObject, event, handles)
% Just to check all variables are loaded as they should 
tdt = handles.figureGUISetupSystem3.UserData.tdt;
varToLoad = handles.figureGUISetupSystem3.UserData.variablesToLoad;
goforit = 1;
for kk=1:length(handles.figureGUISetupSystem3.UserData.variablesToLoad)
    var=varToLoad{kk};
    if ~isempty(tdt.(var))
        varToLoad{kk} = '';
    end
end
if any(~cellfun(@isempty,varToLoad))
    % Play a sound
    beep;
    % Disable bactrace to avoid the trace
    warning('off','backtrace') 
    warning('There are some variables to load before doing this:\n%s\n',sprintf('\t\t%s',varToLoad{:}));
    warning('on','backtrace')
    goforit = 0;
end


% --- Executes on button press in radiobuttonRunCalib.
function radiobuttonRunCalib_Callback(hObject, event, handles)
% hObject    handle to radiobuttonRunCalib (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if there's already a dialog box open
dd = findall(0,'Tag', 'guideRunCalibrationDlg');
if ~isempty(dd)
    figure(dd);
    return;
end

% Hint: get(hObject,'Value') returns toggle state of radiobuttonRunCalib
tdt = handles.figureGUISetupSystem3.UserData.tdt;

% GUI to run the calibrations and save them
if isempty(handles.figureGUISetupSystem3.UserData.tdt.local)
    cInitPath = tdt.masterDirectory;
else 
    cInitPath = tdt.local;
end
filename_left = '';
while isempty(strfind(filename_left, 'left'))
    [filename_left, pathname, filterIndex] = uiputfile(fullfile(cInitPath,'calib_left_freq_resp.mat'),...
        'Select where to save the left calibration m-file (the ''right'' automatically replacing ''left'')');
    if ~filterIndex % pressed 'Cancel'
        set(handles.radioLazy, 'Value' ,1);
        return;
    end
    if isempty(strfind(filename_left, 'left'))
        warning('The chosen filename must contain ''left'' in it!');
    end
end
filename_left = lower(filename_left);
filename_right = strrep(filename_left, 'left', 'right');
filename_left = fullfile(pathname, filename_left);
filename_right = fullfile(pathname, filename_right);
calibdlg(handles, filename_left,filename_right)

% Return if calib_freqs is empty or window still there
if isempty(tdt.calib_freqs)||~isempty(findall(0,'Tag','guideRunCalibrationDlg'))
    return;
end

% Save a plot with the fits
% Approximate calibration at the required frequencies, using the polynomial fit
% Pfit calculated by tdt.load_calibration_single when loading a calibration
f = linspace(0,48818,2000); % Frequencies
curvef_left  = polyval(tdt.calib_polyfit_left,f);
curvef_right = polyval(tdt.calib_polyfit_right,f);
fig = figure('Visible', 'On');
subplot(2,1,1);
plot(f, curvef_left, 'r', tdt.calib_freqs, tdt.calib_left, 'b'); 
title('Left calibration');
ylabel('(in dB SPL)');
%[cPath, cFile, ~] = fileparts(filename_left);
%saveas(fig,fullfile(cPath, [cFile '_withFitCurve', '.pdf']));
subplot(2,1,2);
plot(f, curvef_right,'r', tdt.calib_freqs, tdt.calib_right, 'b'); 
title('Right calibration');
[cPath, ~, ~] = fileparts(filename_right);
xlabel('Frequencies (Hz)');
ylabel('Sound Level Output ');
saveas(fig,fullfile(cPath, ['leftRightCalibrations_withFitCurve', '.pdf']));
saveas(fig,fullfile(cPath, ['leftRightCalibrations_withFitCurve', '.fig']));
fprintf('Saving %s....\n',fullfile(cPath, ['leftRightCalibrations_withFitCurve', '.pdf']));
close(fig);

set(handles.radioLazy, 'Value' ,1);


function calibdlg(handles, filename_left,filename_right)
% Function the open a dialog box and run the calibrations

dd = findall(0,'Tag', 'figureGUISetupSystem3');
pos = get(dd, 'Position');

d = dialog('Tag', 'guideRunCalibrationDlg', ...
    'WindowStyle', 'normal',...
    'Name', 'Run ears calibrations with 120dB noise',...
    'CreateFcn', handles.figureGUISetupSystem3.UserData.CreateFcn,...
    'CloseRequestFcn',handles.figureGUISetupSystem3.UserData.CloseRequestFcn );
set(d, 'Units', get(dd,'Units'));
set(d, 'Position', [pos(1)+pos(3)/2-77/2 pos(2)+pos(4)/2 75 13])
set(d, 'PaperPosition',get(dd,'PaperPosition'))

txt = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(175-210/2) 95 210 40],...
    'String', 'Select an ear');
popup = uicontrol('Parent', d,...
    'Style', 'popup', ...
    'Position', [(175-100/2) 85 100 25],...
    'String', {'left','right'});%,...
    %'Callback', @popup_callback);
autoLoad = uicontrol('Parent', d, ...
    'Position', [(175-230/2) 10 230 40],...
    'Value',1,...
    'Style','checkbox', ...
    'String', 'Check to automatically load this calibration');
btn = uicontrol('Parent', d,...
    'Position', [(175-170/2) 50 170 25],...
    'String', 'Run calibration for this ear', ...
    'Style', 'pushbutton',...
    'Value', 0,...
    'callback', {@popup_callback, handles,popup,autoLoad,filename_left,filename_right});

% Wait until d is deleted
%uiwait(d);


function popup_callback(hObject,event,handles,popup,autoLoad,filename_left,filename_right)
% Run calibration, save according to popup choice
if get(hObject,'Value') == 0
    return;
end
switch get(popup,'Value') % {left, right}
    case 1, filename = filename_left;  tagForCalib.Tag = 'pushLoadCalibLeft';
    case 2, filename = filename_right; tagForCalib.Tag = 'pushLoadCalibRight';
end
% disp(filename); %
try
    % sets the attenuation to 00
    handles.figureGUISetupSystem3.UserData.tdt.headampAtten = 0;
    splcal_trevors(filename);
catch ME
    if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
        warning('The function splcal_trevors.m is not found. addpath the folder that contains it');
    else
        warning('splcal_trevors bugged. Solve it.');
    end
    return
end
if get(autoLoad, 'Value')
    % load calibration
     goGetCalibFileLoad(tagForCalib, event, handles, filename)
end



% --- Executes on button press in radiobuttonPlayPureTones.
function radiobuttonPlayPureTones_Callback(hObject, event, handles)
% hObject    handle to radiobuttonPlayPureTones (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
goforit = checkVariablesAreLoaded(hObject, event, handles);
if ~goforit
    set(handles.radioLazy, 'Value' ,1);
    return;
end


% Hint: get(hObject,'Value') returns toggle state of radiobuttonPlayPureTones
tdt = handles.figureGUISetupSystem3.UserData.tdt;
objPureToneDlg = findall(0,'Tag','PlayPureTonesWindow');
% If there's already a window tagged PlayPureTonesWindow,  put it on front
if ~isempty(objPureToneDlg)
    figure(objPureToneDlg);
    return;
end

pos = get(handles.figureGUISetupSystem3, 'Position');

% Otherwise create one
tdt.checkAndSetHeadAmp;
d = dialog('Position',  [300 300 350 250],...
    'WindowStyle', 'normal', ...
    'Name', 'Playing Pure Tones',...
    'Tag', 'PlayPureTonesWindow',...
    'CreateFcn', handles.figureGUISetupSystem3.UserData.CreateFcn,...
    'CloseRequestFcn',handles.figureGUISetupSystem3.UserData.CloseRequestFcn);

set(d, 'Units', get(handles.figureGUISetupSystem3,'Units'));
set(d, 'Position', [pos(1)+pos(3)/2-77/2 pos(2)+pos(4)/2-20/2 70 20]);
set(d, 'PaperPosition',get(handles.figureGUISetupSystem3,'PaperPosition'))

txt = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(175-210/2) 200 210 40],...
    'Tag', 'txt',...
    'String', 'Play pure tones, draw receptive fields');
nbRep = uicontrol('Parent', d,...
    'Style', 'edit', ...
    'Position', [(240) 188 90 25],...
    'callback', @checkInputIsArray,...
    'Tag', 'nbRep',...
    'String', '1');
nbRepText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 185 210 25],...
    'String', 'Number of repetition per sweep:');
freq = uicontrol('Parent', d,...
    'Style', 'edit', ...
    'Position', [(240) 158 90 25],...
    'callback', @checkInputIsArray,...
    'Tag', 'freq',...
    'String', '2*10.^(1:0.25:4)');
freqText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 155 220 25],...
    'String', 'Frequencies (Matlab array) in Hz:');
loud = uicontrol('Parent', d,...
    'Style', 'edit', ...
    'Position', [(240) 128 90 25],...
    'Tag', 'loud',...
    'callback', @checkInputIsArray,...
    'String', '50:10:90');
loudText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 125 220 25],...
    'String', 'Loudness (Matlab array) in dB:');

saveMat = uicontrol('Parent', d,...
    'Style', 'checkbox', ...
    'Position', [(280) 98 20 25],...
    'Value', 1);
saveMatText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 95 220 25],...
    'String', 'Save in mat-file:');

saveFig = uicontrol('Parent', d,...
    'Style', 'checkbox', ...
    'Position', [(280) 68 20 25],...
    'Value', 1);
saveFigText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 65 220 25],...
    'String', 'Save figure:');

closeFig = uicontrol('Parent', d,...
    'Style', 'checkbox', ...
    'Position', [(280) 38 20 25],...
    'Value', 1);
closeFiggText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 35 220 25],...
    'String', 'Close figure at the end:');

btn = uicontrol('Parent', d,...
    'Position', [(175-170/2) 05 170 25],...
    'String', 'Run!', ...
    'Style', 'pushbutton',...
    'Value', 0,...
    'callback', {@runPlayTones_callback, handles,d,tdt, nbRep,freq,loud,saveMat,saveFig,closeFig});

% d on top
figure(d);

% (Wait until the ui is deleted) NO: replaced by:
% 'CloseRequestFcn', 'closereq; set(findall(0,''Tag'',''radioLazy''),''Value'',1);'
% uiwait(d);
% set(handles.radioLazy, 'Value' ,1);

function checkInputIsArray(hObject, event)
% Check what was given is a string defining a Matlab array
try
    str = evalin('base',get(hObject, 'String'));
catch ME
        set(hObject, 'String', '');
        warning(['The argument you gave (%s) is not interpreted as a '...
            'Matlab expression'], get(hObject,'String'));
        return;
end
if isempty(str)
    set(hObject, 'String', '');
    warning('The argument you gave (%s) isempty', get(hObject,'String'));
    return;
end

if ~isnumeric(str)
    try 
        evalin('base',get(hObject, 'String'))
    catch ME2
        set(hObject, 'String', '');
        warning(['The argument you gave (%s) is not interpreted as a '...
            'Matlab array or a variable in the base workspace'], get(hObject,'String'));
        return;
    end
    set(hObject, 'String', '');
    warning('The argument you gave is not interpreted as a Matlab array');
    return;
end

if strcmp(get(hObject, 'Tag'),'nbRep')
    if ~(length(str)==1 &&str==floor(str))
        error('Number of repetitions should be an integer');
    end
end

function runPlayTones_callback(hObject, event, handles,d, tdt, nbRep,freq,loud,saveMat,saveFig,closeFig)

set(hObject, 'String', 'Running...');

if isempty(tdt.figure)||~isvalid(tdt.figure),  tdt.figure = guideChannelsPlot();   end
tdt_guidata = guidata(tdt.figure);
set(tdt_guidata.radioPlay, 'Value', 1);

%guidata(hObject, handles);
drawnow;

pureTonesStr.nbRep = eval(get(nbRep, 'String')); % 1
pureTonesStr.freq  = eval(get(freq,  'String')); % 2*10.^(1:0.25:4); % Hz
pureTonesStr.loud  = eval(get(loud,  'String')); % 50:10:80; %dB SPL
pureTonesStr.saveMat = get(saveMat,  'Value');
pureTonesStr.saveFig = get(saveFig,  'Value');

[f_puretonogram, pureTones] = tdt.playPureTones(...
    'nbRep', pureTonesStr.nbRep,...
    'freq',  pureTonesStr.freq,...
    'loud',  pureTonesStr.loud,...
    'saveMat', pureTonesStr.saveMat,...
    'saveFig', pureTonesStr.saveFig); % 2*10.^(1:0.15:4),'loud',50:5:80

set(hObject, 'String', 'Run!');
set(tdt_guidata.radioStop, 'Value', 1);

if get(closeFig, 'Value')
    close(f_puretonogram);
end


function [tdt,guidata_tdt, goforit]  = prepareFigureForPlayingSounds(hObject, event, handles)
tdt = handles.figureGUISetupSystem3.UserData.tdt;

% Check headamp
tdt.checkAndSetHeadAmp;

% Hint: get(hObject,'Value') returns toggle state of radiobuttonPlaySpeechFiles
goforit = checkVariablesAreLoaded(hObject, event, handles);
if ~goforit
    set(handles.radioLazy, 'Value' ,1);
    guidata_tdt = -1;
    return;
end

if isempty(tdt.figure)||~isvalid(tdt.figure),  tdt.figure = guideChannelsPlot();   end
guidata_tdt = guidata(tdt.figure);

% guideChannelsPlot should be accessible!
% to change it: guide(fullfile(tdt.masterDirectory, 'scripts/guideChannelsPlot.fig'))
if ~exist('guideChannelsPlot.m', 'file')
    set(handles.radioLazy, 'Value' ,1);
    warning('Add path to function guideChannelsPlot');
end

% Main function: runs through the wav files, play them, plot the channels,
% keep in tdt.listWavFiles the name of the files being created in RS4

% Let's freeze the GUI
eval(handles.figureGUISetupSystem3.UserData.freezeFigureGUISetupSystem3);

% Set runningWavIndex to 0
guidata_tdt.figureguideChannelsPlot.UserData.tdt.runningWavIndex = 0;

drawnow();
set(guidata_tdt.radioPause, 'Value', 1);
set(guidata_tdt.checkFilter, 'Value', 1);
%set(guidata_tdt.checkboxFprintfInfo, 'Value', 0);
figure(tdt.figure);
continueloop = 1;
while continueloop 
    if ~isvalid(guidata_tdt.radioPause)
        continueloop = 0;
    else
        continueloop = get(guidata_tdt.radioPause, 'Value');
    end
    stoploop = get(guidata_tdt.radioStop, 'Value');
    % Prompt to confirm
    if stoploop
        stoploop = guideChannelsPlot('checkWeDoWantToStop', tdt.figure,[],guidata_tdt);
        if stoploop
            goforit = 0;
            tdt.stopTheLoop = 1;
            % Let's unfreeze the GUI
            eval(handles.figureGUISetupSystem3.UserData.unfreezeFigureGUISetupSystem3);
        else % Back to pause
            set(guidata_tdt.radioPause, 'Value', 1)
            continueloop = 1;
        end
    end
    % Delay until we click 'Play in the GUI'
    drawnow;
end


% --- Executes on button press in radiobuttonPlaySpeechFiles.
function radiobuttonPlaySpeechFiles_Callback(hObject, event, handles)
% hObject    handle to radiobuttonPlaySpeechFiles (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Choose experiment to run
cPath = fileparts(get(handles.figureGUISetupSystem3, 'FileName'));
exptFolderName = 'TDTexperimentsScripts';
dirExpt = dir(fullfile(cPath,exptFolderName));
dirExpt(1:2) = [];
nameChoices = {dirExpt.name};


d = dialog('Position',  [300 300 350 250],...
    'WindowStyle', 'normal', ...
    'Name', sprintf('Select Experiment To Run'),...
    'Tag', 'SelectExperimentToRun',...
    'CreateFcn', handles.figureGUISetupSystem3.UserData.CreateFcn,...
    'CloseRequestFcn',handles.figureGUISetupSystem3.UserData.CloseRequestFcn);

pos = get(handles.figureGUISetupSystem3, 'Position');
set(d, 'Units', get(handles.figureGUISetupSystem3,'Units'));
set(d, 'Position', [pos(1)+pos(3)/2-70/2 pos(2)+pos(4)/2-8/2 70 10.5]);
set(d, 'PaperPosition',get(handles.figureGUISetupSystem3,'PaperPosition'))

dirText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(175-370/2) 92 370 45],...
    'FontName','MS Sans Serif',...
    'FontSize', 10,...
    'String', sprintf('Select script or folder containing one \n(by default, script should have same name .m) in:'));
dirText2 = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(175-370/2) 62 370 35],...
    'FontName','MS Sans Serif',...
    'FontSize', 10,...
    'String', sprintf('%s%s',fullfile(cPath,exptFolderName), filesep));

popup = uicontrol('Parent', d,...
    'Style', 'popup', ...
    'Position', [(175-170/2) 45 170 25],...
    'FontSize', 10,...
    'String', nameChoices);
btn = uicontrol('Parent', d,...
    'Position', [(175-110/2) 12 110 25],...
    'String', 'Use this file', ...
    'Style', 'pushbutton',...
    'Value', 0,...
    'callback', {@executeExpt,handles, d,popup,fullfile(cPath,exptFolderName),nameChoices} );




function executeExpt(hObject, event,handles, d, popup,pathToExptFolder,nameChoices)
% Run the given experiment (.m script)
expt = nameChoices{get(popup,'Value')};
exptFile = fullfile( pathToExptFolder, expt );

% The default option is to have a '.m' named as the directory. Otherwise
% prompts the user to enter the script's name
if exist(exptFile, 'file')==2 % checks for files or directories
    % exptFile already the script we want to run
elseif  exist(fullfile(exptFile),              'dir') && ...
        exist(fullfile(exptFile, [expt '.m']), 'file')
    % same script name as folder
    exptFile = fullfile(exptFile, [expt '.m']);
elseif   exist(fullfile(exptFile),              'dir') && ...
        ~exist(fullfile(exptFile, [expt '.m']), 'file')
    % exptFile is only a folder
    [UIxptFile,pathname] = uigetfile(exptFile,'Choose a script to run in this folder');
    if isempty(UIxptFile) || (isnumeric(UIxptFile) && UIxptFile==0)
        return;
    end
    UIxptFile = fullfile(pathname, UIxptFile); 
    if ~strncmp(UIxptFile, exptFile, length(exptFile))
        warning('We enforce that the selected script be within %s\n',exptFile);
        beep; 
        return;
    end
    exptFile = UIxptFile;
end

% Close the dialog box
close(d);

% Reset GUI to Expt
set(handles.radiobuttonPlaySpeechFiles, 'Value', 1);

% Write name in guideSetupSystem3
set(handles.textExptName, 'String', expt)


% Makes the figure ready
[tdt,guidata_tdt, goforit] = prepareFigureForPlayingSounds(hObject, event, handles);
if ~goforit||~isvalid(guidata_tdt.radioPause)
    if tdt.stopTheLoop
        % Already unfreezed GUI within prepareFigureForPlayingSounds
        tdt.stopTheLoop = 0;
    end
    return;
end
tdt.checkAndSetHeadAmp;

% Enable comments in guideChannelsPlot
set(guidata_tdt.editComment, 'Enable', 'On');
set(guidata_tdt.editComment, 'String', '');
setappdata(0,'callbackEffectif',0)


% Check whether RS4isConnected is connected. We assume this doesn't change
% during the experiment
ticRS4 = tic;
setappdata(0, 'RS4isConnected', exist(tdt.RS4IP, 'dir')); 
if toc(ticRS4)>1
    warning('exist(tdt.RS4IP, ''dir'') Takes a long time because RS4 not connected. Consider replacing by 0 manually.');
end


%profile on
%%%%%%%%%%% THIS IS WERE THE EXPERIMENT IS USED %%%%%%%%%%%%%%%%%%%%%%%%
% Run script
run(exptFile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%profile viewer
%profile off

% unfreeze GUI
eval(handles.figureGUISetupSystem3.UserData.unfreezeFigureGUISetupSystem3);

if ~isempty(tdt.figure)&&isvalid(tdt.figure)
    % Disable comments
    set(guidata_tdt.editComment, 'Enable', 'Off');
    % Put on stop
    set(guidata_tdt.radioStop, 'Value', 1);
end

%Save
if ~isempty(tdt.ID)
tdt.saveLightTdt(fullfile(tdt.local,'tdt.mat'), ...
                'figure', 'calib_freqs', 'calib_left', 'calib_right', ...
                'RZ', 'RX', 'zBus', 'PA_left', 'PA_right')
end

set(handles.radioLazy, 'Value' ,1);

% Generate a report if asked to
if get(handles.checkboxGenerate, 'Value')
    pushbuttonGenerate_Callback(hObject, event, handles)
end


% --- Executes on button press in radioNoise.
function radioNoise_Callback(hObject, event, handles)
% hObject    handle to radioNoise (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


dNoise = findall(0,'Tag', 'ExploringNoiseBursts');
if ~isempty(dNoise)
    figure(dNoise);
    return;
end

% Hint: get(hObject,'Value') returns toggle state of radioNoise
goforit = checkVariablesAreLoaded(hObject, event, handles);
if ~goforit
    set(handles.radioLazy, 'Value' ,1);
    return;
end

tdt = handles.figureGUISetupSystem3.UserData.tdt;
tdt.checkAndSetHeadAmp;

dd = findall(0,'Tag', 'figureGUISetupSystem3');
pos = get(dd, 'Position');


d = dialog('Position', [300 300 350 250],...
    'WindowStyle', 'normal',...
    'Name', 'Exploring with noise bursts', ...
    'Tag',  'ExploringNoiseBursts',...
    'CreateFcn', handles.figureGUISetupSystem3.UserData.CreateFcn,...
    'CloseRequestFcn',handles.figureGUISetupSystem3.UserData.CloseRequestFcn);

set(d, 'Units', get(dd,'Units'));
set(d, 'Position', [pos(1)+pos(3)/2-77/2 pos(2)+pos(4)/2-20/2 67 24]); % 19
set(d, 'PaperPosition',get(dd,'PaperPosition'))

txt = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(175-210/2) 250 210 40],... 200 210 40
    'Tag', 'txt',...
    'String', 'Bursts of white noise played during electrode descent');

freqs = uicontrol('Parent', d,...
    'Style', 'edit', ...
    'Position', [(240) 218 90 25],... 25
    'callback', @checkInputIsArray,...
    'Tag', 'freqOrNaN',...
    'String', 'NaN');
freqTxt = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 215 210 25],...
    'String', 'Frequencies (NaN for white noise):');


timeBurst = uicontrol('Parent', d,...
    'Style', 'edit', ...
    'Position', [(240) 188 90 25],... 25
    'callback', @checkInputIsArray,...
    'Tag', 'timeBurst',...
    'String', '100');
timeBurstText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 185 210 25],...
    'String', 'Duration of bursts (in ms):');

timeSil = uicontrol('Parent', d,...
    'Style', 'edit', ...
    'Position', [(240) 158 90 25],...
    'callback', @checkInputIsArray,...
    'Tag', 'timeSil',...
    'String', '200');
timeSilText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 155 220 25],...
    'String', 'Silences between bursts (in ms):');
loud = uicontrol('Parent', d,...
    'Style', 'edit', ...
    'Position', [(240) 128 90 25],...
    'Tag', 'loudNoiseBurst',...
    'callback', @checkInputIsArray,...
    'String', '50:10:70');
loudText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 125 220 25],...
    'String', 'Loudness (in dB):');

% Number of bursts played/to play
numberBurstToPlay = uicontrol('Parent', d,...
    'Style', 'edit', ...
    'Position', [(240) 98 90 25],...
    'Tag', 'numberBurstToPlay',...
    'callback', @checkInputIsArray,...
    'String', '');
numberBurstToPlayText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-210/2) 95 220 25],...
    'String', 'Number of bursts to play (empty or integer):');
numberBurstPlayed = uicontrol('Parent', d,...
    'Style', 'text', ...
    'Position', [(240) 64 90 25],...
    'Tag', 'numberBurstPlayed',...
    'callback', @checkInputIsArray,...
    'String', '0');
numberBurstPlayedText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 65 220 25],...
    'String', 'Number of bursts played:');

% Delete all files that were created during the exploration
deleteRecords = uicontrol('Parent', d,...
    'Style', 'checkbox', ... %pushbutton
    'Position', [277.5 38 170 25],... %[(260) 98 50 25],...
    'Value', 1); %'String', 'del');
deleteRecordsText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(105-230/2) 35 220 25],...
    'String', 'Delete files saved during exploration');

btn = uicontrol('Parent', d,...
    'Position', [(175-170/2) 05 170 25],...
    'String', 'Open figure to play', ...
    'Style', 'pushbutton',...
    'Value', 0,...
    'callback', {@runNoisBursts_callback, handles,d,tdt, freqs,timeBurst,timeSil,loud,deleteRecords});
% Wait until the ui  is deleted (replaced by:
% 'CloseRequestFcn', 'closereq; set(findall(0,''Tag'',''radioLazy''),''Value'',1);'
%uiwait(d);
% set(handles.radioLazy, 'Value' ,1);


function runNoisBursts_callback(hObject, event, handles,d,~,freqs, timeBurst,timeSil,loud,deleteRecords)
tdt = prepareFigureForPlayingSounds(hObject, event, handles);

set(hObject, 'String', 'Running...');
% 
pureTonesStr.type = eval(get(freqs, 'String')); % 'whiteNoise';
pureTonesStr.burstLength = eval(get(timeBurst, 'String')); % 1
pureTonesStr.silenceLength  = eval(get(timeSil,  'String')); % 2*10.^(1:0.25:4); % Hz
pureTonesStr.loud  = eval(get(loud,  'String')); % 50:10:80; %dB SPL
% Keep the choice of deletion, and make the icon unclickable

disp('runNoisBursts_callback: deleteRecords to do'); 
timeBeforeExploration = now;

tdt.playShortBursts(...
    'burstLength', pureTonesStr.burstLength,...
    'silenceLength',  pureTonesStr.silenceLength,...
    'loud',  pureTonesStr.loud,...
    'type',  pureTonesStr.type); % 2*10.^(1:0.15:4),'loud',50:5:80

% dirRS4 = dir(tdt.RS4IP);
% dirRS4 = dirRS4(3:end); % get rid of . and ..
% a = [dirRS4.datenum] > timeBeforeExploration;
% SHOULD BE FINE NOW('Solve naming problem (datenum)');

set(hObject, 'String', 'Run!');

close(d);


% --- Executes on button press in pushbuttonAssigntdt.
function pushbuttonAssigntdt_Callback(hObject, event, handles)
% hObject    handle to pushbuttonAssigntdt (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base', 'tdt', handles.figureGUISetupSystem3.UserData.tdt)
fprintf('Assigning ''tdt'' within base workspace.\n')


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over textTitleSetup.
function textTitleSetup_ButtonDownFcn(hObject, event, handles)
% hObject    handle to textTitleSetup (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Yo')



% --- Executes on button press in pushbuttonSetHeadamp.
function pushbuttonSetHeadamp_Callback(hObject, event, handles)
% hObject    handle to pushbuttonSetHeadamp (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
headampexpect = handles.figureGUISetupSystem3.UserData.tdt.headampAttenExpected;
fprintf('In the window ''Headamp setting'', put ''%d'' in both, click ''Change'' and close.', headampexpect);
fprintf('This won''t be checked. You are the responsible adult here.\n' );
system('Path/to/set_headamp.exe');
handles.figureGUISetupSystem3.UserData.tdt.headampAtten = headampexpect;


% --- Executes on button press in pushbuttonReloadCircuits.
function pushbuttonReloadCircuits_Callback(hObject, event, handles)
% hObject    handle to pushbuttonReloadCircuits (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tdt = handles.figureGUISetupSystem3.UserData.tdt;
fprintf(1, ' - Reloading the circuits to RZ2 and RX8:\n');
tdt.connect('RX8', 'RZ2'); %, 'ZBUS', 'PA_left', 'PA_right');


% --- Executes on button press in pushbuttonSetChannels.
function pushbuttonSetChannels_Callback(hObject, event, handles)
% hObject    handle to pushbuttonSetChannels (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


in = inputdlg('Enter a Matlab Array for the channels you want to record from',...
    'Channels to record from (others are deleted)', 1, {'1:128'});
if isempty(in)
    return
end
try 
    chan = eval(in{1});
catch ME
    disp(ME);
    return;
end
if ~isnumeric(chan) || any(chan-floor(chan))|| max(chan)>128 || min(chan)<1
    fprintf('Enter a proper array of integers between 1 and 128.\n');
    return;
end    
fprintf('Saving channels (others will be deleted after being saved on the RS4):%s.\n',sprintf(' %d', chan(:)));
handles.figureGUISetupSystem3.UserData.tdt.arrayChanRec = chan;
set(handles.textChannelsRecordFrom, 'String', ['Channels to record from: ' in{1}]);


% --- Executes on button press in checkboxGenerate.
function checkboxGenerate_Callback(hObject, event, handles)
% hObject    handle to checkboxGenerate (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxGenerate


% --- Executes on button press in pushbuttonGenerate.
function pushbuttonGenerate_Callback(hObject, event, handles)
% hObject    handle to pushbuttonGenerate (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Boolean to open the report after writing
showReport = 1;

tdt = handles.figureGUISetupSystem3.UserData.tdt;
if ~exist(tdt.local, 'dir')
    beep;
    warning('Cannot generate report if tdt.local hasn''t been set.\n');
    return;
end

if 1==0
    % Will generate a pdf summarising the experiment: calibrations, receptive
    % fields, all comments written, some statistics
    beep;
    warning('guideSetupSystem3: pushbuttonGenerate_Callback: Not coded yet (licensing problem)');
    return;
end

% No mlreportgen.ppt found...
%import mlreportgen.ppt.Presentation
%slidesFile = fullfile(tdt.local, sprintf('summaryExpt_%s.pptx', tdt.ID));
%slides = Presentation(slidesFile);
%winopen(slides);

% in .docx for now because the .html of .pdf were not working
import mlreportgen.dom.*;
%reportFile = fullfile(tdt.local, sprintf('reportExpt_%s', tdt.ID));
reportFile = fullfile(tdt.local, sprintf('reportExpt_%s.docx', tdt.ID));
rep = Document(reportFile, 'docx'); % breaks if document already open

% o = OrderedList({'a' 'bb' 'ccc'});
% append(rep, o);

% TITLE
a = {''};
a{end+1} = '------------------   REPORT  ------------------';
a{end+1} = sprintf('Report generated on %s', datestr(now));



% DESCRIPTION OF THE SETUP
a{end+1} = '------------------   ELECTRODES  ------------------';
a{end+1} = sprintf('Total number of channels used: %d', tdt.nb_channels);
% A sentence per channel
if isfield(handles.figureGUISetupSystem3.UserData,'electrodesConfigs')
    electrodesConfigList = handles.figureGUISetupSystem3.UserData.electrodesConfigList;
    pathToAllConfig = cellfun(@(file)fullfile(tdt.masterDirectory, 'electrodeConfigurations', file),...
        electrodesConfigList, 'unif', false);
    electrodesConfigs = handles.figureGUISetupSystem3.UserData.electrodesConfigs;
    numChannelPerElec = cellfun(@(conf)numel(dlmread(conf, ' ')),pathToAllConfig);
    for kk = 1:length(electrodesConfigs)
        a{end+1} = sprintf(' + Electrode %d (%d channels) is configured with config %s',kk,numChannelPerElec(electrodesConfigs(kk)),...
            electrodesConfigList{electrodesConfigs(kk)} );
    end
end

cellfun(@(d)append(rep,d),a, 'Unif', false);

a = {};
a{end+1} = '------------------   RECEPTIVE FIELDS  ------------------';
rf = dir(fullfile(tdt.local,'receptiveFields*.jpg'));
for kk=1:length(rf)
    cFile = fullfile(tdt.local,rf(kk).name);
    a{end+1} = cFile;
    % load image and set width at 700 px
    img = loadAndResizeImg(hObject, event, cFile, 700);
    append(rep, img);
end
cellfun(@(d)append(rep,d),a, 'Unif', false);


% CALIBRATION FIGURES
append(rep, '------------------   CALIBRATION  ------------------');
sides = {'left', 'right'};
for kk=1:2
    side = sides{kk};
    cFile = fullfile(tdt.local,['calib_' side '_Splines.png']);
    % Trying different combinations. Might break, in which case the used
    % has to change the code accordingly
    if ~exist(cFile, 'file')
         cFile = fullfile(tdt.local,['calib_' side '_Polynomials.png']);
         if ~exist(cFile, 'file')
             warning(['The %s calibration picture %s wasn''t found, with '...
                 'Splines or Polynomials. Might be a problem in new settings'],...
                 side, cFile);
             append(rep, sprintf('Plot of calibration of %s ear wasn''t found.', side));
             beep;
             continue;
         end
    end
    img = loadAndResizeImg(hObject, event, cFile, 700);
    %append(rep, sprintf('Calibration of %s ear, with fit.', side));
    append(rep, img);
end

%calib = {''};
%cellfun(@(d)append(rep,d),calib, 'Unif', false);


% COMMENTS MADE DURING THE EXPERIMENT
append(rep, '------------------   COMMENTS  ------------------');

%tdt_guidata = guidata(tdt.figure);
%CommentsMade = get(tdt_guidata.CommentsMade, 'String');
%CommentsMade = (tdt_guidata.figureguideChannelsPlot.UserData.commentInMemory);
if isprop(tdt, 'listWavFiles')
    if isfield(tdt.listWavFiles, 'comment')
        CommentsMade = {tdt.listWavFiles.comment};
        comments = {''...
            'Comments made during the experiment'...
            CommentsMade(cellfun(@(s)~isempty(s),CommentsMade))};
    else
        comments = {'No field ''comment'' was found in tdt.listWavFiles'};
    end
    cellfun(@(d)append(rep,d),comments, 'Unif', false);
else
    warning('tdt.listWavFiles is expected to be a structure and ''.comment'' contain the comments made during the experiment');
    beep
end


append(rep, '------------------   TDT OBJECT  ------------------');
append(rep, 'Not implementerd yet');

% CLOSING
cellfun(@(d)append(rep,d),{'' '' '' '' '...The Lannisters send their regards'}, 'Unif', false);
close(rep);




% open report
if showReport
    fprintf(' --- open(%s)\n',reportFile);
    %type(reportFile)
    open(reportFile)
end

 

function img = loadAndResizeImg(hObject, event, imgFile, nbPx)
import mlreportgen.dom.*;
% Used to create a DOM image object at nbPx pixels to append
img = Image(imgFile);
wdh = textscan(img.Width,  '%dpx');
hgt = textscan(img.Height, '%dpx');
wdh = double(wdh{1});
hgt = double(hgt{1});
widthChosen = nbPx;%px
ratio = widthChosen/wdh;
img.Width  = sprintf('%dpx', widthChosen); % set image's width
img.Height = sprintf('%dpx', floor(hgt*ratio)); % set image's height to keep ratio



% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figureGUISetupSystem3_WindowButtonDownFcn(hObject, event, handles)
% hObject    handle to figureGUISetupSystem3 (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Get the current value
cButton = get(get(handles.uiUseSystem3, 'SelectedObject'), 'Tag');
% cPoint =  get(handles.figureGUISetupSystem3, 'CurrentPoint');
% cPos =    get(handles.figureGUISetupSystem3, 'Position');
% Return if radioLazy or not recognised, otherwise beep and other window on top
switch cButton
    case 'radioLazy',                   
        return;  % Legit
    case 'radiobuttonRunCalib',         cFig = findall(0,'Tag','guideRunCalibrationDlg');
    case 'radioNoise',                  cFig = findall(0,'Tag','ExploringNoiseBursts');
    case 'radiobuttonPlayPureTones',    cFig = findall(0,'Tag','PlayPureTonesWindow');
    case 'radiobuttonPlaySpeechFiles',  % SelectExperimentToRun has the priority
         cFig = findall(0,'Tag','SelectExperimentToRun');
         if isempty(cFig)
             cFig = findall(0,'Tag','figureguideChannelsPlot');
         end
    otherwise, beep; ...
            fprintf(1, 'figureGUISetupSystem3_WindowButtonDownFcn: Tag %s not recognised.\n', cButton); return;
end

    
% Beep and on top
beep;
figure(cFig);


% --- Executes on mouse motion over figure - except title and menu.
function figureGUISetupSystem3_WindowButtonMotionFcn(hObject, event, handles)
% hObject    handle to figureGUISetupSystem3 (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%figureGUISetupSystem3_WindowButtonDownFcn(hObject, event, handles)


% --- Executes on button press in pushbuttonDummyToLoad.
function pushbuttonDummyToLoad_Callback(hObject, event, handles)
% hObject    handle to pushbuttonDummyToLoad (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tdt = handles.figureGUISetupSystem3.UserData.tdt;

% Dummy button to save time while development by doing everything at once
tdt.calib_left_file  = 'Path/to/calib_left_freq_resp.mat';
tdt.calib_right_file = 'Path/to/calib_left_freq_resp.mat'; 
tdt.depth = 0.00;
tdt.RS4IP = '\\128.243.7.88\data';
tdt.ID = '12_Jul_2016_14_05_40'; 
tdt.local = 'Path/to/12_Jul_2016_14_05_40'; 

handles.figureGUISetupSystem3.UserData.tdt.load_calibration_single(tdt.calib_left_file, 'left');
handles.figureGUISetupSystem3.UserData.tdt.load_calibration_single(tdt.calib_right_file, 'right');


tdt.headampAtten = tdt.headampAttenExpected;
warning('Check headamp at -20. Set it if not!');
fprintf('Set headamp with GUI if not done yet!\n');
set(handles.checkboxCalibLeft, 'Value',1)
set(handles.checkboxCalibRight, 'Value',1)
set(handles.checkboxRS4ID, 'Value',1)
set(handles.editDepth, 'String', '0');
set(handles.checkboxDepth, 'Value',1)
set(handles.checkboxMakeID, 'Value',1)
fprintf(1,'pushbuttonDummyToLoad_Callback: made invisible\n');
set(findall(0,'Tag','pushbuttonDummyToLoad'),'Visible','off')


% --- Executes on button press in pushbuttonSetElectrodeCongig.
function pushbuttonSetElectrodeCongig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSetElectrodeCongig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set the electrode's configuration when plotting the receptive fields

d = dialog('Position',  [300 300 350 250],...
    'WindowStyle', 'normal', ...
    'Name', sprintf('Select Electrodes Configurations'),...
    'Tag', 'SelectElectrodeConfigurations',...
    'CreateFcn',      handles.figureGUISetupSystem3.UserData.CreateFcn,...
    'CloseRequestFcn',handles.figureGUISetupSystem3.UserData.CloseRequestFcn);

pos = get(handles.figureGUISetupSystem3, 'Position');
set(d, 'Units', get(handles.figureGUISetupSystem3,'Units'));
set(d, 'Position', [pos(1)+pos(3)/2-70/2 pos(2)+pos(4)/2-8/2 70 10.5]);
set(d, 'PaperPosition',get(handles.figureGUISetupSystem3,'PaperPosition'))

dirText = uicontrol('Parent', d, ...
    'Style', 'text', ...
    'Position', [(175-370/2-60) 27 370 25],...
    'FontName','MS Sans Serif',...
    'String', sprintf('Select Number of electrodes:'),...
    'Tag', 'numElecText');

% Default number of electrodes
numElec = 2;
if isfield(handles.figureGUISetupSystem3.UserData,'electrodesConfigs')
    numElec = length(handles.figureGUISetupSystem3.UserData.electrodesConfigs);
end

popup = uicontrol('Parent', d,...
    'Style', 'edit', ...
    'Position', [(175-40/2+60) 27 40 25],...
    'String', num2str(numElec), ...
    'Tag', 'selectNumberElectrodes');
btn = uicontrol('Parent', d,...
    'Position', [(175-110/2) 72 110 25],...
    'String', 'Select configurations:', ...
    'Tag', 'SelectConfigurationsPush',...
    'Style', 'pushbutton',...
    'Value', 0,...
    'callback', {@Selectconfigurationforeachelectrode,handles, d,popup} );

function Selectconfigurationforeachelectrode(hObject, eventdata, handles, varargin)
% Open window to select a configuration file for each electrode
% varargin is {d, popup} if we ask through the GUI. If no other argument,
% the function simply iniitialise everything at the default options
listOfElectrodeConfigurations = dir(fullfile(handles.figureGUISetupSystem3.UserData.tdt.masterDirectory, 'electrodeConfigurations'));
listOfElectrodeConfigurations = {listOfElectrodeConfigurations(3:end).name};

% To display with the disp button of GUI, and  the pure tone button
handles.figureGUISetupSystem3.UserData.electrodesConfigList = listOfElectrodeConfigurations;
handles.figureGUISetupSystem3.UserData.electrodesConfigFolder = fullfile(handles.figureGUISetupSystem3.UserData.tdt.masterDirectory, 'electrodeConfigurations');


if ~isempty(varargin)
    d = varargin{1};
    popup = varargin{2};
else 
    % By default, 2 neuroNexusA8x8.txt
    conf = find(strcmp(listOfElectrodeConfigurations, 'neuroNexusA8x8.txt'));
    handles.figureGUISetupSystem3.UserData.electrodesConfigs = [conf; conf];
    return;
end


% check what's given is an integer
numElec = str2double(get(popup, 'String'));
if isnan(numElec) || numElec<1 || numElec-floor(numElec)~=0 || numElec>16
    beep;
    if numElec>16
        fprintf('Selectconfigurationforeachelectrode: An arbitrary limit of 16 electrodes has been set.\n')
    end
    set(popup, 'String', '2');
    return
end

% Give a default setting if already set
if isfield(handles.figureGUISetupSystem3.UserData,'electrodesConfigs')
    popupElecDef = handles.figureGUISetupSystem3.UserData.electrodesConfigs;
    if length(popupElecDef)>numElec
        popupElecDef = popupElecDef(1:numElec);
    end
    if length(popupElecDef)<numElec
        popupElecDef = [popupElecDef; ones(numElec-length(popupElecDef),1)];
    end    % If we kept the strings
    %popupElecDef = struct2cell(handles.figureGUISetupSystem3.UserData.electrodesConfigs);
    %popupElecDef = cellfun(@(x)x{1}, popupElecDef, 'Unif', false); % first elems
else 
    popupElecDef = ones(1,numElec); %arrayfun(@(x)listOfElectrodeConfigurations{1}, 1:numElec, 'Unif', false);
end


set(popup, 'Enable' , 'off');
set(findall(d, 'Tag','SelectConfigurationsPush'), 'Visible', 'off');
numElecText = findall(d, 'Tag', 'numElecText');
set(numElecText, 'Visible', 'off');
% Push button t odisplay the configurations on command window
btn = uicontrol('Parent', d,...
    'Position', [(175-100/2-100) 27 120 25],...
    'String', 'Display Configurations', ...
    'Style', 'pushbutton',...
    'Value', 0,...
    'callback', {@displayAllConfigurations,handles,d,listOfElectrodeConfigurations} );

if numElec>2
    set(d,'Position',get(d,'Position').*[1 1 1 (1+0.22*(numElec-2))])
end
for kk = 1:numElec
    dirText.(sprintf('elec%d',kk)) = uicontrol('Parent', d, ...
        'Style', 'text', ...
        'Position', [(175-370/2-90) 68+30*(numElec-kk) 370 25],...
        'FontName','MS Sans Serif',...
        'String', sprintf('Select config electrode %d', kk));
    % In case we delete some configs after setting some, it would break
    if popupElecDef(kk)>length(listOfElectrodeConfigurations)
        popupElecDef(kk) = 1;
    end
    popupElec.(sprintf('elec%d',kk)) = uicontrol('Parent', d,...
        'Style', 'popup', ...
        'Position', [(175-120/2+60) 73+30*(numElec-kk) 170 25],...
        'String', listOfElectrodeConfigurations, ...
        'Value', popupElecDef(kk));
    
end
btn = uicontrol('Parent', d,...
    'Position', [(175-50/2+120) 28 55 25],...
    'String', 'Set them!', ...
    'Style', 'pushbutton',...
    'Value', 0,...
    'callback', {@setElectrodesGivenPopups,handles,d,popupElec} );

function setElectrodesGivenPopups(hObject, eventdata, handles, d, popupElec)
% Set the electrodes configs in UserData and close d

% To keep the strings:
%    structfun( @(pop)get(pop, 'String'),popupElec, 'UniformOutput', false);
elecVal = structfun( @(pop)get(pop, 'Value'),popupElec);
[goodToGo, totalNbChannels] = checkConfigsAreOK(hObject, eventdata, handles, d, popupElec,elecVal);
if ~goodToGo
    beep;
    %warning('Problem with some electrode configuration');
    return
end

% If everything is fine, keep those in memory
handles.figureGUISetupSystem3.UserData.electrodesConfigs = elecVal;

% Keep in tdt.nb_channels the number of channels this represents
handles.figureGUISetupSystem3.UserData.tdt.nb_channels = totalNbChannels;

close(d);


function [isgud,totalNbChannels] = checkConfigsAreOK(hObject, eventdata, handles, d, popupElec,elecVal)
% Routine to check a channel file is ok (all values appear one, no jump)
isgud = 1;
totalNbChannels = -1;

configs = elecVal; %handles.figureGUISetupSystem3.UserData.electrodesConfigs;
cInd = 1;
Xpos = ones(128,1);
Ypos = ones(128,1);
numElec = length(configs);
numChannelsPerElec = zeros(numElec,1);
for confInd = 1:numElec
    conf = handles.figureGUISetupSystem3.UserData.electrodesConfigList{configs(confInd)};
    cFile = fullfile(handles.figureGUISetupSystem3.UserData.electrodesConfigFolder, ...
        conf);
    try 
        mat = dlmread(cFile, ' ');
        numChannelsPerElec(confInd) = numel(mat);
    % Find the positions
    for jj=1:numel(mat)
        [cX, cY] = find(mat==jj);
        Xpos(cInd) = cX + (confInd-1)/numElec;
        Ypos(cInd) = cY;
        cInd = cInd+1;
    end
    catch ME
        fprintf('Configuration electrodeConfigurations %s has a problem, ', conf);
        switch ME.identifier
            case 'MATLAB:index_assign_element_count_mismatch', 
                str = sprintf('the value %d appears %d times\n',jj, max(length(cX), length(cY)));
            otherwise, 
                str = sprintf('solve it\n');
        end
        fprintf(str);
        isgud = 0;
    end
    mat = mat(:); % simpler
    if ~all(arrayfun(@(jj)any(mat==jj), 1:numel(mat)));
        fprintf('Some values appear to be missing in the conf %s, check it\n',conf);
        isgud = 0;
    end
end
if isgud
    totalNbChannels = sum(numChannelsPerElec);
end


function displayAllConfigurations(hObject, eventdata,handles,d,listOfElectrodeConfigurations)
fprintf('Displaying all electrode configurations in folder %s (Heads are down):','electrodeConfigurations');
for kk=listOfElectrodeConfigurations
    conf = kk{1};
    fprintf('\n');
    fprintf(' - Confguration %s:\n', conf);
    type(fullfile(handles.figureGUISetupSystem3.UserData.tdt.masterDirectory, 'electrodeConfigurations',conf))
    fprintf('\n');
end


% --- Executes on button press in textTitleSetup.
function textTitleSetup_Callback(hObject, event, handles)
% hObject    handle to textTitleSetup (see GCBO)
% event  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = {}; 
str{end+1} = '   /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\                               ';
str{end+1} = '                                                                                      ';
str{end+1} = '    Welcome to the IHR-made GUI to run the TDT system.                                ';
str{end+1} = '  GUI Author: Alban, May 2016.                                                        ';
str{end+1} = '  By default, the GUI is in lazy mode. The 3 ''Play'' actions                         ';
str{end+1} = '  are using a common figure called ''guideChannelsPlot'', explained lower.            ';
str{end+1} = '                                                                                      ';
str{end+1} = '    You have 4 actions available (Panel ''Use System 3''),                            ';
str{end+1} = '  which require some variables (panel ''Variables to load manually''),                ';
str{end+1} = '  and some other actions (Panel ''Other Actions'').                                   ';
str{end+1} = '                                                                                      ';
str{end+1} = '          PANEL Use System 3                                                          ';
str{end+1} = '                                                                                      ';
str{end+1} = ' - To calibrate the system, click on ''Run Calibratios''.                             ';
str{end+1} = '   You''ll be asked to select a folder and name                                       ';
str{end+1} = ' - To play some white noises during the brain exploration,                            ';
str{end+1} = '   click ''Play Noise''. A window will ask for                                        ';
str{end+1} = '   loudness, duration of sounds and silences. The silences will be                    ';
str{end+1} = '   actually longer than what you enter, due to computations within the system.        ';
str{end+1} = '   You can dynamically change those values. You can enter Matlab arrays,              ';
str{end+1} = '   in which case the value used will be randomly chosen among the ones given.         ';
str{end+1} = ' - To obtain receptive fields, click ''Play Pure Tones''. A window will               ';
str{end+1} = '   ask you which frequencies and sound level you want to play                         ';
str{end+1} = '   A big figure will show in real time the activity of the 128 channels for           ';
str{end+1} = '   each selected freqeuncy and level. If you have more than one repetition for each,  '; 
str{end+1} = '   they are summed over. Options allow you to save the figure (.fig, .pdf, .jpg) and the data (.mat))';
str{end+1} = ' - To play the wav files (TIMIT database), click ''Play Speech Files''.               ';
str{end+1} = '   Press ''Play'' in guideChannelsPlot to start playing all files                     ';
str{end+1} = '   (premade structure, see function timit2cell.m for details)                         ';
str{end+1} = '                                                                                      ';
str{end+1} = '          PANEL Variables to load manually                                            ';
str{end+1} = '                                                                                      ';
str{end+1} = ' These variables are necessary to guaranty a weel-behaved GUI.                        ';
str{end+1} = ' They are savedin the tdt object kept within the figure,                              ';
str{end+1} = ' at the properties ''calib_left_file'', ''calib_right_file'',                         ';
str{end+1} = ' ''depth'', , ''RS4IP'', ,''ID'' respectively.                                        ';
str{end+1} = ' - Left and Right Calibration: calibration files made using ''Run Calibrations''      ';
str{end+1} = '   and saved somewhere. A dialog box hhelps you load the calibration you              ';
str{end+1} = ' - RS4 IP: A dialogue box asks you for the RS4''s IP adress,                          ';
str{end+1} = '   which can be found using the touch-screen of the device, panel ''Status''.         ';
str{end+1} = '   want to use (done at the very beginning of each experiment). Required for ''Play'' actions.';
str{end+1} = ' - Electrode Depth: Enter the depth (double) at which the electrode is.               ';
str{end+1} = '   Required for ''Play'' actions''.                                                   ';
str{end+1} = ' - Make ID: Each experiment will create many additional files.                        ';
str{end+1} = '   This identifier (a folder name) will be used to keep them together,                ';
str{end+1} = '   except for the data saved in RS4. By default, ID is the current day and time.      ';
str{end+1} = ' - Display Loaded Variables: Print the value of those variables in the Matlab window. ';
str{end+1} = '                                                                                      ';
str{end+1} = '          PANEL Other Actions                                                         ';
str{end+1} = '                                                                                      ';
str{end+1} = ' - Assign tdt in base workspace: Loads the tdt variable, object of class myTDT        ';
str{end+1} = '   used by the GUIs, to make some changes manually. Since this is a handle,           ';
str{end+1} = '   changes will be applied to the variable used by the GUI.                           ';
str{end+1} = '   Run ''doc myTDT'' for more details about the myTDT class                           ';
str{end+1} = ' - Manually set headamp values: Launches an executable that will                      ';
str{end+1} = '   set the values of the headamp attenuator (top box, above the TDT hardware).        ';
str{end+1} = '   YOU ARE EXPECTED TO SET BOTH TO 20, THIS IS YOUR RESPONSABILITY                    ';
str{end+1} = '   After this, you should see the lights on the box selecting -20 (horizontally and vertically).';
str{end+1} = '                                                                                      ';
str{end+1} = '          FIGURE guideChannelsPlot                                                    ';
str{end+1} = '                                                                                      ';
str{end+1} = '  This figure contains 2 axes; the left one will show what was previously recorded    ';
str{end+1} = '   on all channels given in the top right box. Give it a Matlab array of integers between 1 and 128.';
str{end+1} = '   The right plot show the mean activity calculated for all channels between two moments:';
str{end+1} = '   the minimum time and maximum time are given, in miliseconds, in the ''Measure Activity'' panel.';
str{end+1} = '  The activity being shown is simply the L1 norm of those values (sum(abs(channel))).';
str{end+1} = '  A high-pass filter (defined within the tdt object) can be applied to the signal, always renormalised.';
str{end+1} = '  Some time information can be printed onto the command window by checking the ''Fprintf'' box.';
str{end+1} = '  Press ''Play'' to start playing. If you stopped, you''ll need to reclick on the SetupSystem3 Action panel.';
str{end+1} = '                                                         ';
str{end+1} = '          FAQ                                            ';
str{end+1} = '                                                         ';
str{end+1} = ' HOW CAN I CHANGE THE GUI?';
str{end+1} = ' To change the GUI, depending on what you wish to change, it''ll be one of the two:';
str{end+1} = '    - open(''path\to\the\GUI\guideSetupSystem3.m'');        ';
str{end+1} = '    - guide(''path\to\the\GUI\guideSetupSystem3.fig'');     ';
str{end+1} = '                                                          ';
str{end+1} = ' THIS GUI IS AMAZING, HOW CAN I CONTACT THE DEVELOPPER?';
str{end+1} = ' First, ask Chris Sumner                                    ';
str{end+1} = '                                                             ';
str{end+1} = '                                                             ';
str{end+1} = '                                                             ';
str{end+1} = '  _______  _______  _______     _________          _______    ';
str{end+1} = ' (       )(  ____ )(  ____ \    \__   __/|\     /|(  ____ )   ';
str{end+1} = ' | () () || (    )|| (    \/       ) (   | )   ( || (    )|   ';
str{end+1} = ' | || || || (____)|| |             | |   | (___) || (____)|   ';
str{end+1} = ' | |(_)| ||     __)| |             | |   |  ___  ||     __)   ';
str{end+1} = ' | |   | || (\ (   | |             | |   | (   ) || (\ (      ';
str{end+1} = ' | )   ( || ) \ \__| (____/\    ___) (___| )   ( || ) \ \__   ';
str{end+1} = ' |/     \||/   \__/(_______/    \_______/|/     \||/   \__/   ';
str{end+1} = '                                                            ';
str{end+1} = '                                                            ';
str{end+1} = '   \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/     ';
fprintf('%s\n',str{:});


% --- Executes on key press with focus on textTitleSetup and none of its controls.
function textTitleSetup_KeyPressFcn(hObject, event, handles)
% hObject    handle to textTitleSetup (see GCBO)
% event  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

