function varargout = guideChannelsPlot(varargin)
%GUIDECHANNELSPLOT M-file for guideChannelsPlot.fig
%      GUIDECHANNELSPLOT, by itself, creates a new GUIDECHANNELSPLOT or raises the existing
%      singleton*.
%
%      H = GUIDECHANNELSPLOT returns the handle to a new GUIDECHANNELSPLOT or the handle to
%      the existing singleton*.
%
%      GUIDECHANNELSPLOT('Property','Value',...) creates a new GUIDECHANNELSPLOT using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to guideChannelsPlot_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUIDECHANNELSPLOT('CALLBACK') and GUIDECHANNELSPLOT('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUIDECHANNELSPLOT.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guideChannelsPlot

% Last Modified by GUIDE v2.5 12-Jul-2016 19:09:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guideChannelsPlot_OpeningFcn, ...
                   'gui_OutputFcn',  @guideChannelsPlot_OutputFcn, ...
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


% --- Executes just before guideChannelsPlot is made visible.
function guideChannelsPlot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for guideChannelsPlot
handles.output = hObject;

set(handles.axes1,'XTick',[]);
set(handles.axes1,'YTick',[]);

handles.figureguideChannelsPlot.UserData.commentInMemory = '';

%set(hObject, 'UserData', struct('tdt',0, 'recChannel', 0,'sweep', 0));
tdt = AlbanTDT();
hObject.UserData.tdt = tdt; % initialize it for simplicity
hObject.UserData.recChannel = 0;
hObject.UserData.sweep = tdt.sweep;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guideChannelsPlot wait for user response (see UIRESUME)
% uiwait(handles.figureguideChannelsPlot);


% --- Outputs from this function are returned to the command line.
function varargout = guideChannelsPlot_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Check whether there is something to plot
if isequal(handles.figureguideChannelsPlot.UserData.recChannel, 0) && isempty(varargin)
    return;
end

% We assume those 3 variables are either given, or in the base workspace
varBaseWorksp.tdt        = evalin('base', 'exist(''tdt'',''var'')');
varBaseWorksp.recChannel = evalin('base', 'exist(''recChannel'',''var'')');
varBaseWorksp.sweep      = evalin('base', 'exist(''sweep'',''var'')');

if ~isempty(varargin),                          tdt = varargin{1}; %myTDT object with all info needed
elseif ~isempty(handles.figureguideChannelsPlot.UserData.tdt),  tdt = handles.figureguideChannelsPlot.UserData.tdt;
elseif varBaseWorksp.tdt,                       tdt = evalin('base', 'tdt');
else                                            error('tdt required');
end
if length(varargin) >= 2,                       recChannel = varargin{2};
elseif ~isempty(handles.figureguideChannelsPlot.UserData.recChannel),recChannel = handles.figureguideChannelsPlot.UserData.recChannel;
elseif varBaseWorksp.recChannel,                recChannel = evalin('base', 'recChannel');
else                                            error('recChannel required');
end
if length(varargin) >= 3,                       sweep = varargin{3};
elseif ~isempty(handles.figureguideChannelsPlot.UserData.sweep),sweep = handles.figureguideChannelsPlot.UserData.sweep;
elseif varBaseWorksp.sweep,                     sweep = evalin('base', 'sweep');
else                                            error('sweep required');
end




% Set user data (useful for pause or stop)
handles.figureguideChannelsPlot.UserData.tdt = tdt;
handles.figureguideChannelsPlot.UserData.recChannel = recChannel;
handles.figureguideChannelsPlot.UserData.sweep = sweep;


% This might be problematic. Emergence of the problem: having tdt not
% necessarily constructed when the GUI is made. Nevermind.
tdt.minTimeActivity = str2double(get(handles.editminTimeActivity, 'String'));
tdt.maxTimeActivity = str2double(get(handles.editMaxTime, 'String'));


% Try to read the data as an array of values between 1 and 128
try 
    channelsStr = get(handles.arrayChannelsToPlot,'String');
    switch class(channelsStr)
        case 'cell',  channels =  eval(channelsStr{1});
        case 'char',  channels =  eval(channelsStr);
    end
    
catch ME
    var = get(handles.arrayChannelsToPlot,'String');
    % Check whether this is a variable in base workshop
    if evalin('base', ['exist(''' var ''',''var'')'])
        try
            channels = evalin('base', get(handles.arrayChannelsToPlot,'String'));
        catch ME2
            disp(ME);
            disp(ME2);
            return;
        end
    end
end

if any(isnan(channels))
    fprintf('Some channels given for plotting are NaN\n');
    channels = channels(~isnan(channels));
end

chanMFloor = channels-floor(channels)>0;
if any(chanMFloor)
    fprintf('Some channels given for plotting are not integer\n');
    channels = channels(chanMFloor==0); 
end


% We also check there is a channels_filt_b
if get(handles.checkFilter, 'Value') 
    b = tdt.channels_filt_b; 
    a = tdt.channels_filt_a;  
    recChannel = cellfun(@(x)filter(b,a,x),recChannel,'UniformOutput', false);
end

% We plot recordings in axes1
handles.figureguideChannelsPlot.IntegerHandle = 'On';
handles.figureguideChannelsPlot.HandleVisibility = 'On';
set( handles.figureguideChannelsPlot, 'NumberTitle', 'On');

%handles.figureguideChannelsPlot.Number = 2;
set(0, 'CurrentFigure',  handles.figureguideChannelsPlot );
set(handles.figureguideChannelsPlot, 'CurrentAxes',  handles.axes1 );

XTime = (1:length(recChannel{1}))/tdt.RZ_SamplingFreq;
for kk=1:length(channels)
    chan = channels(kk);
    
    if kk==1
        hold off;
    else 
        hold on;
    end
    
    maxAChan = max(recChannel{chan});
    if maxAChan > 0
        plot(XTime, kk+recChannel{chan}/(2*max(abs(recChannel{chan}))),'DisplayName', num2str(chan));
    else
        plot(XTime, kk+recChannel{chan},'DisplayName', num2str(chan));
    end
end

% If there's a 'sweep', we plot the signal_left
textToDisplay = '';
if isprop(sweep,'text')
   textToDisplay = sweep.text;
end 
if isprop(tdt,'runningWavIndex')
    textToDisplay = sprintf('%s  runningWavIndex=%d/%d', ...
        textToDisplay, tdt.runningWavIndex,length(tdt.listWavFiles));
end
    

if exist('sweep', 'var')
    kk = kk + 1;
    hold on;
    plot(sweep.stimulus_delay/1000 + (1:length(sweep.signal))/sweep.sound_sampling_fq, ...
        kk+sweep.signal_left/(2*max(abs(sweep.signal_left))),'k');
    if isprop(sweep, 'text')
        text(0.05, kk+0.65, textToDisplay);
    end
end

set(handles.axes1,'xlim',[XTime(1) XTime(end)]);
set(handles.axes1,'YTick', 1:length(channels));
set(handles.axes1,'YTickLabel',channels);
set(handles.axes1,'ylim',[0.5 kk+0.8]);

% Plot average activity in axesAllChannelsActivity
% indicesAtWhichWeCalculateActivity
indCalcActMin = floor(max(tdt.minTimeActivity/1000*tdt.RZ_SamplingFreq, 1));
indCalcActMax = floor(min(tdt.maxTimeActivity/1000*tdt.RZ_SamplingFreq, length(recChannel{1})));
arrayWithOnes = zeros(size(recChannel{1}));
arrayWithOnes(indCalcActMin:indCalcActMax) = 1;

% Plot the time limits at which we compute the activity
cMaxMin = max(tdt.minTimeActivity/1000, 0);
plot([cMaxMin cMaxMin],[0 length(channels)+0.5], '-k' )
cMinMax = min(tdt.maxTimeActivity/1000, length(recChannel{1})/tdt.RZ_SamplingFreq);
plot([cMinMax cMinMax],[0 length(channels)+0.5], '-k' )


set(handles.figureguideChannelsPlot, 'CurrentAxes', handles.axesAllChannelsActivity );
plot(cellfun(@(x)norm(x.*arrayWithOnes,1),recChannel), 1:length(recChannel), '.', 'MarkerSize', 12);
ylim([0 length(recChannel)+1]);


% Refresh data 
drawnow();
guidata(hObject, handles);


% Don't check pause or stop while System3 GUI is idle
d = findall(0,'Tag', 'figureGUISetupSystem3');
isIdle = 0;
if ~isempty(d)
   isIdle = get(findall(d,'Tag', 'radioLazy'), 'Value'); 
end


% Pausing
if get(handles.radioPause, 'Value') && ~isIdle
    
    tdt.pauseTheLoop = 1;
    defaultAnswer = 'Continue';
    answer = questdlg('Continue or Stop?', 'Have a break', 'Continue', 'Stop',defaultAnswer);
    tdt.pauseTheLoop = 0;
    
    switch answer
        case 'Continue', set(handles.radioPlay, 'Value', 1);
        case 'Stop',     set(handles.radioStop, 'Value', 1); % Launches the checking next
    end
    
end


% Stopping?
if get(handles.radioStop, 'Value') && ~isIdle
    checkWeDoWantToStop(hObject, eventdata, handles);
end

% Read comments
if ~isempty(get(handles.editComment, 'String')) && ~isIdle
    editComment_Callback(hObject, eventdata, handles);
end

% Refresh GUI
drawnow();


function stopTheLoop = checkWeDoWantToStop(hObject, eventdata, handles)
% pos = get(handles.figureguideChannelsPlot, 'Position');
defaultAnswer = 'Yes'; % can't position it easily
answer = questdlg('Are you sure you want to stop the loop?', 'Check you do want to stop', 'Yes', 'No',defaultAnswer);
stopTheLoop = 0;

switch answer
    case 'Yes'
        % Break the loop
        stopTheLoop = 1;
        handles.figureguideChannelsPlot.UserData.tdt.stopTheLoop = stopTheLoop; 
        runningWavIndex = handles.figureguideChannelsPlot.UserData.tdt.runningWavIndex; 
        fprintf('checkWeDoWantToStop: Stopping ''For'' loop at tdt.runningWavIndex=%d\n',runningWavIndex);
    case 'No'
        stopTheLoop = 0;
        % Uncheck and continue
        set(handles.radioPlay, 'Value', 1);
end


function arrayChannelsToPlot_Callback(hObject, eventdata, handles)
% hObject    handle to arrayChannelsToPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of arrayChannelsToPlot as text
%        str2double(get(hObject,'String')) returns contents of arrayChannelsToPlot as a double
 pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function arrayChannelsToPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to arrayChannelsToPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxPause.
function checkboxPause_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxPause
if get(hObject, 'Value')==0
    return;
end

% pushbutton1_Callback(hObject, eventdata, handles)
% pos = get(handles.figureguideChannelsPlot, 'Position');
% answer = dialog('Position', [(pos(1)+pos(3)/2-250/2) (pos(2)+pos(4)/2-150/2) 250 150],...
%     'Name','Pause dialog Box', 'WindowStyle', 'normal');
% uicontrol('Parent', answer, 'Style', 'text', 'Position', [20 80 210 40], 'String', 'Press me to continue');
% uicontrol('Parent', answer, 'Callback', 'delete(gcf);', 'Position', [85 20 70 25], 'String', 'Close');
% uiwait(gcf)

% Set value back to 0
% set(hObject, 'Value', 0);




% --- Executes on button press in checkFilter.
function checkFilter_Callback(hObject, eventdata, handles)
% hObject    handle to checkFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkFilter
 pushbutton1_Callback(hObject, eventdata, handles)



% --- Executes on button press in checkboxStop.
function checkboxStop_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 

% Only act when clicking on it
if get(handles.checkboxStop, 'Value')==0
    return;
end

% Hint: get(hObject,'Value') returns toggle state of checkboxStop
% pushbutton1_Callback(hObject, eventdata, handles)
% set(handles.checkboxStop, 'Value', 1);
% 
% % uiwait(handles.figureguideChannelsPlot)
% defaultAnswer = 'Yes'; 
% answer = questdlg('Are you sure you want to stop the loop?', 'title', 'Yes', 'No',defaultAnswer);
% 
% switch answer
%     case 'Yes' % Break the loop
%         evalin('base','tdt.stopTheLoop=1;');
%         runningWavIndex = evalin('base', 'tdt.runningWavIndex');
%         fprintf('''For'' loop stoppedd at tdt.runningWavIndex=%d\n',runningWavIndex);
%     case 'No'  % Uncheck and continue
%         set(handles.checkboxStop, 'Value', 0);
% end
% uiresume(handles.figureguideChannelsPlot);


function radioBackOnStopIfIdleGUI(hObject, eventdata, handles)
% To inactivate clicking on pause/play/stop when the GUI is idle

% Default: do nothing
isIdle = 0;

d = findall(0, 'Tag', 'figureGUISetupSystem3');
if ~isempty(d)
    isIdle = get(findall(d, 'Tag', 'radioLazy'), 'Value');
end

if isIdle
    set(handles.radioStop, 'Value', 1);
end


% --- Executes on button press in radioStop.
function radioStop_Callback(hObject, eventdata, handles)
% hObject    handle to radioStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioStop
radioBackOnStopIfIdleGUI(hObject, eventdata, handles)

% --- Executes on button press in radioPause.
function radioPause_Callback(hObject, eventdata, handles)
% hObject    handle to radioPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioPause
radioBackOnStopIfIdleGUI(hObject, eventdata, handles)

% --- Executes on button press in radioPlay.
function radioPlay_Callback(hObject, eventdata, handles)
% hObject    handle to radioPlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioPlay
radioBackOnStopIfIdleGUI(hObject, eventdata, handles)

% --- Executes on button press in checkboxFprintfInfo.
function checkboxFprintfInfo_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxFprintfInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxFprintfInfo



function editMaxTime_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxTime as text
%        str2double(get(hObject,'String')) returns contents of editMaxTime as a double

% Check whether there is something to plot
if isequal(handles.figureguideChannelsPlot.UserData.recChannel, 0)
    return;
end

try
    val = str2double(get(hObject, 'String'));
catch 
    warning('Value not recognised as double')
    set(hObject, 'String',handles.figureguideChannelsPlot.UserData.tdt.maxTimeActivity);
    return;
end

if val<handles.figureguideChannelsPlot.UserData.tdt.minTimeActivity;
    warning('The min time should be smaller than the max.');
    set(hObject, 'String',handles.figureguideChannelsPlot.UserData.tdt.maxTimeActivity);
    return;
end
handles.figureguideChannelsPlot.UserData.tdt.maxTimeActivity = val;
pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editMaxTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editminTimeActivity_Callback(hObject, eventdata, handles)
% hObject    handle to editminTimeActivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editminTimeActivity as text
%        str2double(get(hObject,'String')) returns contents of editminTimeActivity as a double

% Check whether there is something to plot
if isequal(handles.figureguideChannelsPlot.UserData.recChannel, 0)
    return;
end

try
    val = str2double(get(hObject, 'String'));
catch 
    warning('Value not recognised as double')
    set(hObject, 'String',handles.figureguideChannelsPlot.UserData.tdt.minTimeActivity);
    return;
end
if val>handles.figureguideChannelsPlot.UserData.tdt.maxTimeActivity;
    warning('The min time should be smaller than the max.');
    set(hObject, 'String',handles.figureguideChannelsPlot.UserData.tdt.minTimeActivity);
    return;
end
handles.figureguideChannelsPlot.UserData.tdt.minTimeActivity = val;
pushbutton1_Callback(hObject, eventdata, handles)




% --- Executes during object creation, after setting all properties.
function editminTimeActivity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editminTimeActivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function figureguideChannelsPlot_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figureguideChannelsPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

radioLazy = findall(0, 'Tag', 'radioLazy');

if ~isempty(radioLazy)
    set(radioLazy, 'Value', 1);
end



function editComment_Callback(hObject, eventdata, handles)
% hObject    handle to editComment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editComment as text
%        str2double(get(hObject,'String')) returns contents of editComment as a double

% Writes the comment in commentInMemory and removes from edit box
% guideSetupSystem3 resets commentInMemory=''; after playing a .wav

% Only when the comment is enabled
if strcmp(get(handles.editComment, 'Enable'), 'off')
    return;
end

% Switch focus (if on editComment) in order to execute everything, 
% otherwise the user needs to click out of the edit box to make his 
% 'enter' effective.
% PROBLEM!!! cObj.Tag is  editComment even if you clicked elsewhere. It's
% annoying but better than previous option. Leaving this for now.
cObj = get(handles.figureguideChannelsPlot, 'CurrentObject');
if strcmp(get(cObj,'Tag'),'editComment')
    % Based on http://uk.mathworks.com/matlabcentral/answers/8921-pull-string-out-of-edit-text-without-user-hitting-enter
    uicontrol(handles.CommentsMade) %Set focus to another control
    uicontrol(handles.editComment)  %Set it back
end
%disp(handles.editComment);

%get(obj.h.edit__cmd_window,'String') 

str = get(handles.editComment, 'String')';
if isempty(str)
    return;
end
%disp(str)
commentsMade = get(handles.CommentsMade, 'String');

% Looks nicer as a cell of strings
if ischar(commentsMade)
    set(handles.CommentsMade, 'String',{commentsMade})
    commentsMade  = get(handles.CommentsMade, 'String');
end
cComment = sprintf('%s (%s)',str, datestr(now));

% Aggregates comments in memory until they're all saved, and reset from guideSetupSystem3.m 
cInMemory = handles.figureguideChannelsPlot.UserData.commentInMemory;
if isempty(cInMemory)
    commentInMemory = sprintf('%s',     cComment);
else
    commentInMemory = sprintf('%s\n%s', cInMemory, cComment);
end
handles.figureguideChannelsPlot.UserData.commentInMemory = commentInMemory;
set(handles.editComment, 'String', '');

% Concatenate
commentsMade = [commentsMade; cComment];

set(handles.CommentsMade, 'String', commentsMade);

% Set focus back to where it was?
% uicontrol(handles.(get(cObj,'Tag')));

% To tell the GUI
setappdata(0,'callbackEffectif',1);


% --- Executes during object creation, after setting all properties.
function editComment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editComment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on editComment and none of its controls.
function editComment_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to editComment (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return') || strcmp(eventdata.Key, 'enter')
    editComment_Callback(hObject, eventdata, handles)
end



% --- Executes on key press with focus on figureguideChannelsPlot or any of its controls.
function figureguideChannelsPlot_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figureguideChannelsPlot (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% This function is used whenever the focus is on the figure or any of its
% children.

% jEditbox = findjobj(handles.CommentsMade);
% if jEditbox.isFocusOwner()
%     disp('FocusOwner');
%     % set focus on figure?
% end

% If the focus is currently on the textbox, use this 

if ( strcmp(eventdata.Key, 'return') || strcmp(eventdata.Key, 'enter') ) %&&...
        % ~isempty(get(handles.editComment, 'String'))
        editComment_Callback(hObject, eventdata, handles)
end

% 
% % Only use when the comment box is not empty
% if ~isempty(get(handles.editComment, 'String')) && ...
%         ( strcmp(eventdata.Key, 'return') || strcmp(eventdata.Key, 'enter') )
%     %editComment_KeyPressFcn(hObject, eventdata, handles)
% end


% --- Executes during object creation, after setting all properties.
function figureguideChannelsPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figureguideChannelsPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Blocks access to figureGUISetupSystem3
d = findall(0, 'Tag', 'figureGUISetupSystem3');
if ~isempty(d)
    UserData = get(d, 'UserData');
    eval(UserData.CreateFcn);
end

% --- Executes when user attempts to close figureguideChannelsPlot.
function figureguideChannelsPlot_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figureguideChannelsPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% empty the comments
handles.figureguideChannelsPlot.UserData.commentInMemory = '';

% If there is a guideSetupSystem3, set radio back to lazy
dradioLazy = findall(0,'Tag','radioLazy');
if ~isempty(dradioLazy)
    set(dradioLazy,'Value',1)
end
    
% Hint: delete(hObject) closes the figure
% delete(hObject);% replaced by CloseRequestFcn from figureGUISystem3

% Unfreeze figureGUISetupSystem3
d = findall(0, 'Tag', 'figureGUISetupSystem3');
if ~isempty(d)
    UserData = get(d, 'UserData');
    eval(UserData.CloseRequestFcn);
end
