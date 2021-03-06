function RP = Circuit_Loader_device(varargin)
% CIRCUIT_LOADER Loads a *.rcx circuit onto a RP2, returns ActiveX control object
% RP = CIRCUIT_LOADER(connectionType, deviceNumber, circuitPath)
%     User input require if no inputs to function
%     Options for connectionType are 'GB' and 'USB'
%     circuitPath does not require an extension.
%     connectionType defaults to 'GB'
%     deviceNumber defaults to 1
%     Note: code must be modified to work with non-RP2 devices or with *.rco files
%
%     Alban: Adapted from original, to accept type of device (RX8)

if nargin == 4
  if ~ismember(varargin{4},{'RX8','RZ2','RP2'})
    error('Given device (%s) is not recognised. Add this device in list of Circuit_Loader_device.m',varargin{4});
  end
  connectDevEnd = varargin{4};
else
  connectDevEnd = 'RX8'; % ConnectRP2 originally
end
connectDev = ['Connect' connectDevEnd];

if nargin >= 3
  
  connectionType = varargin{1};
  deviceNumber = varargin{2};
  circuitPath = varargin{3};
  
elseif nargin == 1
  
  connectionType = 'GB';
  deviceNumber = 1;
  circuitPath = varargin{1};
  
elseif nargin == 0
  
  % path - set this to wherever the examples are stored
  path = 'Path/to/RP_files';
  
  connectionType = input('Enter the type of connection (USB or GB):  ','s');
  connectionType = upper(connectionType);
  
  % Error check for correct connection
  if ~(strcmp(connectionType,'USB')||strcmp(connectionType,'GB'))
    connectionType = 'GB';
    disp('   Device connection = GB');
  end
  
  deviceNumber = input('Enter the device number:  ');
  
  % Error check for correct device number
  if (~isnumeric(deviceNumber) || deviceNumber < 1)
    deviceNumber = 1;
    disp('   Device number = 1');
  end
  
  % Show available circuits
  disp(' ');
  disp(['path: ' path]);
  dir(path)
  
  circuitPath = input('Enter the name of the circuit:  ','s');
  circuitPath = strcat(path,circuitPath);
  
else
  error('Invalid number of arguments.');
end

% Error check circuit file path
if size(strfind(circuitPath,'.rcx')) == 0
  circuitPath = strcat(circuitPath,'.rcx');
end

% Error check for existing file
fileExists=(exist(circuitPath,'file'));
if fileExists==0
  disp('   File doesnt exist'); return;
end

% Load circuit onto device and run
RP = actxcontrol('RPco.x',[5 5 26 26]); 

RP.(connectDev)(connectionType, deviceNumber); % Connects device via USB or GB given the proper device number
RP.Halt; % Stops any processing chains running on device
RP.ClearCOF; % Clears all the buffers and circuits on that device
disp(['Loading ' circuitPath '...']);
RP.LoadCOF(circuitPath); % Loads circuit. Bif SerialBuf size leads to error
RP.Run; % Starts circuit

status=double(RP.GetStatus); % Gets the status
if bitget(status,1)==0; % Checks for connection
  disp(['Error connecting to ' connectDevEnd]); return;
elseif bitget(status,2)==0; % Checks for errors in loading circuit
  disp('Error loading circuit'); return;
elseif bitget(status,3)==0 % Checks for errors in running circuit
  disp('Error running circuit'); return;
else
  disp('Circuit loaded and running'); return;
end
end

