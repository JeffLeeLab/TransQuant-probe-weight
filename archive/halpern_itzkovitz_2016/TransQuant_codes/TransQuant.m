function varargout = TransQuant(varargin)
% TransQuant M-file for TransQuant.fig
%
%      TransQuant, by itself, creates a new TransQuant or raises the existing
%      singleton*.
%
%      H = TransQuant returns the handle to a new TransQuant or the handle to
%      the existing singleton*.
%
%      TransQuant('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TransQuant.M with the given input arguments.
%
%      TransQuant('Property','Value',...) creates a new TransQuant or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TransQuant_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TransQuant_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TransQuant

% Last Modified by GUIDE v2.5 09-Nov-2015 17:05:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @TransQuant_OpeningFcn, ...
    'gui_OutputFcn',  @TransQuant_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
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

% --- Executes just before TransQuant is made visible.
function TransQuant_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TransQuant (see VARARGIN)

% Choose default command line output for TransQuant
handles.output = hObject;
% compatibility with R2014
xx=ver;nms={xx(:).Name};ind=find(strcmp(lower(nms),'matlab'));
if isempty(findstr(lower(xx(1).Release),'2012')) & isempty(findstr(lower(xx(1).Release),'2011')) & isempty(findstr(lower(xx(1).Release),'2011'))  & isempty(findstr(lower(xx(1).Release),'2013'))
    hObject.UserData=[];
end
% This sets up the initial page
if strcmp(get(hObject,'Visible'),'off')
    % set windows title
    set(gcf,'Name','TransQuant: Software for extracting transcription and degradation rates')
    plot(1)
    axis off;
    % message
    title({'TransQuant - Single molecule FISH estimation of',' transcription and degradation rates'},'FontSize',20,'Color','b');drawnow;
    %text(0.3,1.8,{'transquant is a MATLAB package with GUI for RNA FISH image analysis,',...
    %   },...
    %    'HorizontalAlignment','left')
end


% initialization of user dadta
Initialization(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = TransQuant_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- initialize UserData
function Initialization(hObject, eventdata, handles)
UserData=get(handles.figure1,'UserData');
%UserData=guidata(handles.figure1)
UserData.stack_range=[1 1000];% stack filter [min,max]

%status: indicator of progress
UserData.status.loaded = 0;
UserData.status.projected = 0;
UserData.status.enhanced = 0;
UserData.status.cropped = 0 ;
UserData.status.segmented = 0;
UserData.status.counted = 0;
UserData.status.saved = 0;

UserData.status.changed = 0; % data changed, need save now

% laplacian filter
UserData.N= 15;%
UserData.sigma = 1.5;%


% Size Threshold to filter dots
UserData.SIZE_THRESHOLD_FILTER=25;

% LOG weight factor
UserData.SHARP_FACTOR=1;

% when displaying dot images filter outliers (0/1)
UserData.FILTER_OUTLIER=0;

% LOG saturation Zscore 
UserData.ZTHRESH=-1;

% Fast mode option
UserData.FAST_MODE=0;

% # of intensity thresholds when counting dots
UserData.n_thresholds = 100;

%auto-thresholding parameters
UserData.auto_thresholding.width=5;
UserData.auto_thresholding.offset = 10;

%default projection method
UserData.projection_method='max';

% default enhance method
UserData.enhance_method = 'imadjust';

% cell information
UserData.cell = [];
UserData.divide_lines={};
UserData.nuc_contour=[];
UserData.dapi_areas=[];
UserData.dapi_areas2=[];
UserData.cell_subset=[];
UserData.pclass2=[];
UserData.pclass3=[];
% struct array
% cell.label: a character.
% cell.edge: N-by-2 matrix, each row is a point. together is a polygon.
% cell.center: [x,y], calculated as the mean of cell.edge
% cell.area: area of the cell defined by cell.edge

% dot information
UserData.dots = [];
UserData.nuclear_dots = [];
% for column
% 1-3: 3D position
% 4: channel
% 5: cell
% 6: outline projection

% original image file, multiple stacks
UserData.I = [];

% projected images
UserData.I2 = [];

% enhanced/cropped image
UserData.I3 = [];

% scale compared to cropped image
UserData.scale =[1 1];

% size of original image
UserData.image_size = [];

% save path and file name
UserData.save_path = [];
UserData.save_name = [];

% crop region
UserData.BW = [];
UserData.RECT = [];

% reference image
UserData.reference = [];

% skeleton for relabel cells
UserData.skeleton = [];

% channels that are counted
UserData.channel_name = [];
UserData.nuclear_channel_name = [];


% dots summary
UserData.num_dot_cell_channel = [];

% each row is a cell
% each column is a channel

% normalized dots summary
UserData.num_dot_cell_channel_normalized = [];

% outline for projecting dots
UserData.outline = [];

% neighboring cells
UserData.neighbors = [];
% each row contains the indices of two cells that are defined as neighbors

UserData.version = 2;

set(handles.figure1,'UserData',UserData);
%guidata(handles.figure1,'UserData');
% enable all commands, responds to menu
% --------------------------------------------------------------------
function enable_all_Callback(hObject, eventdata, handles)
% hObject    handle to enable_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.previous,'Enable','on');% toolbar button, previous stack
set(handles.next,'Enable','on');% toolbar button, next stack
set(handles.projection,'Enable','on');% toolbar button, projection
%set(handles.crop,'Enable','on');% toolbar button, crop image
set(handles.save,'Enable','on');% toolbar button, save
set(handles.save_menu,'Enable','on'); % menu item,
%set(handles.ProjectionMenu,'Enable','on');% menu item, projection
set(handles.open_saved,'Enable','on');% menu item, open saved results
set(handles.image,'Enable','on');% menu item, image enhance
set(handles.image_enhance,'Enable','on');% toolbar button, image enhance
set(handles.segment,'Enable','on');%toolbar button, segment cells
set(handles.count,'Enable','on');% toolbar, count dots
set(handles.count_dots,'Enable','on') %menu
set(handles.Analyze,'Enable','on');% menu, analyze
set(handles.compare_neighbor,'Enable','on');%menu,analyze->neighbor cell->
set(handles.checkbox3,'Visible','on')% dots
set(handles.checkbox4,'Visible','on')% cells
set(handles.checkbox5,'Visible','on')% outline
set(handles.checkbox6,'Visible','on')% projection
set(handles.resetcelllabel,'Enable','on')% menu
set(handles.plot_counts,'Enable','on')%menu
set(handles.neighbor_cell,'Enable','on')%menu
set(handles.detect_bg_dots,'Enable','on');%menu
set(handles.uipanel1,'Visible','on')%show channel
set(handles.plot_cryptogram,'Enable','on');%menu

% enable components under different conditions to ensure the right order of
% actions
% --------------------------------------------------------------------
function enable_components(UserData,handles)
if UserData.status.loaded
    set(handles.previous,'Enable','on');
    set(handles.next,'Enable','on');
    set(handles.projection,'Enable','on');
    %set(handles.ProjectionMenu,'Enable','on');
    set(handles.open_saved,'Enable','on');
    set(handles.image,'Enable','on')
    set(handles.Analyze,'Enable','on');
end
if UserData.status.projected
    set(handles.crop,'Enable','on');
    set(handles.save_menu,'Enable','on');
    set(handles.image,'Enable','on');
    set(handles.image_enhance,'Enable','on');
    set(handles.count,'Enable','on');
    set(handles.count_dots,'Enable','on');
    set(handles.segment,'Enable','on');
end
if UserData.status.cropped
    set(handles.crop,'Enable','off');
    set(handles.projection,'Enable','off');
    set(handles.ProjectionMenu,'Enable','off');
end

if UserData.status.counted
    set(handles.image_enhance,'Enable','off');
    set(handles.image,'Enable','off');
end
if UserData.status.changed
    set(handles.save,'Enable','on');
else
    set(handles.save,'Enable','off');
end

if ~isempty(UserData.dots) | ~isempty(UserData.nuclear_dots)
    set(handles.detect_bg_dots,'Enable','on');
    set(handles.uipanel1,'Visible','on')
    set(handles.checkbox3,'Visible','on')
    set(handles.plot_cryptogram,'Enable','on');
    set(handles.Analyze,'Enable','on');
    if ~isempty(UserData.cell)
        set(handles.plot_counts,'Enable','on')
        set(handles.neighbor_cell,'Enable','on')
    else
        set(handles.plot_counts,'Enable','off')
        set(handles.neighbor_cell,'Enable','off')
    end
else
    set(handles.checkbox3,'Visible','off')
end

if isempty(UserData.cell)
    set(handles.checkbox4,'Visible','off')
else
    set(handles.checkbox4,'Visible','on')
    set(handles.resetcelllabel,'Enable','on')
end

if isempty(UserData.outline)
    set(handles.checkbox5,'Visible','off')
    set(handles.checkbox6,'Visible','off')
else
    set(handles.checkbox5,'Visible','on')
    set(handles.checkbox6,'Visible','on')
end

if isempty(UserData.neighbors)
    set(handles.compare_neighbor,'Enable','off');
else
    set(handles.compare_neighbor,'Enable','on');
end


% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ask user to confirm for closing the software
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    % exit if canceled
    return;
end

% if confirmed, run the following code

UserData=get(handles.figure1,'UserData');

if UserData.status.changed % if some changes have been made but not saved
    choice = questdlg('Do you want to save your data?',...
        'Save your data',...
        'Yes',...
        'No',...
        'Cancel',...
        'Yes');
    switch choice
        case 'Cancel'
            return;
        case 'No'
            delete(handles.figure1)
        case 'Yes'
            
            % if never saved
            if isempty(UserData.save_name)
                
                % remember current directory
                current_path = pwd;
                s = regexp(UserData.path, '\\|\/', 'split');
                default_file_name = [s{end-1} '.' UserData.file_index '.ima'];
                
                % set default path as UserData.path
                cd(UserData.path);
                [save_name, save_path] = uiputfile({'*.ima';'*.*'},'Please specify a file name for automatic save service',default_file_name);
                
                cd(current_path);
                
                if save_name==0 return; end %canceled
                
                UserData.save_name = save_name;
                UserData.save_path = save_path;
                UserData.status.saved = 1;
                
            else
                set(handles.figure1,'UserData',UserData);
                save([UserData.save_path UserData.save_name],'UserData','-mat');
            end
            delete(handles.figure1)
    end
end

% open new image file
% --------------------------------------------------------------------
function open_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current data
UserData=get(handles.figure1,'UserData');

% if a file is already opened, warn the user
if UserData.status.loaded
    selection = questdlg('Are you sure? This will discard current data.','Warnings!', ...
        'Yes','No','Yes');
    if strcmp(selection,'No')
        return;
    end
end

% open file choosing dialog
if isfield(UserData,'name')
    [UserData.name, UserData.path, filter_index] = uigetfile('*.tif','Open an image file (.tif)',[UserData.path UserData.name]);
else
    [UserData.name, UserData.path, filter_index] = uigetfile('*.tif','Open an image file (.tif)');
end
% debug
UserData.name

% exit if no file is choosed
if UserData.name==0 return; end

% clear figure
cla
legend('off')
title('')

% save file path and name
set(handles.figure1,'UserData',UserData);

% initialization of Non-GUI generated user data: status and parameters
Initialization(hObject, eventdata, handles);

% get all data, particularly GUI-generated data
UserData=get(handles.figure1,'UserData');

% enable/disable component according to status, doesn't change UserData
enable_components(UserData,handles);

% extract file index
% xxx001.tif, extract 001
%UserData.file_index=UserData.name((end-6):(end-4));
ind1=findstr(UserData.name,'_');
ind2=findstr(UserData.name,'.');
UserData.file_index=UserData.name((ind1+1):(ind2-1));

% show message panel
set(handles.uipanel3,'Visible','on');drawnow;
% display message
set(handles.text1,'String',['Loading file ' UserData.name ', please wait ...']);drawnow;

% loading image file
%UserData.I=parse_stack([UserData.path UserData.name],UserData.stack_range(1),UserData.stack_range(2));

% COMMENT - changed on 5/1/2014 by SHALEV to normalize for bleaching before
% max-projection of the reference image
UserData.I=parse_stack_bleach_normalize([UserData.path UserData.name],UserData.stack_range(1),UserData.stack_range(2),UserData.FILTER_OUTLIER);

UserData.image_size = size(UserData.I);
if length(UserData.image_size)==2,
    UserData.image_size(3)=1;
end
% update stack_range
if UserData.stack_range(2) > UserData.image_size(3)
    UserData.stack_range(2) = UserData.image_size(3);
end

% update windows name
set(handles.figure1,'Name',UserData.name);

% display the first stack
UserData.current_stack = 1;
imagesc(UserData.I(:,:,UserData.current_stack));colormap gray;axis square;axis off;

set(handles.text1,'String',{'Stack 1','(check other stacks by mouse scroll wheel or arrows on toolbar)'})

% update status and other
UserData.status.loaded = 1;
UserData.status.changed = 1;

% update component status
enable_components(UserData,handles);

% save data
set(handles.figure1,'UserData',UserData);

%debug
UserData.status


% open an image from menu
% --------------------------------------------------------------------
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open_ClickedCallback(hObject, eventdata, handles)


% open a saved results only.
% --------------------------------------------------------------------
function open_saved_Callback(hObject, eventdata, handles)
% hObject    handle to open_saved (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current data
UserData=get(handles.figure1,'UserData');

% if a file is already opened, warn the user
if UserData.status.loaded
    selection = questdlg('Are you sure? This will discard current data.','Warnings!', ...
        'Yes','No','Yes');
    if strcmp(selection,'No')
        return;
    end
end


% open file choosing dialog
if ~isempty(UserData.save_name)
    [name, path, filter_index] = uigetfile('*.ima','Open saved results',[UserData.save_path UserData.save_name]);
else
    [name, path, filter_index] = uigetfile('*.ima','Open saved results');
end
if name==0 return; end

% load file
load([path name],'-mat');


UserData.path = path; % UserData.path is where the program looks for the other files (dot counting, segmentation reference)
UserData.save_path = path; % UserData.save_path is where the program saves the ima file when save is clicked.
UserData.save_name = name;

if ~isfield(UserData,'SIZE_THRESHOLD_FILTER'),
    UserData.SIZE_THRESHOLD_FILTER=25;
end

if ~isfield(UserData,'SN_THRESHOLD_FILTER'),
    UserData.SN_THRESHOLD_FILTER=2;
end

if ~isfield(UserData,'SHARP_FACTOR'),
    UserData.SHARP_FACTOR=1;
end

if ~isfield(UserData,'ZTHRESH'),
    UserData.ZTHRESH=-1;
end


if ~isfield(UserData,'FAST_MODE'),
    UserData.FAST_MODE=0;
end

if ~isfield(UserData,'FILTER_OUTLIER'),
    UserData.FILTER_OUTLIER=1;
end

% set data to gui
set(handles.figure1,'UserData',UserData);

% for debug, check status
UserData.status

% clear figure
cla;
legend('off')
title('')

% resize cropped image
if ~isempty(UserData.image_size)
    I3 = imresize(UserData.I2, UserData.image_size(1:2));
else
    I3=UserData.I2;
end

% show image
imagesc(I3);
axis off;axis square;colormap gray;

% change windows title
set(handles.figure1,'Name',name);

set(handles.uipanel3,'Visible','on')

% show message
set(handles.text1,'String',name)

% for data saved via previous version
if ~isfield(UserData,'version')
    warndlg('this data is saved via a previous version of the software. It will be transformed to the current version without loss of data.',...
        'Version change',...
        'modal');
    if isfield(UserData,'scale')
        UserData.scale=1./UserData.scale;
        if isfield(UserData,'dots')
            if ~isempty(UserData.dots)
                UserData.dots(:,1)=UserData.dots(:,1)/UserData.scale(1);
                UserData.dots(:,2)=UserData.dots(:,2)/UserData.scale(2);
            end
        end
        if isfield(UserData,'cell')
            if ~isempty(UserData.cell)
                for i=1:length(UserData.cell)
                    UserData.cell(i).center = UserData.cell(i).center./UserData.scale;
                    UserData.cell(i).edge(:,1)=UserData.cell(i).edge(:,1)./UserData.scale(1);
                    UserData.cell(i).edge(:,2)=UserData.cell(i).edge(:,2)./UserData.scale(2);
                end
            end
        end
    end
    if ~isfield(UserData,'name')
        UserData.name = name;
    else
        UserData.name =[];
    end
    
    UserData.version = 2;
end

if isfield(UserData,'border')
    UserData.reference = UserData.border;
end
if isfield(UserData,'reference')
    if ~isempty(UserData.reference)
        if ~isempty(UserData.image_size)
        UserData.reference = imresize(UserData.reference,UserData.image_size(1:2));
        end
    end
else
    UserData.reference = [];
end

if isfield(UserData.status,'changed')
    UserData.status.changed = 1;
else
    UserData.status.changed = 0;
end

if ~isfield(UserData,'outline')
    UserData.outline = [];
end

if ~isfield(UserData,'neighbors')
    UserData.neighbors = [];
end

% display dots if counted
if ~isfield(UserData,'dots')
    UserData.dots = [];
end

% display dots if counted
if ~isfield(UserData,'nuclear_dots')
    UserData.nuclear_dots = [];
end
plot_all_dots(UserData);

% display cells if segmented
if isfield(UserData,'cell')
    if ~isempty(UserData.cell)
        % assign labels if none assigned, for compatability with previous
        % version
        if ~isfield(UserData.cell,'label')
            UserData.cell = set_cell_label(UserData.cell);
        end
        plot_all_cell(UserData.cell);
    end
else
    UserData.cell=[];
end

axis off;axis square;colormap gray;

UserData.status.changed = 1;
% update component status
enable_components(UserData,handles);
% save data
set(handles.figure1,'UserData',UserData);

% open a reference image
% --------------------------------------------------------------------
function open_image_Callback(hObject, eventdata, handles)
% hObject    handle to open_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');

% open file choosing dialog
[UserData.refname, UserData.refpath] = uigetfile('*.tif','Open a reference image file (.tif)');

% exit if no file is choosed
if UserData.refname==0 return; end

% load image file
UserData.reference=parse_stack([UserData.refpath UserData.refname],UserData.stack_range(1),UserData.stack_range(2),UserData.FILTER_OUTLIER);

set(handles.figure1,'UserData',UserData);
panel(handles)


% save, respond to the button on the toolbar
% --------------------------------------------------------------------
function save_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to save_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');
set(handles.text1,'String','Saving data ...');drawnow
assign_dots_to_cells(handles);
UserData=get(handles.figure1,'UserData');
%NormTotalCount_Callback_no_plot(hObject, eventdata, handles)
UserData=get(handles.figure1,'UserData');
% if this is the first time to save, let user specify a file name
if isempty(UserData.path),
    return;
end
if isempty(UserData.save_name)
    
    % remember current directory
    current_path = pwd;
    s = regexp(UserData.path, '\\|\/', 'split');
    default_file_name = [s{end-1} '.' UserData.file_index '.ima'];
    
    % set default path as UserData.path
    cd(UserData.path);
    [save_name, save_path] = uiputfile({'*.ima';'*.*'},'Please specify a file name for automatic save service',default_file_name);
    
    cd(current_path);
    
    if save_name==0 return; end %canceled
    
    UserData.save_name = save_name;
    UserData.save_path = save_path;
    UserData.status.saved = 1;
end


% do not save the original stack image, it's huge! Instead keep the summed
% image over Z sections - this will be used for Dapi measurements of cells
I = UserData.I;
%UserData.I = [];
UserData.I = sum(UserData.I,3);
save([UserData.save_path UserData.save_name],'UserData','-mat');
UserData.I=I;


set(handles.figure1,'UserData',UserData);

figure(handles.figure1)
set(handles.text1,'String','Data saved');

UserData.status.changed = 0;
enable_components(UserData,handles)


% debug
UserData=get(handles.figure1,'UserData');
UserData.status

% --------------------------------------------------------------------
function save_menu_Callback(hObject, eventdata, handles)
% hObject    handle to save_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
save_ClickedCallback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function saveas_Callback(hObject, eventdata, handles)
% hObject    handle to saveas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');

% remember current directory
current_path = pwd;
s = regexp(UserData.path, '\\|\/', 'split');
default_file_name = [s{end-1} '.' UserData.file_index '.ima'];

% set default path as UserData.path
cd(UserData.path);
[save_name, save_path] = uiputfile({'*.ima';'*.*'},'Please specify a file name for automatic save service',default_file_name);

cd(current_path);

if save_name==0 return; end %canceled

UserData.save_name = save_name;
UserData.save_path = save_path;
UserData.status.saved = 1;
set(handles.figure1,'UserData',UserData);
enable_components(UserData,handles);
UserData.I = [];

set(handles.figure1,'UserData',UserData);
save([UserData.save_path UserData.save_name],'UserData','-mat');

% debug
UserData=get(handles.figure1,'UserData');
UserData.status


% for batch analysis, not implemented
% --------------------------------------------------------------------
function open_dir_Callback(hObject, eventdata, handles)
% hObject    handle to open_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = uigetdir;

if path == 0
    return;
end

files = dir(fullfile(path,'*.ima'));

UserDataArray=[];
n_file = length(files)
for i=1:n_file
    % load the file
    tmp = load(fullfile(path,files(i).name),'-mat');
    UserDataArray=[UserDataArray tmp.UserData];
    %UserDataArray(i).save_name
end


% Count dots 
% --------------------------------------------------------------------
function count_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

count_type=0;
UserData=get(handles.figure1,'UserData');

% find all channels
all_channel=all_channel_names(UserData.path,UserData.file_index);

count_type=2;

USE_LOG=2;
UserData.SN=[];


% let user choose channels to be counted
[selected,v] = listdlg('Name','Channel filter',...
    'PromptString','Select the channels you want to count the dots:',...
    'ListSize',[300 100],...
    'ListString',all_channel,'SelectionMode','single');
if isempty(selected)
    return;
end

new_channels = all_channel(selected); % user selected channels

% check whether some of the channels have already been counted
% if so, ask the user whether to count again or not

counted_channels = [];% to store the indices of channels that are already counted
counted_channels_in_channel_name = []; % indices in channel_name

% if already counted some channels
if ~isempty(UserData.channel_name)
    if ~strcmp(UserData.channel_name,new_channels)
        h = errordlg('You can only recount the same channel, if you want to count a new channel open a new project (you can still use ''load the segmented cells''))');
        return;
    end
end
if ~isempty(UserData.channel_name)
    for i=1:length(new_channels)
        for j=1:length(UserData.channel_name)
            % if there is overlab between new channel and counted channel
            if strcmp(new_channels{i},UserData.channel_name{j})
                % add to indices of counted channels
                counted_channels=[counted_channels i];
                counted_channels_in_channel_name = [counted_channels_in_channel_name j];
            end
        end
    end
    % if some channels are counted, let user to decide
    if ~isempty(counted_channels)
        choice = questdlg([new_channels(counted_channels) ' have been counted already, do you like to count them AGAIN?'],...
            'Channels already counted',...
            'Yes',...
            'No',...
            'Cancel',...
            'No');
        switch choice
            case 'No'
                % remove counted channels from the list of channels to be
                % counted
                new_channels(counted_channels)=[];
            case 'Yes'
                UserData.channel_name(counted_channels_in_channel_name) = [];
                % also need to update dots!!
                for k = 1:numel(counted_channels_in_channel_name)
                    UserData.dots(UserData.dots(:,4)==counted_channels_in_channel_name(k),:)=[];
                end
            case 'Cancel'
                return;
        end
        
    end
    UserData.channel_name = [UserData.channel_name new_channels];
    
else
    UserData.channel_name = [new_channels];    
end


% Create the combined channel name
UserData.full_channel_name={};
counter=1;
for i=1:length(UserData.channel_name),
    UserData.full_channel_name{counter}=UserData.channel_name{i};
    counter=counter+1;
end
for i=1:length(UserData.nuclear_channel_name),
    UserData.full_channel_name{counter}=['nuc-' UserData.nuclear_channel_name{i}];
    counter=counter+1;
end


set(handles.text1,'String','Starting counting dots, please wait ...')

% start counting dots in each channel
col = 'rgb'; % color for dots in different channel

if count_type==2,
    channel_start = length(UserData.channel_name) - length(new_channels)+1;
elseif count_type==1,
    channel_start = length(UserData.nuclear_channel_name) - length(new_channels)+1;
end

%debug
channel_start
new_channels
UserData.full_channel_name

% for each channel that need to count
for j=channel_start:length(UserData.channel_name),
    % load channel image
    filename = [UserData.channel_name{j} '_' UserData.file_index '.tif'];
    set(handles.text1,'String',['loading ' filename ' ...']);
    ims=parse_stack_bleach_normalize([UserData.path filename],UserData.stack_range(1),UserData.stack_range(2),UserData.FILTER_OUTLIER);
    
    % correct the channel according to the relevant shift
    
    %         if ~isempty(findstr(lower(UserData.channel_name{j}),'cy')),
    %             ims=correct_shift(ims,UserData.Cy5_shift);
    %         elseif ~isempty(findstr(lower(UserData.channel_name{j}),'a594')),
    %             ims=correct_shift(ims,UserData.A594_shift);
    %         elseif ~isempty(findstr(lower(UserData.channel_name{j}),'alexa')),
    %             ims=correct_shift(ims,UserData.A594_shift);
    %         elseif ~isempty(findstr(lower(UserData.channel_name{j}),'tmr')),
    %             ims=correct_shift(ims,UserData.TMR_shift);
    %         end
    %
    % crop image
    set(handles.text1,'String','Cropping each stack ...');drawnow
    if ~isempty(UserData.BW)
        ims=gui_crop_stack_poly(ims,UserData.RECT,UserData.BW);
    end
    
    
    %        UserData.cell_intensities{j}=extract_cell_intensity(UserData,ims);
    % calculate the number of dots for different thresholds
    
    [thresholdfn,thresholds,ims2,ims3]=calculate_file_threshold_function(handles,ims,UserData.BW,UserData.N,UserData.sigma,UserData.n_thresholds,UserData.ZTHRESH);
    [num_dots,locations,threshold,pixel_list,pixel_values,SN]=threshold_selection(handles,thresholds,thresholdfn,ims,ims2,UserData.BW,...
        UserData.auto_thresholding.width,...
        UserData.auto_thresholding.offset,...
        UserData.projection_method,...
        UserData.enhance_method,ims3,USE_LOG,UserData.SIZE_THRESHOLD_FILTER,UserData.SHARP_FACTOR,UserData.FAST_MODE);
    
    % if there are some dots in the channel
    if num_dots > 0
        % correct the coordinates
        locations(:,1) = locations(:,1)./UserData.scale(1);
        locations(:,2) = locations(:,2)./UserData.scale(2);
        
        for i=1:length(pixel_list),
            temp=pixel_list(i).PixelList;
            temp(:,1) = temp(:,1)./UserData.scale(1);
            temp(:,2) = temp(:,2)./UserData.scale(2);
            pixel_list(i).PixelList=temp;
        end
        
        if size(UserData.dots,2)==5 % dots have already been assigned to cells
            UserData.dots = [UserData.dots; [locations zeros(num_dots,1)+j NaN*ones(num_dots,1)]];
            UserData.dot_pixel_list{j}=pixel_list;
            UserData.dot_pixel_values{j}=pixel_values;
        else
            UserData.dots = [UserData.dots; [locations zeros(num_dots,1)+j]];
            UserData.dot_pixel_list{j}=pixel_list;
            UserData.dot_pixel_values{j}=pixel_values;
        end
        
        % ADD the signal to noise data
        
        UserData.SN=[UserData.SN;SN];        
             
        % save data
        set(handles.figure1,'UserData',UserData);
        save_ClickedCallback(hObject, eventdata, handles);
        
    else
        set(handles.text1,'String',['No dots in ' filename])
    end
end

UserData.status.counted  = 1;
set(handles.figure1,'UserData',UserData);

figure(handles.figure1);


plot_all_dots(UserData);

% add legend for channels
% first need to shut off the legend for cells if any
hs = findobj('-regexp','Tag','plot_cell');
if ~isempty(hs)
    for i=1:length(hs)
        set(get(get(hs(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
    end
end
legend(UserData.full_channel_name,'FontSize',14)

b=2;


% Add in the intensity of the dots
set(handles.text1,'String','Extracting dot intensities');pause(0.2)
for i=1:length(UserData.dot_pixel_list)
    UserData.dot_intensities{i}=zeros(size(UserData.dot_pixel_list{i},1),1);
    for j=1:length(UserData.dot_pixel_list{i})
        UserData.dot_intensities{i}(j)=sum_max_stack(UserData.dot_pixel_list{i}(j).PixelList,UserData.dot_pixel_values{i}(j).PixelValues);
    end
end
set(handles.text1,'String','Done');pause(0.2)
UserData.status.changed = 1;
enable_components(UserData,handles)

set(handles.figure1,'UserData',UserData);
save_ClickedCallback(hObject, eventdata, handles);


% segment cells by clicking polygons
% --------------------------------------------------------------------
function segment_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to segment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');

if isempty(UserData.reference) % first time
    %if 1,
    %debug
    disp('first time to segment')
    UserData.status
    
    % let user to choose a channel that helps segmentation
    all_channel=all_channel_names(UserData.path,UserData.file_index);
    [selected,v] = listdlg('Name','Channel filter',...
        'PromptString','Select a channel to aid region segmentation:',...
        'ListSize',[300 100],...
        'SelectionMode','single',... % select only one channel
        'ListString',all_channel);
    
    % debug
    selected
    
    % if none selected, exit
    if isempty(selected)
        return;
    end
    
    % the selected channel name
    reference_channel =all_channel(selected);
    
    UserData.reference_filename= [reference_channel{1} '_' UserData.file_index '.tif'];
    
    set(handles.text1,'String',['Loading ' UserData.reference_filename ', please wait ...'])
    
    % load the channel
    ims=parse_stack([UserData.path UserData.reference_filename],UserData.stack_range(1),UserData.stack_range(2),UserData.FILTER_OUTLIER);
    dapi_ims=ims;
    
    % z-projection
    UserData.reference = zprojection(ims,UserData.projection_method);
    if ~isempty(UserData.BW)
        UserData.reference=gui_crop_stack_poly(UserData.reference,UserData.RECT,UserData.BW);
    end
    
    %debug
    size(UserData.reference)
    
    if ~isempty(UserData.image_size)
        UserData.reference = imresize(UserData.reference, UserData.image_size(1:2));
    end
    UserData.reference=enhance_image(UserData.reference,UserData.enhance_method);
    
    if isempty(UserData.cell)
        UserData.cell=[];
    end
else
    ims=parse_stack([UserData.path UserData.reference_filename],UserData.stack_range(1),UserData.stack_range(2),UserData.FILTER_OUTLIER);
    dapi_ims=ims;
end

% change the position of the new figure window
screen_size = get(0, 'ScreenSize')
position = get(handles.figure1, 'Position')
href = figure;

set(gcf, 'Position', [position(3)+10 position(2) screen_size(3)-position(3) position(4)] );
imagesc(uint16(UserData.reference));
axis off;axis square; colormap gray;
set(handles.text1,'String','Refer to the other image when segmenting regions');

set(handles.uipanel1,'Visible','on')

% debug
size(UserData.reference)

figure(handles.figure1)
set(handles.text1,'String','Click the region border');
plot_all_cell(UserData.cell);
figure(href)
plot_all_cell(UserData.cell);

% start segmenting cells
go_on = 1;
hold on;
while go_on
    figure(handles.figure1)
    UserInput=1;
    while UserInput > 0
        UserInput = waitforbuttonpress;      % Wait for click
    end
    SelectionType = get(gcf,'SelectionType')
    CurrentPoint = get(gca,'CurrentPoint')
    if strcmp(SelectionType,'alt') % right click
        figure(handles.figure1)
        choice = questdlg('I want to ...',...
            'Tell me what you want to do:',...
            'stop',...
            'delete the cell',...
            'stop');
        switch choice
            case 'delete the cell'
                % find the cell
                found = which_cell(CurrentPoint(1,1:2),UserData.cell);
                if found==0
                    % the point is not inside any cell
                    set(handles.text1,'String','You need to right-click INSIDE the cell!');
                else
                    % delete the cell
                    UserData.cell(found)=[];
                    % update cell labels
                    UserData.cell = set_cell_label(UserData.cell);
                    % save data
                    set(handles.figure1,'UserData',UserData);
                    % debug
                    UserData.cell
                    
                    % update figure
                    figure(handles.figure1)
                    % first delete old plot
                    h = findobj('-regexp','Tag','plot_cell');
                    delete(h)
                    h = findobj('-regexp','Tag','plot_nuclear_outline');
                    delete(h)
                    % then plot all cells
                    plot_all_cell(UserData.cell);
                    
                    % update reference figure
                    figure(href)
                    plot_all_cell(UserData.cell);
                end
                
            otherwise % stop
                UserData.status.segmented = 1;
                set(handles.figure1,'UserData',UserData);
                set(handles.text1,'String',['You have selected ' num2str(length(UserData.cell)) ' cells' ])
                go_on =0;
        end
    else
        % segment cells
        figure(handles.figure1)
        % user clicks a polygon
        poly=getline(handles.figure1,'closed');
        
        if size(poly,1)>2   % ignore when less than 3 points are clicked
            % calculate the cell area
            area = polyarea(poly(:,1),poly(:,2));
            if area > 100 % ignore when the polygon is too small
                % update cell information
                newcell.edge = poly;
                newcell.center = mean(newcell.edge);
                newcell.area = area;
                newcell.label = num2str(length(UserData.cell)+1);
                newcell.dapi = NaN;
                newcell.extended_edge=NaN;
                newcell.internal_edge=NaN;
                newcell.extended_edge_new=NaN;
                newcell.rim_edge=NaN;
                if isfield(UserData.cell,'TS_list')
                    newcell.TS_list=[];
                end
                
                if isfield(UserData.cell,'channel')
                    UserData.cell=rmfield(UserData.cell,'channel');
                end
                UserData.cell
                newcell
                newcell.edge
                UserData.cell = [UserData.cell newcell];
                set(handles.figure1,'UserData',UserData);
                
                %debug
                %length(cell)
                %cell(length(cell))
                
                % plot all cells
                plot_all_cell(UserData.cell);
                figure(href)
                plot_all_cell(UserData.cell);
                title(['You have selected ' num2str(length(UserData.cell)) ' regions' ],'FontSize',18,'Color','b');
                
                figure(handles.figure1)
                set(handles.text1,'String',['Please click to crop the next region, or right-click to delete region or exit']);
            end
        end
    end
end
% end of segmenting cells

% close reference image
close(href)

% update / save data

UserData.status.changed = 1;
UserData.status.segmented  = 1;
enable_components(UserData,handles)

set(handles.figure1,'UserData',UserData);

set(handles.text1,'String',['saving data ...'])
save_ClickedCallback(hObject, eventdata, handles);
set(handles.text1,'String','Data saved')



% --------------------------------------------------------------------
function ProjectionMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ProjectionMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% do default z-projection: responds to the button on the toolbar
% can also specify the projection methods
% --------------------------------------------------------------------
function projection_ClickedCallback(hObject, eventdata, handles,method)
% hObject    handle to projection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');

%debug
disp('perform projection')
UserData.status
UserData.stack_range

if UserData.status.loaded % image loaded
    if nargin <4
        method = UserData.projection_method;
    end
    method
    
    I2 = zprojection(UserData.I(:,:,UserData.stack_range(1):UserData.stack_range(2)),method);
    size(I2)
    cla
    imagesc(I2)
    axis off;axis square;colormap gray;
    set(handles.text1,'String','try other projection/enhance method, or start to crop image')
    UserData.I2=I2;
    UserData.status.projected = 1;
    
    UserData.status.changed = 1;
    enable_components(UserData,handles)
    set(handles.figure1,'UserData',UserData);
    save_ClickedCallback(hObject, eventdata, handles);
    
end

% z-projection: standard deviation, respond to menu item:projection->std
% --------------------------------------------------------------------
function Proj_STD_Callback(hObject, eventdata, handles)
% hObject    handle to Proj_STD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projection_ClickedCallback(hObject, eventdata, handles,'std');

% z-projection: coefficient variation, responds to menu item:
% projection->cv
% --------------------------------------------------------------------
function Proj_CV_Callback(hObject, eventdata, handles)
% hObject    handle to Proj_CV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
projection_ClickedCallback(hObject, eventdata, handles,'cv');

% z-projection: max, rsponds to menu item:projection->max
% --------------------------------------------------------------------
function Proj_Max_Callback(hObject, eventdata, handles)
% hObject    handle to Proj_Max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projection_ClickedCallback(hObject, eventdata, handles,'max');

% z-projection, average, responds to menu item: projection->average
% --------------------------------------------------------------------
function Proj_Mean_Callback(hObject, eventdata, handles)
% hObject    handle to Proj_Mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

projection_ClickedCallback(hObject, eventdata, handles,'mean');

% set stack filter: responds to menu item: load_segmented_cells->stack filter
% --------------------------------------------------------------------
function setting_stack_filter_Callback(hObject, eventdata, handles)
% hObject    handle to ProjSetting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');

% prompt an input dialog
prompt = {'Only include stacks ','to '};
dlg_title = 'Stack filtering';
num_lines = 1;
def = {num2str(UserData.stack_range(1)),num2str(UserData.stack_range(2))};% default input

go_on = 1;
while go_on
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        return;
    end
    input1 =  str2num(answer{1});
    input2 =  str2num(answer{2});
    if input1 <1 | input2<1
        clicked=questdlg('Both numbers should be positive!','Oops!','Retry','Cancel','Don''t touch me!','Retry')
        switch clicked
            case 'Retry'
            otherwise
                go_on = 0;
        end
    elseif (input1 > input2)
        clicked=questdlg('The first number should not be larger than the second!','Oops!','Retry','Cancel','Don''t touch me!','Retry')
        switch clicked
            case 'Retry'
            otherwise
                go_on = 0;
        end
    elseif UserData.status.loaded
        max_stack = UserData.image_size(3);
        if input2 <= max_stack% the only right input
            UserData.stack_range=[input1 input2];
            set(handles.figure1,'UserData',UserData);
            go_on = 0;
        else
            clicked=questdlg('The second number is larger than the total number of stacks in current image!','Oops!','Retry','Cancel','Don''t touch me!','Retry')
            switch clicked
                case 'Retry'
                otherwise
                    go_on = 0;
            end
        end
    else
        go_on = 0;
    end
end

UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);

% use mouse wheel to explore stacks
% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');
if UserData.status.loaded
    UserData.current_stack =  UserData.current_stack + eventdata.VerticalScrollCount;
    if UserData.current_stack > UserData.image_size(3)
        UserData.current_stack = 1;
    elseif UserData.current_stack <1
        UserData.current_stack = UserData.image_size(3);
    end
    figure(handles.figure1)
    imagesc(UserData.I(:,:,UserData.current_stack));colormap gray;axis square;axis off;
    set(handles.text1,'String',['stack ' num2str(UserData.current_stack)]);
    drawnow;
end
set(handles.figure1,'UserData',UserData);

UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);


% --- Executes during object creation, after load_segmented_cells all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate axes1


% --- assign dots to cell
% ------------------------------------------------------------------------
function assign_dots_to_cells(handles)
UserData=get(handles.figure1,'UserData');
n_channel = length(UserData.channel_name);
n_nuclear_channel = length(UserData.nuclear_channel_name);

n_cell = length(UserData.cell);

set(handles.text1,'String','Assigning dots to cell ...');drawnow;
if size(UserData.dots,1)>1,
    UserData.dots(:,5)=NaN;
end

if ~isempty(UserData.dots),
    if ~isempty(UserData.dots),
        for i=1:n_cell % for each cell
            %if isfield(UserData.cell(1),'rim_edge'),
            %   edges=UserData.cell(i).rim_edge;
            if isfield(UserData.cell(i),'extended_edge_new') ...
                    & ~isnan(UserData.cell(i).extended_edge_new),
                edges=UserData.cell(i).extended_edge_new;
            else
                edges=UserData.cell(i).edge;
            end
            %in = inpoly(UserData.dots(:,1:2),UserData.cell(i).edge);
            in = inpoly(UserData.dots(:,1:2),edges);
            UserData.dots(in,5) = i;
            % if a dot is in multiple cells, it will be assigned to the one
            % with larger label
        end
    end
end


if ~isempty(UserData.nuclear_dots),
    if ~isempty(UserData.nuclear_dots),
        for i=1:n_cell % for each cell
            in = inpoly(UserData.nuclear_dots(:,1:2),UserData.cell(i).edge);
            UserData.nuclear_dots(in,5) = i;
            % if a dot is in multiple cells, it will be assigned to the one
            % with larger label
        end
    end
end

% statistics
UserData.num_dot_cell_channel = zeros(n_cell,n_channel);
UserData.num_nuclear_dot_cell_channel = zeros(n_cell,n_nuclear_channel);
for i=1:n_cell
    for j=1:n_channel
        if ~isempty(UserData.dots),
            UserData.num_dot_cell_channel(i,j) = sum(UserData.dots(:,4)==j & UserData.dots(:,5)==i);
        end
    end
    for j=1:n_nuclear_channel,
        if ~isempty(UserData.nuclear_dots),
            UserData.num_nuclear_dot_cell_channel(i,j) = sum(UserData.nuclear_dots(:,4)==j & UserData.nuclear_dots(:,5)==i);
        end
    end
end

% Finally for each cell calculate its summed intensity
cell_intensity=zeros(length(UserData.cell),1);


set(handles.text1,'String','Dots have been assigned to cells');drawnow;
UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);


% --------------------------------------------------------------------
function CountNormalization_Callback(hObject, eventdata, handles)
% hObject    handle to CountNormalization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function image_Callback(hObject, eventdata, handles)
% hObject    handle to image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MergeChannels_Callback(hObject, eventdata, handles)
% hObject    handle to MergeChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% go to previous stack
% --------------------------------------------------------------------
function previous_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData=get(handles.figure1,'UserData');
if UserData.status.loaded
    UserData.current_stack =  UserData.current_stack-1;
    if UserData.current_stack <1
        UserData.current_stack = UserData.image_size(3);
    end
    imagesc(UserData.I(:,:,UserData.current_stack));colormap gray;axis square;axis off;
    set(handles.text1,'String',['stack ' num2str(UserData.current_stack)]);
    drawnow;
end

UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);

% go to next stack
% --------------------------------------------------------------------
function next_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData=get(handles.figure1,'UserData');
if UserData.status.loaded
    UserData.current_stack =  UserData.current_stack+1;
    if UserData.current_stack > UserData.image_size(3)
        UserData.current_stack = 1;
    end
    imagesc(UserData.I(:,:,UserData.current_stack));colormap gray;axis square;axis off;
    set(handles.text1,'String',['stack ' num2str(UserData.current_stack)]);
    drawnow;
end

UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);





% let user specify the parameters for auto-thresholding
% --------------------------------------------------------------------
function setting_auto_thresholding_Callback(hObject, eventdata, handles)
% hObject    handle to setting_auto_thresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData=get(handles.figure1,'UserData');

prompt = {'Set the size of sliding window (odd integer, >=3)','Set the offset punishment for number of dots (positive number)'};
dlg_title = 'Set parameters for automatic thresholding';
num_lines = 1;
def = {num2str(UserData.auto_thresholding.width),num2str(UserData.auto_thresholding.offset)};
go_on = 1;
while go_on
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        return;
    end
    input1 =  floor(str2num(answer{1})/2)*2+1;
    input2 =  str2num(answer{2});
    if input1 <3 | input2<=0
        clicked=questdlg('Invalid input!','Oops!','Revise','Cancel','Don''t touch me!','Revise')
        switch clicked
            case 'Cancel'
                go_on = 0;
            otherwise
        end
    else
        UserData.auto_thresholding.width = input1;
        UserData.auto_thresholding.offset = input2;
        set(handles.figure1,'UserData',UserData);
        go_on = 0;
    end
end

UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);




% do image enhancement from toolbar
% --------------------------------------------------------------------
function image_enhance_ClickedCallback(hObject, eventdata, handles, method)
% hObject    handle to image_enhance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData=get(handles.figure1,'UserData');
if UserData.status.projected % image projected
    if nargin <4
        method = UserData.enhance_method;
    end
    UserData.I3 = UserData.I2;
    UserData.I2= enhance_image(UserData.I2,method);
    %imagesc(UserData.I3)
    imagesc(UserData.I2)
    axis off;axis square;colormap gray;
    set(handles.text1,'String','Choose another enhance/projection method, or start to crop image')
    UserData.status.enhanced = 1;
    set(handles.figure1,'UserData',UserData);
    enable_components(UserData,handles);
end
UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);

% do image enhancemetn from menu
% --------------------------------------------------------------------
function histeq_Callback(hObject, eventdata, handles)
% hObject    handle to histeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image_enhance_ClickedCallback(hObject, eventdata, handles, 'histeq')

% --------------------------------------------------------------------
function adapthisteq_Callback(hObject, eventdata, handles)
% hObject    handle to adapthisteq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image_enhance_ClickedCallback(hObject, eventdata, handles, 'adapthisteq')

% --------------------------------------------------------------------
function imadjust_Callback(hObject, eventdata, handles)
% hObject    handle to imadjust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image_enhance_ClickedCallback(hObject, eventdata, handles, 'imadjust')

% change default image enhance method
% --------------------------------------------------------------------
function setting_enhance_histeq_Callback(hObject, eventdata, handles)
% hObject    handle to setting_enhance_histeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData=get(handles.figure1,'UserData');
UserData.enhance_method = 'histeq';
hs = findobj('-regexp','Tag','setting_enhance');
for i=1:numel(hs)
    set(hs(i),'Checked','off');
end
set(hObject,'Checked','on');
UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);


% --------------------------------------------------------------------
function setting_enhance_adapthisteq_Callback(hObject, eventdata, handles)
% hObject    handle to setting_enhance_adapthisteq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData=get(handles.figure1,'UserData');
UserData.enhance_method = 'adapthisteq';
hs = findobj('-regexp','Tag','setting_enhance');
for i=1:numel(hs)
    set(hs(i),'Checked','off');
end
set(hObject,'Checked','on');
UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);

% --------------------------------------------------------------------
function setting_enhance_imadjust_Callback(hObject, eventdata, handles)
% hObject    handle to setting_enhance_imadjust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData=get(handles.figure1,'UserData');
UserData.enhance_method = 'imadjust';
hs = findobj('-regexp','Tag','setting_enhance');
for i=1:numel(hs)
    set(hs(i),'Checked','off');
end
set(hObject,'Checked','on');
UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);


% count from menu
% --------------------------------------------------------------------
function count_Callback(hObject, eventdata, handles)
% hObject    handle to count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
count_ClickedCallback(hObject, eventdata, handles)


% adjust brightness/contrast
% --------------------------------------------------------------------
function contrast_Callback(hObject, eventdata, handles)
% hObject    handle to contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');

h = figure(handles.figure1);

hh = imcontrast(h);

position1 = get(h, 'Position');
position2 = get(hh, 'Position');
set(hh, 'Position', [position1(1)+position1(3)+10 position1(2)+position1(4)-500 300 400]);

uiwait

% save data to workspace
% --------------------------------------------------------------------
function save_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to save_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData=get(handles.figure1,'UserData');
save
msgbox('Your data has been saved to current workspace. You can type ''load'' in the command line to check it.','Data saved','modal');




% --------------------------------------------------------------------
function reset_label_channel_peak_Callback(hObject, eventdata, handles)
% hObject    handle to reset_label_channel_peak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% reset image
% --------------------------------------------------------------------
function reset_image_Callback(hObject, eventdata, handles)
% hObject    handle to reset_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData=get(handles.figure1,'UserData');
imagesc(UserData.I(:,:,UserData.current_stack));colormap gray;axis square;axis off;
set(handles.text1,'String','you can now project, enhance or crop image')
set(handles.figure1,'UserData',UserData);


% to show dots or not
% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hs = findobj('-regexp','Tag','plot_cell')
hs = [hs; findobj('-regexp','Tag','plot_outline')]
hs = [hs ;findobj('-regexp','Tag','plot_dot_proj')];
hs = [hs ;findobj('-regexp','Tag','plot_neighbor')];
hs = [hs ;findobj('-regexp','Tag','plot_dot_bg')];
if ~isempty(hs)
    for i=1:length(hs)
        t=get(hs(i));
        if isfield(t,'Annotation')
            set(get(get(hs(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
        end
    end
end
if ~isempty(hs)
    for i=1:length(hs)
        t=get(hs(i));
        if isfield(t,'Annotation')
            set(get(get(hs(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
        end
    end
end
legend('off')
UserData=get(handles.figure1,'UserData');
% Hint: get(hObject,'Value') returns toggle state of checkbox3
if ~isempty(UserData.dots) | ~isempty(UserData.nuclear_dots)
    if (get(hObject,'Value') == get(hObject,'Max'))
        % Checkbox is checked-take approriate action
        UserData=get(handles.figure1,'UserData');
        [ndot,ndot_nuc]=plot_all_dots(UserData);
        if ndot>0 | ndot_nuc>0,
            if isfield(UserData,'full_channel_name')
                legend(UserData.full_channel_name,'FontSize',14);
            end
        end
    else
        % Checkbox is not checked-take approriate action
        hd = findobj('-regexp','Tag','plot_dot');
        delete(hd)
        set(handles.checkbox6,'Value',0);
    end
elseif (get(hObject,'Value') == get(hObject,'Max'))
    set(handles.text1,'String','No dots!!')
end


% to show cells or not
% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');
% Hint: get(hObject,'Value') returns toggle state of checkbox4
if ~isempty(UserData.cell)
    if (get(hObject,'Value') == get(hObject,'Max'))
        % Checkbox is checked-take approriate action
        UserData=get(handles.figure1,'UserData');
        plot_all_cell(UserData.cell);%plot_nuclear_outline(UserData)
%        plot_all_divide_lines(UserData);
    else
        % Checkbox is not checked-take approriate action
        hd = findobj('-regexp','Tag','plot_cell');
        delete(hd)
        hd = findobj('-regexp','Tag','plot_divide_lines');
        delete(hd)
        hd = findobj('-regexp','Tag','plot_nuclear_outline');
        delete(hd)
    end
elseif (get(hObject,'Value') == get(hObject,'Max'))
    set(handles.text1,'String','No cells!!')
end

% to show which channel, main or reference
% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
hs = findobj('-regexp','Tag','plot_cell')
hs = [hs; findobj('-regexp','Tag','plot_outline')]
hs = [hs ;findobj('-regexp','Tag','plot_dot_proj')];
hs = [hs ;findobj('-regexp','Tag','plot_neighbor')];
hs = [hs ;findobj('-regexp','Tag','plot_dot_bg')];
hs = [hs ;findobj('-regexp','Tag','plot_divide_lines')];
hs = [hs ;findobj('-regexp','Tag','plot_nuclear_outline')];
if ~isempty(hs)
    for i=1:length(hs)
        set(get(get(hs(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
    end
end
UserData=get(handles.figure1,'UserData');
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton2'
        % Code for when radiobutton1 is selected.
        if ~isempty(UserData.reference)
            imagesc(UserData.reference)
            axis off;axis square;colormap gray;
            plot_all_cell(UserData.cell);
        else
            set(handles.text1,'String','No reference images loaded!');drawnow
        end
    case 'radiobutton1'
        % Code for when radiobutton2 is selected.
        imagesc(UserData.I2)
        axis off;axis square;colormap gray;
        plot_all_cell(UserData.cell);
end

if ~isempty(UserData.cell)
    if (get(handles.checkbox4,'Value') == get(handles.checkbox4,'Max'))
        % Checkbox is checked-take approriate action
        plot_all_cell(UserData.cell);plot_nuclear_outline(UserData)
    else
        % Checkbox is not checked-take approriate action
        hd = findobj('-regexp','Tag','plot_cell');
        delete(hd)
        hd = findobj('-regexp','Tag','plot_divide_lines');
        delete(hd)
        hd = findobj('-regexp','Tag','plot_nuclear_outline');
        delete(hd)
    end
elseif (get(handles.checkbox4,'Value') == get(handles.checkbox4,'Max'))
    set(handles.text1,'String','No cells!!')
end

if ~isempty(UserData.dots)
    if (get(handles.checkbox3,'Value') == get(handles.checkbox3,'Max'))
        % Checkbox is checked-take approriate action
        if plot_all_dots(UserData)
            if isfield(UserData,'channel_name')
                legend(UserData.channel_name,'FontSize',14);
            end
        end
    else
        % Checkbox is not checked-take approriate action
        hd = findobj('-regexp','Tag','plot_dot');
        delete(hd)
    end
elseif (get(handles.checkbox3,'Value') == get(handles.checkbox3,'Max'))
    set(handles.text1,'String','No dots!!')
end


% --- Executes during object creation, after load_segmented_cells all properties.
function checkbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function plot_counts_Callback(hObject, eventdata, handles)
% hObject    handle to plot_counts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function count_dots_Callback(hObject, eventdata, handles)
% hObject    handle to count_dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
count_ClickedCallback(hObject, eventdata, handles)


% detect background dots
% --------------------------------------------------------------------
function detect_bg_dots_Callback(hObject, eventdata, handles)
% hObject    handle to detect_bg_dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');


b4=UserData.SIZE_THRESHOLD_FILTER;

cy5_channel=NaN;
for kk=1:length(UserData.channel_name)
    if ~isempty(findstr(lower(UserData.channel_name{kk}),'cy5'))
        cy5_channel=kk;
    end
end

suggested_size_threshold='100';


prompt = {'Remove dots larger than this size'};
dlg_title = 'Set size threshold for removal';
num_lines = 2;
%def = {'0','0','1',suggested_size_threshold};
def = {suggested_size_threshold};

go_on = 1;
while go_on
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        return;
    end    
    sz_threshold =  str2num(answer{1});    
    go_on=0;
end

set(handles.text1,'String','Finding potential background dots ...');drawnow;

quest_show=0;
if ~isfield(UserData,'fiji_xy')
    quest_show=1;    
    UserData.fiji_xy=[];
end
if isempty(UserData.fiji_xy)
    quest_show=1;
end
if quest_show
     questdlg('If using TS it is recommended to load fiji_xy first to prevent their removal','Warning','ok','ok');
end
[bgid] = find_background_dots(UserData.dots,UserData.dot_pixel_list,sz_threshold,UserData.fiji_xy);
bgds=UserData.dots(bgid,:);
hbg = plot(bgds(:,1),bgds(:,2),'yo','MarkerSize',10,'Tag','plot_dot_bg');
set(handles.text1,'String',{[num2str(numel(bgid)) ' background dots found (circled)'], ['at threshold =' num2str(sz_threshold)]});drawnow;
n_bg_dots = size(bgds,1);
n_dots = size(UserData.dots,1);
if n_bg_dots
    msg =[num2str(n_bg_dots) ' of ' num2str(n_dots) ' in total are bigger than threshold size.']
    choice = questdlg(msg,...
        'Remove potential background dots',...
        'Remove',...
        'Keep them',...
        'Remove');
    if ~strcmp(choice,'Keep them')
        UserData.dots(bgid,:)=[];
        UserData.SN(bgid,:)=[];
        UserData=refine_dots(UserData,bgid);
        delete(hbg)
        hd = findobj('Tag','plot_dot');
        delete(hd)
        plot_all_dots(UserData);
        set(handles.text1,'String','Potential background dots removed');drawnow;
    else
        delete(hbg)
        hd = findobj('Tag','plot_dot');
        delete(hd)
        plot_all_dots(UserData);
        set(handles.text1,'String','Potential background dots removed');drawnow;
    end
end

UserData.status.changed = 1;
% update component status
enable_components(UserData,handles);
% save data
set(handles.figure1,'UserData',UserData);
save_ClickedCallback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function uipanel1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function extract_TS_Callback(hObject, eventdata, handles)
% hObject    handle to extract_TS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');
if isfield(UserData,'dot_pixel_values')
    set(handles.text1,'String','Extracting TS data');pause(0.2)
    [UserData]=extract_TS_fiji(UserData);
    % output TS file
    fid=fopen([UserData.save_path UserData.save_name '_TS_results.txt'],'w');
    fprintf(fid,'Total #mRNA=%d\tTotal transcription rate=%.2f mRNA/hr\tDegradation rate (delta)=%.2f 1/hr\tmRNA lifetime=%.2f hr\tTotal # Pol2 on target=%.2f\n', ...
       UserData.Xst,UserData.beta,UserData.delta, log(2)/UserData.delta,UserData.Mtot);
    fprintf(fid,'TS x\tTS y\tTS z\tIntensity\t# Pol2\tTranscription rate([mRNA/hr]\n');
    for k=1:size(UserData.TS,1)
        fprintf(fid,'%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n',UserData.TS(k,1),UserData.TS(k,2),UserData.TS(k,3),UserData.TS(k,4),UserData.TS(k,5),UserData.TS(k,6));
    end
    fprintf(fid,'\n');
    % calculate the summed transcription rate per cell
    summed_TS_cell=zeros(length(UserData.cell),1);
    if ~isempty(UserData.TS)
        for k=1:length(UserData.cell)
            indin=find(inpoly(UserData.TS(:,1:2),UserData.cell(k).edge));
            if isempty(indin),
                summed_TS_cell(k)=0;
            else
                summed_TS_cell(k)=sum(UserData.TS(indin,end));
            end
        end
    end
    fprintf(fid,'Cell\t#mRNA\tTranscription rate [mRNA/hr]\n');
    for k=1:size(UserData.num_dot_cell_channel),
        fprintf(fid,'%d\t%d\t%.2f\n',k,UserData.num_dot_cell_channel(k,end),summed_TS_cell(k));
    end
    fclose(fid);
    set(handles.text1,'String','Done');    
end
set(handles.figure1,'UserData',UserData);
save_ClickedCallback(hObject, eventdata, handles);



% --------------------------------------------------------------------
function Analyze_Callback(hObject, eventdata, handles)
% hObject    handle to Analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function SET_LOG_SHARP_FACTOR_Callback(hObject, eventdata, handles)
% hObject    handle to SET_LOG_SHARP_FACTOR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');

prompt = {'LOG sharpning factor'};
dlg_title = 'Set factor to weigh LOG-filtered image in threshold selection';
num_lines = 1;
def = {num2str(UserData.SHARP_FACTOR)};
go_on = 1;
while go_on
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        return;
    end
    UserData.SHARP_FACTOR =  str2num(answer{1});
    
    go_on=0;
end

UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);


% --------------------------------------------------------------------
function USE_FAST_MODE_Callback(hObject, eventdata, handles)
% hObject    handle to USE_FAST_MODE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


UserData=get(handles.figure1,'UserData');

choice = questdlg('Fast Mode counting','Would you like to use fast counting mode?','Yes');
if strcmp(choice,'Yes'),
    UserData.FAST_MODE=1;
end
if strcmp(choice,'No'),
    UserData.FAST_MODE=0;
end

UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);


% --------------------------------------------------------------------
function Saturate_LOG_image_Callback(hObject, eventdata, handles)
% hObject    handle to Saturate_LOG_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UserData=get(handles.figure1,'UserData');

prompt = {'Saturate values in LOG image above Zscore X (-1 default means this is not done)'};
dlg_title = 'Set Zscore above which pixel values in the LOG image will be saturated (improves dynamic range of threshold selection)';
num_lines = 1;
def = {num2str(UserData.ZTHRESH)};
go_on = 1;
while go_on
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        return;
    end
    UserData.ZTHRESH =  str2num(answer{1});
    
    go_on=0;
end

UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);


% --------------------------------------------------------------------
function filter_image_outlier_Callback(hObject, eventdata, handles)
% hObject    handle to filter_image_outlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');

prompt = {'Filter outliers when displaying dot images to count?'};
dlg_title = 'Filter outliers when displaying dot images to count?(0,1)';
num_lines = 1;
def = {num2str(UserData.FILTER_OUTLIER)};
go_on = 1;
while go_on
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if isempty(answer)
        return;
    end
    UserData.FILTER_OUTLIER =  str2num(answer{1});
    
    go_on=0;
end

UserData.status.changed = 1;
enable_components(UserData,handles)
set(handles.figure1,'UserData',UserData);




% --------------------------------------------------------------------
function load_fiji_validated_Callback(hObject, eventdata, handles)
% hObject    handle to load_fiji_validated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');


[FijiName,PathName,FilterIndex] = uigetfile([UserData.save_path '*.txt'],'Choose Fiji validated TS coordinates'); 

UserData.fiji_xy=load([PathName FijiName]);

set(handles.figure1,'UserData',UserData);
save_ClickedCallback(hObject, eventdata, handles);

set(handles.figure1,'UserData',UserData);


% --------------------------------------------------------------------
function cmpute_probe_weight_Callback(hObject, eventdata, handles)
% hObject    handle to cmpute_probe_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');
set(handles.text1,'String','Computing probe weight factor, please wait (this may take a minute) ...');
[UserData.W,UserData.L,UserData.probe_spread_profile]=compute_correction_factor();
set(handles.figure1,'UserData',UserData);
save_ClickedCallback(hObject, eventdata, handles);
set(handles.text1,'String','Done');
set(handles.figure1,'UserData',UserData);




% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Output_dot_counts_Callback(hObject, eventdata, handles)
% hObject    handle to Output_dot_counts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

UserData=get(handles.figure1,'UserData');

fid=fopen([UserData.save_path UserData.save_name '_dot_counts.txt'],'w');
fprintf(fid,'Cell\tCell Area (pixel squared)\t# stacks analyzed\t#mRNA\n');
num_stacks=diff(UserData.stack_range);
for k=1:size(UserData.num_dot_cell_channel,1),
    fprintf(fid,'%d\t%.0f\t%d\t%d\n',k,UserData.cell(k).area,num_stacks,UserData.num_dot_cell_channel(k,end));
end
fprintf(fid,'\nAll dots\n');
fprintf(fid,'X\tY\tZ\tCell\n');
for i=1:size(UserData.dots)
    fprintf(fid,'%.1f\t%.1f\t%.1f\t%d\n',UserData.dots(i,1),...
        UserData.dots(i,2),UserData.dots(i,3),UserData.dots(i,end));
end
fclose(fid);
set(handles.text1,'String','Done');
set(handles.figure1,'UserData',UserData);
save_ClickedCallback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function load_segmented_cells_Callback(hObject, eventdata, handles)
% hObject    handle to load_segmented_cells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


UserData=get(handles.figure1,'UserData');
UserData_temp=UserData;

[name, path, filter_index] = uigetfile('*.ima','Open ima file with segmented cells to transfer');
if name==0
    return;
end

% load file
load([path name],'-mat');

UserData_temp.cell=UserData.cell;
UserData=UserData_temp;

h = findobj('-regexp','Tag','plot_cell');
delete(h)
h = findobj('-regexp','Tag','plot_nuclear_outline');
delete(h)
plot_all_cell(UserData.cell);
assign_dots_to_cells(handles);

set(handles.figure1,'UserData',UserData);
