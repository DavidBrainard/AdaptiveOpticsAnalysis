function varargout = AutoAOMontagingGUI(varargin)
% AUTOAOMONTAGINGGUI MATLAB code for AutoAOMontagingGUI.fig
%      AUTOAOMONTAGINGGUI, by itself, creates a new AUTOAOMONTAGINGGUI or raises the existing
%      singleton*.
%
%      H = AUTOAOMONTAGINGGUI returns the handle to a new AUTOAOMONTAGINGGUI or the handle to
%      the existing singleton*.
%
%      AUTOAOMONTAGINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AUTOAOMONTAGINGGUI.M with the given input arguments.
%
%      AUTOAOMONTAGINGGUI('Property','Value',...) creates a new AUTOAOMONTAGINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AutoAOMontagingGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AutoAOMontagingGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AutoAOMontagingGUI

% Last Modified by GUIDE v2.5 15-Dec-2015 14:47:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AutoAOMontagingGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AutoAOMontagingGUI_OutputFcn, ...
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


% --- Executes just before AutoAOMontagingGUI is made visible.
function AutoAOMontagingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AutoAOMontagingGUI (see VARARGIN)

% Choose default command line output for AutoAOMontagingGUI
handles.output = hObject;

%exim= imread('D:\MCW\Normative_Metrics_2014\data\Normal Data\AD_10253\AD_10253_OS_seeded_0.5176umppx_trim_flat.tif');
%exim= imread('C:\Users\Min\Dropbox (Aguirre-Brainard Lab)\AOSLOImageProcessing\ConstructMontage\InputImages_Set1\CO_20140408_NIR_OS_0012_ref_126_lps_8_lbss_8_sr_n_50_cropped_5.tif');
% addpath(genpath('./private'))
% exim= imread('C:\Users\Min\Documents\Research\AdaptiveOpticsMosaic\NewDataSet2_DF_2015_8_26\CS_13213_20150302_OS_1p00_Montage_DarkField.bmp');
% axes(handles.canvas);
% imagesc(exim); colormap gray; axis off;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AutoAOMontagingGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AutoAOMontagingGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function filemenu_Callback(hObject, eventdata, handles)
% hObject    handle to filemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in imageList.
function imageList_Callback(hObject, eventdata, handles)
% hObject    handle to imageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns imageList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imageList

index_selected = get(handles.imageList,'Value');
axes(handles.canvas);
img = imread(fullfile(handles.imgfolder_name,handles.imageFile_names{index_selected}));
imagesc(img); colormap gray; axis equal; axis off;


% --- Executes during object creation, after setting all properties.
function imageList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in imagefolder.
function imagefolder_Callback(hObject, eventdata, handles)
% hObject    handle to imagefolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.imgfolder_name = uigetdir;
set(handles.selectFolderText, 'String', handles.imgfolder_name) ;
Allfiles = dir(strcat(handles.imgfolder_name,'\*.tif'));
Allfiles = {Allfiles.name};
confocal = sort(Allfiles(~cellfun(@isempty, strfind(Allfiles, '_confocal_'))));

handles.imageFile_names = confocal;
set(handles.imageList,'String',handles.imageFile_names,...
	'Value',1)
guidata(hObject, handles);


% --- Executes on button press in selectPosFile.
function selectPosFile_Callback(hObject, eventdata, handles)
% hObject    handle to selectPosFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile(fullfile(pwd,'*.xlsx'));
handles.postionFile_name = fullfile(PathName,FileName);
set(handles.posFileText, 'String', handles.postionFile_name) ;
guidata(hObject, handles);

% --- Executes on selection change in montageList.
function montageList_Callback(hObject, eventdata, handles)
% hObject    handle to montageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns montageList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from montageList

index_selected = get(handles.montageList,'Value');
axes(handles.canvas);
img = imread(fullfile(handles.outputFolder_name,handles.combinedFile_names{index_selected}));
imagesc(img); colormap gray; axis equal; axis off;


% --- Executes during object creation, after setting all properties.
function montageList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to montageList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in montageAll.
function montageAll_Callback(hObject, eventdata, handles)
% hObject    handle to montageAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tic
handles.combinedFile_names = AOMosiacAllMultiModal(handles.imgfolder_name,handles.postionFile_name,'A1:C59',handles.outputFolder_name,[1 2 3], 0 );
toc

set(handles.montageList,'String',handles.combinedFile_names,...
	'Value',1)
img=imread(fullfile(handles.outputFolder_name,handles.combinedFile_names{1}));
axes(handles.canvas);
imagesc(img); colormap gray; axis equal; axis off;
guidata(hObject, handles);



% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 close all;


% --- Executes on button press in outputFolder.
function outputFolder_Callback(hObject, eventdata, handles)
% hObject    handle to outputFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.outputFolder_name = uigetdir;
set(handles.outputFolderText, 'String', handles.outputFolder_name) ;
guidata(hObject, handles);
