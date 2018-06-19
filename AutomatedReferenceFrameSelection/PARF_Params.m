function varargout = PARF_Params(varargin)
% PARF_PARAMS MATLAB code for PARF_Params.fig
%      PARF_PARAMS, by itself, creates a new PARF_PARAMS or raises the existing
%      singleton*.
%
%      H = PARF_PARAMS returns the handle to a new PARF_PARAMS or the handle to
%      the existing singleton*.
%
%      PARF_PARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARF_PARAMS.M with the given input arguments.
%
%      PARF_PARAMS('Property','Value',...) creates a new PARF_PARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PARF_Params_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PARF_Params_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PARF_Params

% Last Modified by GUIDE v2.5 19-Jun-2018 09:59:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PARF_Params_OpeningFcn, ...
                   'gui_OutputFcn',  @PARF_Params_OutputFcn, ...
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


% --- Executes just before PARF_Params is made visible.
function PARF_Params_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PARF_Params (see VARARGIN)

% Choose default command line output for PARF_Params
handles.output = hObject;

handles.NUM_REF_OUTPUT = str2double( get(handles.num_ref_output,'String') );
handles.STRIP_SIZE = str2double( get(handles.parf_strip_size,'String') );
handles.BAD_STRIP_THRESHOLD = str2double( get(handles.bad_strip_thresh,'String') );
handles.MIN_NUM_FRAMES_PER_GROUP = str2double( get(handles.min_frames_per_group,'String') );

tabledata = get(handles.modality_table,'Data');
handles.MODALITIES = tabledata(:,1)';
handles.MODALITY_WEIGHTS = cell2mat(tabledata(:,2))';

handles.LPS = str2double( get(handles.lps,'String') );
handles.LBSS = str2double( get(handles.lbss,'String') );
handles.OVERLAP = str2double( get(handles.overlap,'String') );
handles.NUM_FRAMES = str2double( get(handles.num_frames,'String') );
handles.THRESHOLD = str2double( get(handles.threshold,'String') );
handles.OUTPUT_TIFS = get(handles.output_tifs,'Value')==1;
handles.OUTPUT_AVIS = get(handles.output_vids,'Value')==1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PARF_Params wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PARF_Params_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



output = struct();
if ~handles.cancelled
    output.NUM_REF_OUTPUT = handles.NUM_REF_OUTPUT;
    output.STRIP_SIZE = handles.STRIP_SIZE;
    output.BAD_STRIP_THRESHOLD = handles.BAD_STRIP_THRESHOLD;
    output.MIN_NUM_FRAMES_PER_GROUP = handles.MIN_NUM_FRAMES_PER_GROUP;

    tabledata = get(handles.modality_table,'Data');
    output.MODALITIES = tabledata(:,1)';
    output.MODALITY_WEIGHTS = cell2mat(tabledata(:,2))';

    output.LPS = handles.LPS;
    output.LBSS = handles.LBSS;
    output.OVERLAP = handles.OVERLAP;
    output.NUM_FRAMES = handles.NUM_FRAMES;
    output.THRESHOLD = handles.THRESHOLD;
    output.OUTPUT_TIFS = handles.OUTPUT_TIFS;
    output.OUTPUT_AVIS = handles.OUTPUT_AVIS;
end

% Get default command line output from handles structure
varargout{1} = output;
close(handles.figure1);



function num_ref_output_Callback(hObject, eventdata, handles)
% hObject    handle to num_ref_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

value = str2double(get(hObject,'String'));

if isnan(value)    
    set(hObject,'String', num2str(handles.NUM_REF_OUTPUT))
else
    handles.NUM_REF_OUTPUT = value;
end
% Hints: get(hObject,'String') returns contents of num_ref_output as text
%        str2double(get(hObject,'String')) returns contents of num_ref_output as a double


% --- Executes during object creation, after setting all properties.
function num_ref_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_ref_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function parf_strip_size_Callback(hObject, eventdata, handles)
% hObject    handle to parf_strip_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject,'String'));

if isnan(value)    
    set(hObject,'String', num2str(handles.STRIP_SIZE))
else
    handles.STRIP_SIZE = value;
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of parf_strip_size as text
%        str2double(get(hObject,'String')) returns contents of parf_strip_size as a double


% --- Executes during object creation, after setting all properties.
function parf_strip_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parf_strip_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bad_strip_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to bad_strip_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject,'String'));

if isnan(value)    
    set(hObject,'String', num2str(handles.BAD_STRIP_THRESHOLD))
else
    handles.BAD_STRIP_THRESHOLD = value;
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of bad_strip_thresh as text
%        str2double(get(hObject,'String')) returns contents of bad_strip_thresh as a double


% --- Executes during object creation, after setting all properties.
function bad_strip_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bad_strip_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min_frames_per_group_Callback(hObject, eventdata, handles)
% hObject    handle to min_frames_per_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject,'String'));

if isnan(value)    
    set(hObject,'String', num2str(handles.MIN_NUM_FRAMES_PER_GROUP))
else
    handles.MIN_NUM_FRAMES_PER_GROUP = value;
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of min_frames_per_group as text
%        str2double(get(hObject,'String')) returns contents of min_frames_per_group as a double


% --- Executes during object creation, after setting all properties.
function min_frames_per_group_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_frames_per_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lbss_Callback(hObject, eventdata, handles)
% hObject    handle to lbss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject,'String'));

if isnan(value)    
    set(hObject,'String', num2str(handles.LBSS))
else
    handles.LBSS = value;
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of lbss as text
%        str2double(get(hObject,'String')) returns contents of lbss as a double


% --- Executes during object creation, after setting all properties.
function lbss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lps_Callback(hObject, eventdata, handles)
% hObject    handle to lps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject,'String'));

if isnan(value)    
    set(hObject,'String', num2str(handles.LPS))
else
    handles.LPS = value;
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of lps as text
%        str2double(get(hObject,'String')) returns contents of lps as a double


% --- Executes during object creation, after setting all properties.
function lps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function overlap_Callback(hObject, eventdata, handles)
% hObject    handle to overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject,'String'));

if isnan(value)    
    set(hObject,'String', num2str(handles.OVERLAP))
else
    handles.OVERLAP = value;
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of overlap as text
%        str2double(get(hObject,'String')) returns contents of overlap as a double


% --- Executes during object creation, after setting all properties.
function overlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_frames_Callback(hObject, eventdata, handles)
% hObject    handle to num_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject,'String'));

if isnan(value)    
    set(hObject,'String', num2str(handles.NUM_FRAMES))
else
    handles.NUM_FRAMES = value;
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of num_frames as text
%        str2double(get(hObject,'String')) returns contents of num_frames as a double


% --- Executes during object creation, after setting all properties.
function num_frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject,'String'));

if isnan(value)    
    set(hObject,'String', num2str(handles.THRESHOLD))
else
    handles.THRESHOLD = value;
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double


% --- Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in output_vids.
function output_vids_Callback(hObject, eventdata, handles)
% hObject    handle to output_vids (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');

if isnan(value)    
    set(hObject,'Value', handles.OUTPUT_AVIS)
else
    handles.OUTPUT_AVIS = value;
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of output_vids


% --- Executes on button press in output_tifs.
function output_tifs_Callback(hObject, eventdata, handles)
% hObject    handle to output_tifs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');

if isnan(value)    
    set(hObject,'Value', handles.OUTPUT_TIFS)
else
    handles.OUTPUT_TIFS = value;
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of output_tifs


% --- Executes on button press in beginbutton.
function beginbutton_Callback(hObject, eventdata, handles)
% hObject    handle to beginbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cancelled = false;
guidata(hObject, handles);

uiresume(handles.figure1);


% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cancelled = true;
guidata(hObject, handles);

uiresume(handles.figure1);

