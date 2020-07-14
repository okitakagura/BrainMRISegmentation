function varargout = MainGUI(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MainGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MainGUI_OutputFcn, ...
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


% --- Executes just before MainGUI is made visible.
function MainGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MainGUI (see VARARGIN)

% Choose default command line output for MainGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MainGUI wait for user response (see UIRESUME)
% uiwait(handles.MainFig);

% ����ʱ����ʾ
set(handles.timestr,'string',datestr(now,0));
htimer = timer('StartDelay',1,'TimerFcn',...
    'htimestr=findall(0,''tag'',''timestr'');set(htimestr,''string'',datestr(now,0));',...
    'Period',1,'ExecutionMode','fixedSpacing','tag','showtime');
start(htimer);



% --- Outputs from this function are returned to the command line.
function varargout = MainGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function hslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes during object creation, after setting all properties.
function hedit_fname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hedit_fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function hedit_heb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hedit_heb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function hedit_xuhao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hedit_xuhao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function hedit_zhishu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hedit_zhishu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function hedit_zuixiao_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hedit_zuixiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function hedit_cishu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hedit_cishu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function hedit_leibie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hedit_leibie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function hedit_alfa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hedit_alfa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function hmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function hedit_info_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hedit_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function hslider_Callback(hObject, eventdata, handles)
% hObject    handle to hslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider




% --- Executes on button press in btnOpen.
function btnOpen_Callback(hObject, eventdata, handles)
% hObject    handle to btnOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.mat', '��ͼ�������ļ�');
if filename~=0
    set(gcf,'Pointer','watch'); % �ı����״̬Ϊ�ȴ�״̬
    drawnow;
    set(handles.hedit_fname,'string',filename);
   % fid=fopen([pathname filename],'r');
   % F=fread(fid,'uint8');
   % fclose(fid);
    load([pathname filename]);
    F = brainweb;
    set(hObject,'UserData',F);
   % set(handles.hedit_xuhao,'string',num2str(3));
    slicenum = get(handles.hedit_xuhao,'UserData');    %���
    I=F{slicenum};
%     I=F(1+181*181*(3-1):181*181*3);
%     I=reshape(I,181,181);
%     I=imrotate(I,90);
    axes(handles.hyuansiaxes);
    himage = imshow(I,[]);
    set(gcf,'NextPlot','add'); % This is to cinquer the bug of imshow
    % in image toolbox ver 5.0.0 and 5.0.1.
    set(handles.hyuansiaxes,'UserData',I);
    set(gcf,'Pointer','arrow');
    set(handles.hslider,'enable','off'); % ������������
    set(handles.btnsave,'enable','off');    % ������水ť������ 
else
    set(handles.hedit_fname,'string','��δ��');
end


function hedit_heb_Callback(hObject, eventdata, handles)
% hObject    handle to hedit_heb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hedit_heb as text
%        str2double(get(hObject,'String')) returns contents of hedit_heb as a double
str = get(hObject,'string');
data = str2num(str);
if isempty(data)  % ���������Ч�Լ��
    errordlg('�������Ϊ��ֵ��','��������');
    set(hObject,'BackgroundColor','r');
else
    set(hObject,'BackgroundColor','w');
    set(hObject,'UserData',data);
end


function hedit_xuhao_Callback(hObject, eventdata, handles)
% hObject    handle to hedit_xuhao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hedit_xuhao as text
%        str2double(get(hObject,'String')) returns contents of hedit_xuhao as a double
str = get(hObject,'string');
data = str2num(str);
if isempty(data)  % ���������Ч�Լ��
    errordlg('�������Ϊ��ֵ��','��������');
    set(hObject,'BackgroundColor','r');
elseif data<1  || data~=round(data) % �����ų�����Χ����������
    errordlg('�������Ϊ��������','��������');
    set(hObject,'BackgroundColor','r');
else
    set(hObject,'BackgroundColor','w');
    % �������������ʾԭͼ
    F = get(handles.btnOpen,'UserData');
    if isempty(F)
        msgbox('��δ��ͼ�������ļ�','��������');
        return;
    else
        I=F{data};
%         I=F(1+181*181*(data-1):181*181*data);
%         I=reshape(I,181,181);
%         I=imrotate(I,90);
        axes(handles.hyuansiaxes);
        himage = imshow(I,[]);
        set(gcf,'NextPlot','add'); % This is to cinquer the bug of imshow 
                        % in image toolbox ver 5.0.0 and 5.0.1. 
        set(handles.hyuansiaxes,'UserData',I);
        set(handles.hedit_xuhao,'UserData',data);
     end
end

function hedit_zhishu_Callback(hObject, eventdata, handles)
% hObject    handle to hedit_zhishu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hedit_zhishu as text
%        str2double(get(hObject,'String')) returns contents of hedit_zhishu as a double
str = get(hObject,'string');
data = str2num(str);
if isempty(data)  % ���������Ч�Լ��
    errordlg('�������Ϊ��ֵ��','��������');
    set(hObject,'BackgroundColor','r');
elseif data<=1 % ���ָ��������1
    errordlg('����������1','��������');
    set(hObject,'BackgroundColor','r');
else    
    set(hObject,'BackgroundColor','w');
    set(hObject,'UserData',data);
end



function hedit_zuixiao_Callback(hObject, eventdata, handles)
% hObject    handle to hedit_zuixiao (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hedit_zuixiao as text
%        str2double(get(hObject,'String')) returns contents of hedit_zuixiao as a double
str = get(hObject,'string');
data = str2num(str);
if isempty(data)  % ���������Ч�Լ��
    errordlg('�������Ϊ��ֵ��','��������');
    set(hObject,'BackgroundColor','r');
elseif data<0 % ���ָ��������1
    errordlg('����������0','��������');
    set(hObject,'BackgroundColor','r');
else    
    set(hObject,'BackgroundColor','w');
    set(hObject,'UserData',data);
end



function hedit_cishu_Callback(hObject, eventdata, handles)
% hObject    handle to hedit_cishu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hedit_cishu as text
%        str2double(get(hObject,'String')) returns contents of hedit_cishu as a double
str = get(hObject,'string');
data = str2num(str);
if isempty(data)  % ���������Ч�Լ��
    errordlg('�������Ϊ��ֵ��','��������');
    set(hObject,'BackgroundColor','r');
elseif data<1 || data~=round(data) % �����ų�����Χ����������
    errordlg('�������Ϊ���ڵ���1��������','��������');
    set(hObject,'BackgroundColor','r');
else    
    set(hObject,'BackgroundColor','w');
    set(hObject,'UserData',data);
    
end


function hedit_leibie_Callback(hObject, eventdata, handles)
% hObject    handle to hedit_leibie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hedit_leibie as text
%        str2double(get(hObject,'String')) returns contents of hedit_leibie as a double
str = get(hObject,'string');
data = str2num(str);
if isempty(data)  % ���������Ч�Լ��
    errordlg('�������Ϊ��ֵ��','��������');
    set(hObject,'BackgroundColor','r');
elseif data<2 || data~=round(data) % �����ų�����Χ����������
    errordlg('�������Ϊ���ڵ���2��������','��������');
    set(hObject,'BackgroundColor','r');
else    
    set(hObject,'BackgroundColor','w');
    set(hObject,'UserData',data);    
end


function hedit_alfa_Callback(hObject, eventdata, handles)
% hObject    handle to hedit_alfa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of hedit_alfa as text
%        str2double(get(hObject,'String')) returns contents of hedit_alfa as a double
% --- Executes during object creation, after setting all properties.
% --- Executes on button press in btnrun.
str = get(hObject,'string');
data = str2num(str);
if isempty(data)  % ���������Ч�Լ��
    errordlg('�������Ϊ��ֵ��','��������');
    set(hObject,'BackgroundColor','r');
elseif data<0 % ���ָ��������1
    errordlg('����������0','��������');
    set(hObject,'BackgroundColor','r');
else    
    set(hObject,'BackgroundColor','w');
    set(hObject,'UserData',data);
end


function btnrun_Callback(hObject, eventdata, handles)
% hObject    handle to btnrun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ���Ľ���ؼ�״̬
set(handles.hedit_info,'enable','on');  % ��Ϣ����༭�����
set(handles.btnsave,'enable','on');    % ������水ť����

expo = get(handles.hedit_zhishu,'UserData');    % ָ��m
max_iter = get(handles.hedit_cishu,'UserData'); % ����������
min_impro = get(handles.hedit_zuixiao,'UserData');  % Ŀ�꺯����С�ı���
cluster_n = get(handles.hedit_leibie,'UserData');   % ������Ŀ
alfa=get(handles.hedit_alfa,'UserData');    % ����Ȩ��
display = get(handles.chkbx,'UserData');    % ÿ�ε����Ƿ������Ϣ
algorithm = get(handles.hmenu,'UserData');  % ѡ����㷨����
data_load = get(handles.hyuansiaxes,'UserData'); % ��ȡԭʼͼ������

% %��ȡ��׼�ָ�ͼ��
slicenum = get(handles.hedit_xuhao,'UserData');    %���
mark = Mark('brainweb/phantom_1.0mm_normal_crisp.rawb',slicenum);
read = readrawb('brainweb/t1_icbm_normal_1mm_pn0_rf0.rawb',slicenum);
[n1,n2] = size(read);
[row,col] = size(read);
read_new = read;
for i = 1:row
    for j = 1:col
        if mark(i,j) == 0
            read_new(i,j)=0;
        end
    end
end
read_new=imrotate(read_new, 90); 
real_label=imrotate(mark, 90);
real_count1 = 0;
real_count2 = 0;
real_count3 = 0;
real_count4 = 0;
for x=1:n2
    for y=1:n1
        if real_label(x,y) == 0
            real_count1 = real_count1 + 1;
        elseif real_label(x,y) == 1
            real_count2 = real_count2 + 1;
        elseif real_label(x,y) == 2
            real_count3 = real_count3 + 1;
        elseif real_label(x,y) == 3
            real_count4 = real_count4 + 1;
        end
    end
end
% data_load = read_new;


data2 = data_load;
[a b]=size(data2);%ͼ�����Ĵ�С
data=data2(:);%ת��Ϊ������
data=double(data);
data_n = size(data, 1); % ���data�ĵ�һά(rows)��,����������
in_n = size(data, 2);   % ���data�ĵڶ�ά(columns)����������ֵ����

%obj_fcn = zeros(max_iter, 1);	% ��ʼ���������obj_fcn
if algorithm==1 % FCMͼ��ָ�
    [U, center,obj_fcn] = FCMClust2(data_load, cluster_n, expo, max_iter, min_impro, display);   
elseif algorithm==2% FCM_Sͼ��ָ�
       [U, center,obj_fcn] = FCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);
elseif algorithm==3% KFCMͼ��ָ�
       [U, center,obj_fcn] = KFCMClust2brainweb(data_load, cluster_n, expo, max_iter, min_impro, display);
elseif algorithm==4% KFCM_Sͼ��ָ�
       [U, center,obj_fcn] = KFCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);
elseif algorithm==5% MyKFCM_Sͼ��ָ�
       [U, center,obj_fcn] = MyKFCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);
end

fm = U.^expo;
V_pc = sum(fm(:))/data_n;
V_xieben = XieBeniInverted(U, expo, center, data);

% �ָ�����ʾ
data=data';
wholeG=zeros(size(data));
maxU=max(U);
x=[1,2,3,4];
[center,id] = sort(center);
x=x(id);
for k=1:cluster_n
    indexk=(U(k,:)==maxU);
    count{k} = sum(indexk(:));
    Ik = indexk.*data;
    Ik = reshape(Ik,a,b);
    result{k} = Ik;%result{k}��¼�ڣ����е����ؼ���λ����Ϣ
    wholeG(indexk) = x(k);%wholeG��¼�ָ������ͼ�����Ϣ
end
wholeG=reshape(wholeG,a,b);%������ת��Ϊa*b��С�ľ���
[label_new,accuracy]=succeed(real_label,cluster_n,wholeG);
wholeG = reshape(label_new,a,b);

set(handles.btnrun,'UserData',result);

% ��ʾ��һ��
axes(handles.hzhengtiaxes);
imshow(result{x(2)},[]);
set(gcf,'NextPlot','add');

% ��ʾ�ڶ���
axes(handles.h1leiaxes);
imshow(result{x(3)},[]);
set(gcf,'NextPlot','add');

% ��ʾ������
axes(handles.h2leiaxes);
imshow(result{x(4)},[]);
set(gcf,'NextPlot','add');

set(handles.hslider,'enable','on'); % ����������

infostr =sprintf('accuracy =  %f\n\t V_pc = %f\n\t V_xieben = %f\n',accuracy,V_pc,V_xieben);
set(handles.infortext,'String',infostr);



% --- Executes on button press in btnhelp.
function btnhelp_Callback(hObject, eventdata, handles)
% hObject    handle to btnhelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.hedit_info,'enable','on');  % ��Ϣ����༭�����
set(handles.btnsave,'enable','on');    % ������水ť����
infostr =" ";
set(handles.infortext,'String',infostr);
expo = get(handles.hedit_zhishu,'UserData');    % ָ��m
max_iter = get(handles.hedit_cishu,'UserData'); % ����������
min_impro = get(handles.hedit_zuixiao,'UserData');  % Ŀ�꺯����С�ı���
cluster_n = get(handles.hedit_leibie,'UserData');   % ������Ŀ
alfa=get(handles.hedit_alfa,'UserData');    % ����Ȩ��
display = get(handles.chkbx,'UserData');    % ÿ�ε����Ƿ������Ϣ
algorithm = get(handles.hmenu,'UserData');  % ѡ����㷨����
F = get(handles.btnOpen,'UserData');
v_p1 = 0; 
v_p2 = 0; 
v_p3 = 0; 
v_p4 = 0;
n_start = 10;
n_end =154;
for slicenum = n_start:n_end
    data_load = F{slicenum};
    data2 = data_load;
[a b]=size(data2);%ͼ�����Ĵ�С
data=data2(:);%ת��Ϊ������
data=double(data);
data_n = size(data, 1); % ���data�ĵ�һά(rows)��,����������
in_n = size(data, 2);   % ���data�ĵڶ�ά(columns)����������ֵ����

if algorithm==1 % FCMͼ��ָ�
    [U, center,obj_fcn] = FCMClust2(data_load, cluster_n, expo, max_iter, min_impro, display);   
elseif algorithm==2% FCM_Sͼ��ָ�
       [U, center,obj_fcn] = FCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);
elseif algorithm==3% KFCMͼ��ָ�
       [U, center,obj_fcn] = KFCMClust2brainweb(data_load, cluster_n, expo, max_iter, min_impro, display);
elseif algorithm==4% KFCM_Sͼ��ָ�
       [U, center,obj_fcn] = KFCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);
elseif algorithm==5% MyKFCM_Sͼ��ָ�
       [U, center,obj_fcn] = MyKFCM_S(data_load, cluster_n, expo, max_iter, min_impro, display,alfa);
end
%ͼ��ָ�
data=data';
wholeG=zeros(size(data));
maxU=max(U);
x=[1,2,3,4];
[center,id] = sort(center);
x=x(id);
for k=1:cluster_n
    indexk=(U(k,:)==maxU);
    count{k} = sum(indexk(:));
    Ik = indexk.*data;
    Ik = reshape(Ik,a,b);
    result{k} = Ik;%result{k}��¼�ڣ����е����ؼ���λ����Ϣ
    wholeG(indexk) = x(k);%wholeG��¼�ָ������ͼ�����Ϣ
end
wholeG=reshape(wholeG,a,b);%������ת��Ϊa*b��С�ľ���

% ��ʾԭͼ��
axes(handles.hyuansiaxes);
imshow(data_load,[]);
set(gcf,'NextPlot','add');

% ��ʾ��һ��
axes(handles.hzhengtiaxes);
imshow(result{x(2)},[]);
set(gcf,'NextPlot','add');

% ��ʾ�ڶ���
axes(handles.h1leiaxes);
imshow(result{x(3)},[]);
set(gcf,'NextPlot','add');

% ��ʾ������
axes(handles.h2leiaxes);
imshow(result{x(4)},[]);
set(gcf,'NextPlot','add');
s_p2 = count{x(2)} * 0.01;
s_p3 = count{x(3)} * 0.01;
s_p4 = count{x(4)} * 0.01;
if slicenum == n_start || slicenum == n_end
    v_p2 = v_p2+s_p2;
    v_p3 = v_p3+s_p3;
    v_p4 = v_p4+s_p4;
else
    v_p2 = v_p2+2*s_p2;
    v_p3 = v_p3+2*s_p3;
    v_p4 = v_p4+2*s_p4;
end
end
V_p2 = v_p2*0.1/2;
V_p3 = v_p3*0.1/2;
V_p4 = v_p4*0.1/2;

infostr =sprintf('�Լ�Һ��� = %f\n\t ������� = %f\n\t ������� = %f\n',V_p2,V_p3,V_p4);
set(handles.infortext,'String',infostr);


% --- Executes on button press in btnsave.
function btnsave_Callback(hObject, eventdata, handles)
% hObject    handle to btnsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

result = get(handles.btnrun,'UserData');

ButtonName=questdlg('������Ϊ������ʽ?', ...
                       '�ָ�������', ...
                       'JPGͼƬ','MAT�ļ�','MAT�ļ�');
switch ButtonName,
     case 'MAT�ļ�',
         [fname,path] = uiputfile('*.mat','����Ϊ');
         if path==0  % ȡ���ļ��������
             return;
         end
         save([path fname],'result','-v6');         
     case 'JPGͼƬ',
         directoryname = uigetdir;
         if directoryname == 0 % canceled
             return;
         end
         pwdir = pwd; % ��õ�ǰĿ¼
         cd(directoryname); % ת��ѡ��Ŀ¼
         fnametmp = datestr(now,31);
         fnametmp = strrep(fnametmp,':','-');   % �ļ����в��ܹ���:��
         fnametmp = strrep(fnametmp,' ','_');
         for k = 1:length(result)-1
             fnamestr = [fnametmp 'Kind',num2str(k),'.jpg'];
             fcell ={fnamestr};
             tmpM = result{k};
             imwrite(tmpM,fcell{1},'jpg');
         end
         fnamestr = [fnametmp,'Whole.jpg'];
         fcell ={fnamestr};
         tmpM = result{end};
         imwrite(tmpM,fcell{1},'jpg');
         cd(pwdir);
end % switch                   
 



% --- Executes on button press in btnquit.
function btnquit_Callback(hObject, eventdata, handles)
% hObject    handle to btnquit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% selection = questdlg(['�˳� ' get(handles.MainFig,'Name') '?'],...
%                      ['�˳� ...'],...
%                      '��','��','��');
% if strcmp(selection,'��')
%     return;
% end
htimer = timerfind('tag','showtime');
stop(htimer);
delete(htimer);
delete(handles.MainFig);



% --- Executes on selection change in hmenu.
function hmenu_Callback(hObject, eventdata, handles)
% hObject    handle to hmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns hmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from hmenu

val = get(hObject,'Value');
set(hObject,'UserData',val);
% if  val == 1
%     set(handles.hbtxt,'Visible','off');
%     set(handles.hedit_heb,'Visible','off');
% elseif val == 2
%     set(handles.hbtxt,'Visible','on');
%     set(handles.hedit_heb,'Visible','on');
% else
%     msgbox('�����ܳ��ֵİ�');
% end


% --- Executes on button press in chkbx.
function chkbx_Callback(hObject, eventdata, handles)
% hObject    handle to chkbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbx

val = get(hObject,'Value');
set(hObject,'UserData',val);
if val==0
    set(handles.hedit_info,'string',{' ';'�������Ϣ';' '});
else
    set(handles.hedit_info,'string',{' ';'��δ���зָ����';' '});
end


% --------------------------------------------------------------------
function hinfosave_Callback(hObject, eventdata, handles)
% hObject    handle to hinfosave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

infostr = get(handles.hedit_info,'string');
[fname,path] = uiputfile('*.txt','����Ϊ');
if path==0  % ȡ���ļ��������
    return;
end
fid = fopen([path fname],'w');
if fid == -1    % ���ܹ����ļ�
    msgbox({'���ܴ��ļ�' fname},'�ļ��������','error');
    return;
end

% ������Ϣ
fprintf(fid,'\r\n=====%s=====\r\n\r\n',datestr(now,31));
for k = 1:size(infostr,1)
    fprintf(fid,'%s',infostr(k,:));
    fprintf(fid,'\r\n');
end
fclose(fid);



function xieBeni = XieBeniInverted(U, uExpoent, center, data)
%XIEBENIINVERTED Implementation of the measure of cluster validation of the Xie-Beni.

% Calculate Xie-Beni Inverted cluster validation method

dist = sum(sum((U.^uExpoent) .* (distfcm(center, data).^2)));
minimum = intmax;
for i = 1 : size(center, 1)-1
    minimumAux = min(distfcmfp(center(i, :), center(i+1:end, :)).^2);
    if(minimum > minimumAux)
        minimum = minimumAux;
    end
end
xieBeni = dist / (size(data, 1) .* minimum);


function out = distfcm(center, data)
% �������������������ĵľ���
% ���룺
%   center     ---- ��������
%   data       ---- ������
% �����
%   out        ---- ����
out = zeros(size(center, 1), size(data, 1));
for k = 1:size(center, 1), % ��ÿһ����������
    % ÿһ��ѭ��������������㵽һ���������ĵľ���
   % out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2),2)',1);
     out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
end


function out = distfcmfp(center, data)
%DISTFCMFP Distance measure in fuzzy c-mean clustering with focal point.
%	OUT = DISTFCMFP(CENTER, DATA) calculates the Euclidean distance
%	between each row in CENTER and each row in DATA, and returns a
%	distance matrix OUT of size M by N, where M and N are row
%	dimensions of CENTER and DATA, respectively, and OUT(I, J) is
%	the distance between CENTER(I,:) and DATA(J,:).
out = zeros(size(center, 1), size(data, 1));
% fill the output matrix
if size(center, 2) > 1,
    for k = 1:size(center, 1),
        out(k, :) = sqrt(sum(((data-ones(size(data, 1), 1)*center(k, :)).^2)'));
    end
else	% 1-D data
    for k = 1:size(center, 1),
        out(k, :) = abs(center(k)-data)';
    end
end
