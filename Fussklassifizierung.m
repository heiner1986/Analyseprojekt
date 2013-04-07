function varargout = Fussklassifizierung(varargin)
% FUSSKLASSIFIZIERUNG MATLAB code for Fussklassifizierung.fig
%      FUSSKLASSIFIZIERUNG, by itself, creates a new FUSSKLASSIFIZIERUNG or raises the existing
%      singleton*.
%
%      H = FUSSKLASSIFIZIERUNG returns the handle to a new FUSSKLASSIFIZIERUNG or the handle to
%      the existing singleton*.
%
%      FUSSKLASSIFIZIERUNG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FUSSKLASSIFIZIERUNG.M with the given input arguments.
%
%      FUSSKLASSIFIZIERUNG('Property','Value',...) creates a new FUSSKLASSIFIZIERUNG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Fussklassifizierung_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Fussklassifizierung_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Fussklassifizierung

% Last Modified by GUIDE v2.5 20-Mar-2013 18:29:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Fussklassifizierung_OpeningFcn, ...
                   'gui_OutputFcn',  @Fussklassifizierung_OutputFcn, ...
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


% --- Executes just before Fussklassifizierung is made visible.
function Fussklassifizierung_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Fussklassifizierung (see VARARGIN)

% Choose default command line output for Fussklassifizierung
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Fussklassifizierung wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Fussklassifizierung_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Dateiauswahl.
function Dateiauswahl_Callback(hObject, eventdata, handles)
% hObject    handle to Dateiauswahl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global filename;
global pathname;

[filename,pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
          '*.*','All Files' });
set(handles.Dateiname,'string',filename);



function Dateiname_Callback(hObject, eventdata, handles)
% hObject    handle to Dateiname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dateiname as text
%        str2double(get(hObject,'String')) returns contents of Dateiname as a double


% --- Executes during object creation, after setting all properties.
function Dateiname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dateiname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global filename;
global pathname;

format short

%% Daten laden
load DATA.mat

%% Bild einlesen
pic_org = imread([pathname filename]);

%% Parameter
D_size = size(pic_org);

%% Fuß in RGB Einzelfarben
R = pic_org;
G = pic_org;
B = pic_org;

R(:,:,2:3) = zeros;
G(:,:,1) = zeros;
G(:,:,3) = zeros;
B(:,:,1:2) = zeros;

%% Fuß in weiß darstellen

schwelle = 20;
pic_sw = pic_org;

for i = 1:D_size(1)
    for j = 1:D_size(2)
     
        if pic_sw(i,j,3) > schwelle;

            pic_sw(i,j,1) = 255;
            pic_sw(i,j,2) = 255;
            pic_sw(i,j,3) = 255;  

        else
            pic_sw(i,j,1) = 0;
            pic_sw(i,j,2) = 0;
            pic_sw(i,j,3) = 0;
        end
    end
end

%% Ränder glätten
pic_sw_smooth = pic_sw;

for col = 2:( D_size(2)-1 )
    for row = 2: ( D_size(1)-1 )
        summe = 0;
        idx       = zeros(1,8);
        idx(1)    = pic_sw_smooth( ( row-1),(col-1),1 );
        idx(2)    = pic_sw_smooth( ( row  ),(col-1),1 );
        idx(3)    = pic_sw_smooth( ( row+1),(col-1),1 );
        idx(4)    = pic_sw_smooth( ( row-1),(col  ),1 );
        idx(5)    = pic_sw_smooth( ( row+1),(col  ),1 );
        idx(6)    = pic_sw_smooth( ( row-1),(col+1),1 );
        idx(7)    = pic_sw_smooth( ( row  ),(col+1),1 );
        idx(8)    = pic_sw_smooth( ( row+1),(col+1),1 );
        
        summe = sum(idx);
        
        if summe < 1020
            pic_sw_smooth(row,col,1) = 0;
            pic_sw_smooth(row,col,2) = 0;
            pic_sw_smooth(row,col,3) = 0;
        end
    end
end

%% Histogramm der Einzelnen Farben nur auf Fuß
HIST= zeros(255,3);

for col = 1:D_size(2)
    for row = 1:D_size(1)  
        for i = 1:254
            if pic_sw_smooth(row,col,3) == 255 && pic_sw_smooth(row,col,2) == 255 && pic_sw_smooth(row,col,3) == 255 
                if pic_org(row,col,1) == i
                HIST(i+1,1) =   HIST(i+1,1)+1;
                elseif pic_org(row,col,2) == i
                HIST(i+1,2) =   HIST(i+1,2)+1;
                elseif pic_org(row,col,3) == i
                HIST(i+1,3) =   HIST(i+1,3)+1;
                end
            end
        end  
    end
end

%% Median und Schwellen für Belastete Flaeche

HIST(:,2) = smooth(HIST(:,2),5);

median = 1;
summe =0;
while summe < (sum(HIST(:,2))/2)
    summe = summe+HIST(median,2);
    median = median+1;
end

% Maximas Gruen
max_first = max( HIST(1:median,2) );
[val_max_first idx_max_first] = find(HIST(:,2) == max_first);

max_second = findpeaks( HIST(median:end,2) );
[val_max_second_temp idx_max_second] = find(HIST(:,2) == max_second(1));

val_max_second = max(val_max_second_temp);

% Bereich um Masimum Blau
max_blue = max( HIST(:,3) );
[val_max_blue_temp idx_max_blue] = find(HIST(:,3) == max_blue);

val_max_blue=max(val_max_blue_temp);

% Minimum zwischen Maximas Gruen
minmin = min( HIST(val_max_first:val_max_second,2) );
[val_min_temp idx_min] = find(HIST(val_max_first:val_max_second,2) == minmin);

val_min = max(val_min_temp) + val_max_first - 1;


%% Fuß Original und Belastungsflächen Blau darstellen
schwelle_bel_vorne = val_max_first ;
schwelle_bel_hinten = val_max_second;

pic_org_bel = pic_sw_smooth;

fenster_um_min = 15;
fenster_um_max = 5;

for i = 1:D_size(1)
    for j = 1:D_size(2)
        if pic_sw_smooth(i,j,1) == 255
            
            if ( pic_org(i,j,2) > val_min + fenster_um_min ) 
                
                pic_org_bel(i,j,1) = 0;
                pic_org_bel(i,j,2) = 0;
                pic_org_bel(i,j,3) = 255 ;
                
            elseif (pic_org(i,j,2) >= val_min - fenster_um_min && pic_org(i,j,2) <= val_min + fenster_um_min) || ( pic_org(i,j,3) > val_max_blue )
                
                pic_org_bel(i,j,1) = 255;
                pic_org_bel(i,j,2) = 0;
                pic_org_bel(i,j,3) =0 ;
            
            else
                pic_org_bel(i,j,1) = pic_org(i,j,1);
                pic_org_bel(i,j,2) = pic_org(i,j,2);
                pic_org_bel(i,j,3) = pic_org(i,j,3);     

            end
        end
    end
end


%% Classify Belastungsflächen

blue   = find(pic_org_bel(:,:,3) == 255);
origin = find(pic_org_bel == pic_org );
red    = find(pic_org_bel(:,:,1) == 255);

l_blue   = length(blue);
l_origin = length(origin);
l_red    = length(red);

xtraining_size  = l_blue + l_origin;
x_size          = l_red;   

pic_org_bel_vec     = reshape(pic_org_bel, D_size(1)*D_size(2), 3);
pic_org_vec         = reshape(pic_org,D_size(1)*D_size(2),3);
pic_smooth_vec      = reshape(pic_sw_smooth, D_size(1)*D_size(2), 3);

xtraining   = zeros(xtraining_size,3);
x           = zeros(x_size,3);

xgroup = zeros(xtraining_size,1);

n = 1;
o = 1;

for m = 1:( D_size(1)*D_size(2) )
        if pic_smooth_vec(m,1) == 255
            
            if  pic_org_bel_vec(m,3) == 255
                
            xtraining(n,:) = pic_org_vec(m,:);
            xgroup(n) = 1;
            
            n = n+1;
            
            elseif  pic_org_bel_vec(m) == pic_org_vec(m)
                
            xtraining(n,:) = pic_org_vec(m,:);
            xgroup(n) = 2;
            
            n = n+1;
            
            elseif pic_org_bel_vec(m,1) == 255
            
            x(o,:) = pic_org_vec(m,:);
            o = o+1;    
            end
        end
end


[xclass,xerr,xPosterior,xlogp,xcoeff] = classify(x,xtraining,xgroup);

pic_finish_vec = zeros(D_size(1)*D_size(2),1,'uint8');

n = 1;
for m = 1:( D_size(1)*D_size(2) )
        if pic_smooth_vec(m,1) == 255
            
            if  pic_org_bel_vec(m,1) == 255
                
                if xclass(n) == 1
                    
                    pic_finish_vec(m,1) = 0;
                    pic_finish_vec(m,2) = 0;
                    pic_finish_vec(m,3) = 255;
                    
                    n = n+1;
                    
                elseif xclass(n) == 2
                
                    pic_finish_vec(m,1) = pic_org_vec(m,1);
                    pic_finish_vec(m,2) = pic_org_vec(m,2);
                    pic_finish_vec(m,3) = pic_org_vec(m,3);
                    
                    n = n+1;
                end

            else
                    
                   pic_finish_vec(m,1) = pic_org_bel_vec(m,1);
                   pic_finish_vec(m,2) = pic_org_bel_vec(m,2);
                   pic_finish_vec(m,3) = pic_org_bel_vec(m,3);
                   
            end
        end
end

pic_finish_temp = reshape(pic_finish_vec,D_size(1),D_size(2),3);

%% Ränder glätten
pic_finish = zeros(D_size(1), D_size(2), 3, 'uint8');

for col = 2:( D_size(2)-1 )
    for row = 2: ( D_size(1)-1 )
          if pic_finish_temp(row,col,1) > 0 && pic_finish_temp(row,col,2) > 0 && pic_finish_temp(row,col,3) > 0

            idx       = zeros(1,8);
            idx(1)    = pic_finish_temp( ( row-1),(col-1),3 );
            idx(2)    = pic_finish_temp( ( row  ),(col-1),3 );
            idx(3)    = pic_finish_temp( ( row+1),(col-1),3 );
            idx(4)    = pic_finish_temp( ( row-1),(col  ),3 );
            idx(5)    = pic_finish_temp( ( row+1),(col  ),3 );
            idx(6)    = pic_finish_temp( ( row-1),(col+1),3 );
            idx(7)    = pic_finish_temp( ( row  ),(col+1),3 );
            idx(8)    = pic_finish_temp( ( row+1),(col+1),3 );

            summe = sum(idx);

            if pic_finish_temp(row,col,3) == 255 && summe > 1530
                pic_finish(row,col,1) = 0;
                pic_finish(row,col,2) = 0;
                pic_finish(row,col,3) = 255;
                
            else
                pic_finish(row,col,1) = pic_finish_temp(row,col,1);
                pic_finish(row,col,2) = pic_finish_temp(row,col,2);
                pic_finish(row,col,3) = pic_finish_temp(row,col,3);
            end
            if pic_finish_temp(row,col,3) < 255  && summe > 1275

                pic_finish(row,col,1) = 0;
                pic_finish(row,col,2) = 0;
                pic_finish(row,col,3) = 255;

            else    
                pic_finish(row,col,1) = pic_finish_temp(row,col,1);
                pic_finish(row,col,2) = pic_finish_temp(row,col,2);
                pic_finish(row,col,3) = pic_finish_temp(row,col,3);
            end
          else
                pic_finish(row,col,1) = pic_finish_temp(row,col,1);
                pic_finish(row,col,2) = pic_finish_temp(row,col,2);
                pic_finish(row,col,3) = pic_finish_temp(row,col,3);
          end
    end
end

%% Links-Rechts

[spalte_toe_temp zeile_toe] = find(pic_finish(:,:,1) ~= 0);

[zeile_toe_temp spalte_toe ] = find( pic_finish( min(zeile_toe(:)) ,:,1)~= 0,1,'first');

if spalte_toe >= (D_size(2)/2)
    RECHTS_LINKS = 1;
else
    RECHTS_LINKS = 2;
end


%% Pathologie

%% Achse Vorfuß
% Oben links
start_oben    = round(D_size(1)*0.3);
end_oben     = round(D_size(1)*0.5);

hits_links_oben = zeros( D_size(1), 1);
hits_links_oben( 1:start_oben ) = 999;
hits_links_oben( end_oben:end ) = 999;

for i = start_oben : end_oben
    j = 1;
    while j < D_size(2) && pic_finish(i,j,3) ~= 255
                   
           hits_links_oben(i) = hits_links_oben(i)+1;
           j = j+1;
     end
end
[spalte_links_oben zeile_links_oben] = min(hits_links_oben);

% Oben rechts
hits_rechts_oben = zeros( D_size(1), 1);
hits_rechts_oben( 1:start_oben ) = 999;
hits_rechts_oben( end_oben:end ) = 999;

for i = start_oben : end_oben
    j = D_size(2);
    while j > 1 && pic_finish(i,j,3) ~= 255
                   
           hits_rechts_oben(i) = hits_rechts_oben(i)+1;
           j = j-1;
     end
end
[spalte_rechts_oben_temp zeile_rechts_oben] = min(hits_rechts_oben);

spalte_rechts_oben = D_size(2) - spalte_rechts_oben_temp;

differenz_zeilen = abs(zeile_rechts_oben-zeile_links_oben);



%% Achse Rückfuß

if zeile_rechts_oben <= zeile_links_oben
% Unten links
start_unten    = round(D_size(1)*0.9);

hits_links_unten = zeros( D_size(1), 1);
hits_links_unten( 1:start_unten ) = 999;

for i = start_unten:D_size(1) 
    j = 1;
    while j < D_size(2) && pic_finish(i,j,3) ~= 255
                   
           hits_links_unten(i) = hits_links_unten(i)+1;
           j = j+1;
     end
end
[spalte_links_unten zeile_links_unten] = min(hits_links_unten);

% Unten rechts
hits_rechts_unten = zeros( D_size(1), 1);
hits_rechts_unten( 1:end ) = 999;

for i = (zeile_links_unten-differenz_zeilen ):zeile_links_unten
    j = D_size(2);           
    hits_rechts_unten(i) = 0;
    while  j > 1 && pic_finish(i,j,3) ~= 255
           hits_rechts_unten(i) = hits_rechts_unten(i)+1;
           j = j-1;
     end
end
[spalte_rechts_unten_temp zeile_rechts_unten] = min(hits_rechts_unten);
spalte_rechts_unten = D_size(2) - spalte_rechts_unten_temp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if zeile_rechts_oben > zeile_links_oben
% Unten rechts

start_unten    = round(D_size(1)*0.9);
hits_rechts_unten = zeros( D_size(1), 1);
hits_rechts_unten( 1:start_unten ) = 999;

for i = start_unten:D_size(1) 
    j = D_size(2);
    while  j > 1 && pic_finish(i,j,3) ~= 255
                   
           hits_rechts_unten(i) = hits_rechts_unten(i)+1;
           j = j-1;
     end
end
[spalte_rechts_unten_temp zeile_rechts_unten] = min(hits_rechts_unten);
spalte_rechts_unten = D_size(2) - spalte_rechts_unten_temp;

% Unten links


hits_links_unten = zeros( D_size(1), 1);
hits_links_unten( 1:end ) = 999;

for i = (zeile_rechts_unten-differenz_zeilen ): zeile_rechts_unten
    j = 1;           
    hits_links_unten(i) = 0;     
    while j < D_size(2) && pic_finish(i,j,3) ~= 255
  
           hits_links_unten(i) = hits_links_unten(i)+1;
           j = j+1;
     end
end
[spalte_links_unten zeile_links_unten] = min(hits_links_unten);

end

%% Mittelachse
m_z_l = round(mean((zeile_links_oben+zeile_links_unten)/2));
m_s_l = round(mean((spalte_links_oben+spalte_links_unten)/2));

m_z_r = round(mean((zeile_rechts_oben+zeile_rechts_unten)/2));
m_s_r = round(mean((spalte_rechts_oben+spalte_rechts_unten)/2));

m_z = round( (m_z_l+m_z_r)/2);

%% Zurodnung Pathologie

[v1 v2] = find(pic_finish(zeile_links_oben:zeile_rechts_oben,:,3) == 255);       
[v3 v4] = find(pic_finish(zeile_rechts_oben:zeile_links_oben,:,3) == 255);

[m1 m2] = find(pic_finish(m_z,:,3) == 255);

[h1 h2] = find(pic_finish(zeile_links_unten:zeile_rechts_unten,:,3) == 255);       
[h3 h4] = find(pic_finish(zeile_rechts_unten:zeile_links_unten,:,3) == 255);
% [m3 m4] = find(pic_finish(m_z_r:m_z_l,:,3) == 255);

% Belastung Mittelfuß
if isempty(m1) && isempty(m2)
BEL_MITTELFUSS = 0;
else
BEL_MITTELFUSS = (max(m2) - min(m2));

x1 = min(m2);
[zeile_m1 none] = find(pic_finish(m_z,x1,3) == 255);
y = m_z + zeile_m1(1);

x2 = max(m2);
end

% Belastung Vorfuß
if isempty(v1)
BEL_VORFUSS = (max(v4) - min(v4)); 
elseif isempty(v3)
BEL_VORFUSS = (max(v2) - min(v2));
elseif v1 == v3
BEL_VORFUSS = (max(v4) - min(v4));    
end


ZUORDUNG_FLAECHENVERHAELTNIS(1,2) = (BEL_MITTELFUSS / BEL_VORFUSS)*100;
%% plot
% Fuss unbelastet

if RECHTS_LINKS == 2
    side = 'Rechts';
else
    side = 'Links';
end

h = subplot(1,2,1);
cla
image(pic_org) 
axis image
axis off
title('Original','FontSize',16)

% Fuss belastet
h = subplot(1,2,2);
cla
image(pic_finish)
axis image
hold on
axis off
title('Belastet','FontSize',16)

arrow([spalte_links_oben zeile_links_oben],[spalte_rechts_oben zeile_rechts_oben],'Length',10,'Ends',[1 2])
if isempty(m1) && isempty(m2)
else
arrow([x1 y],[x2 y],'Length',10,'Ends',[1 2])
end

%% classify

ndata_new = ndata;
ZUORDUNG_FLAECHENVERHAELTNIS_cut    = ZUORDUNG_FLAECHENVERHAELTNIS;

xx(1,1) = ZUORDUNG_FLAECHENVERHAELTNIS(1,2);

idx = find(ismember(text(:,1), filename)==1);
ndata_new(idx,:) = [];

ZUORDUNG_FLAECHENVERHAELTNIS_cut(1,:) = [];

xxgroup         = ndata_new(:,2);

xxtraining(:,1) = ZUORDUNG_FLAECHENVERHAELTNIS_cut(:,2);

[xxclass,xxerr,xxPosterior(1,1:3),xxlogp,xxcoeff] = classify(xx,xxtraining,xxgroup);

if xxclass == 1
   vorschlag = 'Hohlfusstyp';
elseif xxclass == 2
   vorschlag = 'Platt-/Senkfuss';
elseif xxclass == 3
   vorschlag = 'Normalfuss/Anderer Typ';
end

set(handles.seite,'string', side)

set(handles.verhaeltnis,'string',num2str(ZUORDUNG_FLAECHENVERHAELTNIS(1,2),'%.2f'))

set(handles.phohl,  'string', num2str(xxPosterior(1,1)*100,'%.2f') )
set(handles.pplatt, 'string', num2str(xxPosterior(1,2)*100,'%.2f') )
set(handles.pnormal,'string', num2str(xxPosterior(1,3)*100,'%.2f') )

set(handles.fusstyp,'string', vorschlag )

function verhaeltnis_Callback(hObject, eventdata, handles)
% hObject    handle to verhaeltnis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of verhaeltnis as text
%        str2double(get(hObject,'String')) returns contents of verhaeltnis as a double


% --- Executes during object creation, after setting all properties.
function verhaeltnis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to verhaeltnis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function phohl_Callback(hObject, eventdata, handles)
% hObject    handle to phohl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phohl as text
%        str2double(get(hObject,'String')) returns contents of phohl as a double


% --- Executes during object creation, after setting all properties.
function phohl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phohl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pplatt_Callback(hObject, eventdata, handles)
% hObject    handle to pplatt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pplatt as text
%        str2double(get(hObject,'String')) returns contents of pplatt as a double


% --- Executes during object creation, after setting all properties.
function pplatt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pplatt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pnormal_Callback(hObject, eventdata, handles)
% hObject    handle to pnormal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pnormal as text
%        str2double(get(hObject,'String')) returns contents of pnormal as a double


% --- Executes during object creation, after setting all properties.
function pnormal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pnormal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fusstyp_Callback(hObject, eventdata, handles)
% hObject    handle to fusstyp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fusstyp as text
%        str2double(get(hObject,'String')) returns contents of fusstyp as a double


% --- Executes during object creation, after setting all properties.
function fusstyp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fusstyp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seite_Callback(hObject, eventdata, handles)
% hObject    handle to seite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seite as text
%        str2double(get(hObject,'String')) returns contents of seite as a double


% --- Executes during object creation, after setting all properties.
function seite_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
