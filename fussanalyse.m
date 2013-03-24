clear all;
close all;
clc;

%% Daten Einlesen
[ndata, text, alldata] = xlsread('list.xlsx');

%% Ordner für Bilder erstellen
% mkdir('Bilder2')

%% Programm-/Plotauswahl
% Programm
hauptprogramm       = 1;
test_schwellenwerte = 0;
klassifizierung     = 0;
darstellung         = 0;

%% Auswahl Plot
original                            = 1;
rgb_einzeln                         = 1;
histogramm_einzelfarben             = 1;
pic_sw_und_smoothed                 = 1;
histogramm_fuss                     = 1;
subplot_roh_belastet_klassifiziert  = 1;
fussteilung                         = 1;

start_pic   = 1;
end_pic     = 78;
%% $$$$$$$Hauptprogramm$$$$$$$
if hauptprogramm == 1
HIST_RGB                        = zeros(78,255,3);
BEL_MITTELFUSS                  = zeros(78,1);
BEL_VORFUSS                     = zeros(78,1);
ZUORDUNG_FLAECHENVERHAELTNIS    = zeros(78,4);
RECHTS_LINKS                    = zeros(78,1);


for picture = start_pic:end_pic
    
        pic_org = imread(['FootPower/',text{picture,1}]);

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

%% Histogramm der einzelnen Farben
for col = 1:D_size(2) % Schleife für Pixelspalten
    for row = 1:D_size(1) % Schleife für Pixelzeilen
        for i = 1:254 % Schleife für Helligkeitswerte   
            if pic_org(row,col,1) == i
            HIST_RGB(picture,i+1,1) =   HIST_RGB(picture,i+1,1)+1;
            elseif pic_org(row,col,2) == i
            HIST_RGB(picture,i+1,2) =   HIST_RGB(picture,i+1,2)+1;
            elseif pic_org(row,col,3) == i
            HIST_RGB(picture,i+1,3) =   HIST_RGB(picture,i+1,3)+1;
            end
        end
    end
end

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
HIST = zeros(255,3);
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
            

%% Median
% median = 1;
% summe =0;
% while summe < (sum(HIST(:,2))/2)
%     summe = summe+HIST(median,2);
%     median = median+1;
% end
% 
% % Maximas
% max_first = max( HIST(1:median,2) );
% [val_max_first idx_max_first] = find(HIST(:,2) == max_first);
% 
% max_second = max( HIST(median:end,2) );
% [val_max_second idx_max_second] = find(HIST(:,2) == max_second);
% 
% % Minimum zwischen Maximas
% minmin = min( HIST(val_max_first:val_max_second,2) );
% [val_min_temp idx_min] = find(HIST(val_max_first:val_max_second,2) == minmin);
% 
% val_min = val_min_temp + val_max_first - 1;

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
    RECHTS_LINKS(picture) = 1;
else
    RECHTS_LINKS(picture) = 2;
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

% % Unten rechts
% hits_rechts_unten = zeros( D_size(1), 1);
% hits_rechts_unten( 1:start_unten ) = 999;
% 
% for i = start_unten:D_size(1) 
%     j = D_size(2);
%     while  j > 1 && pic_finish(i,j,3) ~= 255
%                    
%            hits_rechts_unten(i) = hits_rechts_unten(i)+1;
%            j = j-1;
%      end
% end
% [spalte_rechts_unten_temp zeile_rechts_unten] = min(hits_rechts_unten);
% 
% spalte_rechts_unten = D_size(2) - spalte_rechts_unten_temp;

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
BEL_MITTELFUSS(picture,1) = 0;
else
BEL_MITTELFUSS(picture,1) = (max(m2) - min(m2));

x1 = min(m2);
[zeile_m1 none] = find(pic_finish(m_z,x1,3) == 255);
y = m_z + zeile_m1(1);

x2 = max(m2);
end

% Belastung Vorfuß
if isempty(v1)
BEL_VORFUSS(picture,1) = (max(v4) - min(v4)); 
elseif isempty(v3)
BEL_VORFUSS(picture,1) = (max(v2) - min(v2));
elseif v1 == v3
BEL_VORFUSS(picture,1) = (max(v4) - min(v4));    
end

% Belastung Vorfuß
if isempty(h1)
BEL_RUECKFUSS(picture,1) = (max(h4) - min(h4)); 
elseif isempty(h3)
BEL_RUECKFUSS(picture,1) = (max(h2) - min(h2));
elseif h1 == h3
BEL_RUECKFUSS(picture,1) = (max(h4) - min(h4));    
end

ZUORDUNG_FLAECHENVERHAELTNIS(picture,2) = (BEL_MITTELFUSS(picture,1) / BEL_VORFUSS(picture,1))*100;
ZUORDUNG_FLAECHENVERHAELTNIS(picture,3) = (BEL_MITTELFUSS(picture,1) / BEL_RUECKFUSS(picture,1))*100;
ZUORDUNG_FLAECHENVERHAELTNIS(picture,4) = (BEL_RUECKFUSS(picture,1) / BEL_VORFUSS(picture,1))*100;



%% $$$$$$$Plot$$$$$$$ 
if darstellung == 1

if original == 1
figure;
image(pic_org)
axis image
axis off
title('Originalbild','FontSize',16)
end

if rgb_einzeln == 1
figure;
subplot(1,3,1)
image(pic_org)
hold on
axis image
axis off
image(R)
title('Rot','FontSize',16)

subplot(1,3,2)
image(pic_org)
hold on
axis image
axis off
image(G)
title('Grün','FontSize',16)

subplot(1,3,3)
image(pic_org)
hold on
axis image
axis off
image(B)
title('Blau','FontSize',16)
end

if histogramm_einzelfarben == 1
figure
hold on
plot(HIST_RGB(:,1),'r')
plot(HIST_RGB(:,2),'g')
plot(HIST_RGB(:,3),'b')

title('Histogramm Einzelfarben','FontSize',16)
xlabel('Farbwert','FontSize',12)
ylabel('Häufigkeit','FontSize',12)
end

if pic_sw_und_smoothed == 1
figure
subplot(1,2,1)
image(pic_sw) 
axis image
axis off
title('Abgrenzung Roh','FontSize',16)

subplot(1,2,2)
image(pic_sw_smooth) 
axis off
axis image
title('Abgrenzung Verbessert','FontSize',16)

end

if histogramm_fuss == 1
figure
hold on
% plot(HIST(:,1),'r')
plot(HIST(:,2),'g')
% plot(HIST(:,3),'b')

plot(median,HIST(median,2),'or')

plot(val_max_first,HIST(val_max_first,2),'ob')
plot(val_max_second,HIST(val_max_second,2),'ob')
plot(val_min,HIST(val_min,2),'ok')

title('Histogramm Grünwert','FontSize',16)
xlabel('Farbwert','FontSize',12)
ylabel('Häufigkeit','FontSize',12)
legend('Verteilung Farbwert','Median','Erstes Maximum','Zweites Maximum','Minimum')
end

if subplot_roh_belastet_klassifiziert == 1
figure
subplot(1,2,1)
image(pic_org_bel) 
axis image
axis off
title('Unklassifiziert','FontSize',16)

subplot(1,2,2)
image(pic_finish_temp)
axis image
axis off
title('Klassifiziert','FontSize',16)
% Fuss unbelastet
figure(picture+10)
subplot(1,2,1)
image(pic_org) 
axis image
axis off
title('Originalbild','FontSize',16)
% Fuss belastet
figure(picture+10)
subplot(1,2,2)
image(pic_finish)
axis image
hold on
axis off
title('Belastungsflächen','FontSize',16)
arrow([spalte_links_oben zeile_links_oben],[spalte_rechts_oben zeile_rechts_oben],'Length',10,'Ends',[1 2])
% arrow([spalte_links_unten zeile_links_unten],[spalte_rechts_unten zeile_rechts_unten],'Length',10,'Ends',[1 2])
if isempty(m1) && isempty(m2)
else
arrow([x1 y],[x2 y],'Length',10,'Ends',[1 2])
end

end

if fussteilung == 1
% Fuss unbelastet
figure
image(pic_finish) 
axis image
axis off
title('Einteilung Fußareale','FontSize',16)
hold on
% Vorfussachse
plot(spalte_links_oben,zeile_links_oben,'og')
plot(spalte_rechts_oben,zeile_rechts_oben,'og')
plot([spalte_links_oben spalte_rechts_oben],[zeile_links_oben zeile_rechts_oben],'-r')

% Rueckfussachse
plot(spalte_links_unten,zeile_links_unten,'og')
plot(spalte_rechts_unten,zeile_rechts_unten,'og')
plot([spalte_links_unten spalte_rechts_unten],[zeile_links_unten zeile_rechts_unten],'-r')

% Verbindungslinien Vorfuss-Rueckfuss
plot([spalte_links_oben spalte_links_unten],[zeile_links_oben zeile_links_unten],'-r');
plot([spalte_rechts_oben spalte_rechts_unten],[zeile_rechts_oben zeile_rechts_unten],'-r');

% Mittelachse
plot(m_s_l,m_z_l,'og')
plot(m_s_r,m_z_r,'og')
plot([m_s_l m_s_r],[m_z_l m_z_r ],'-r');
end
% saveas(gcf,['Bilder2\',num2str(picture)])

end
end
% save('ZUORDUNG_FLAECHENVERHAELTNIS')
end

%% $$$$$$$Auswertung$$$$$$$
%% Auswertung / Test Schwellwerte$$$$$$$
if test_schwellenwerte == 1
load ZUORDUNG_FLAECHENVERHAELTNIS.mat
% Ergebnismatrizen
Ergebnis_gesamt = zeros(5000,3);
Ergebnis_hohl   = zeros(5000,3);
Ergebnis_senk   = zeros(5000,3);
Ergebnis_andere = zeros(5000,3);

% Zähler
w = 1;
nr = 2;
% Schwellwerte Start-Stop
Schwelle_Hohlfuss_end = 45;
Schwelle_Senkfuss_start = 55;

% Anzahl Fuesse und einzelne Fussarten
n_ges = length(find(isnan(ndata(:,1)) == 0) );
n_hohl = length(find(ndata(:,1) == 1) );
n_senk = length(find(ndata(:,1) == 2) );
n_andere = length(find(ndata(:,1) == 3) );

Schwelle_Hohlfuss = 36;
Schwelle_Senkfuss = 55;
% for Schwelle_Hohlfuss = Schwelle_Hohlfuss_end:-1:5
%     
% for Schwelle_Senkfuss = Schwelle_Senkfuss_start:1:100

%% Zuordnung Pathologie
for picture = start_pic: end_pic

if  (ZUORDUNG_FLAECHENVERHAELTNIS(picture,2) <= Schwelle_Hohlfuss)
ZUORDUNG_FLAECHENVERHAELTNIS(picture,1) = 1;
elseif  ZUORDUNG_FLAECHENVERHAELTNIS(picture,2) >= Schwelle_Senkfuss
ZUORDUNG_FLAECHENVERHAELTNIS(picture,1) = 2;
elseif ZUORDUNG_FLAECHENVERHAELTNIS(picture,2) > Schwelle_Hohlfuss && ZUORDUNG_FLAECHENVERHAELTNIS(picture,2) < Schwelle_Senkfuss
ZUORDUNG_FLAECHENVERHAELTNIS(picture,1) = 3;
end

end


%% Trefferquote
% % Ergebnisse Schwelle Hohl
% right_hohl = zeros(78,1);
% for iii = start_pic: end_pic
%     if ZUORDUNG_FLAECHENVERHAELTNIS(iii,1) == ndata(iii,nr) && ndata(iii,nr) == 1
%         right_hohl(iii,1) = 1;
%     else
%         right_hohl(iii,1) = 0;
%     end
% end
% 
% sum_right_hohl = sum( right_hohl(:,1) );
% Ergebnis_hohl(w,1)  =(sum_right_hohl/n_hohl)*100;
% Ergebnis_hohl(w,2)  = Schwelle_Hohlfuss;
% Ergebnis_hohl(w,3)  = Schwelle_Senkfuss;
% 
% % Ergebnisse Schwelle Senk
% right_senk = zeros(78,1);
% for iii = start_pic: end_pic
%     if ZUORDUNG_FLAECHENVERHAELTNIS(iii,1) == ndata(iii,nr) && ndata(iii,nr) == 2
%         right_senk(iii,1) = 1;
%     else
%         right_senk(iii,2) = 0;
%     end
% end
% 
% sum_right_senk = sum( right_senk(:,1) );
% Ergebnis_senk(w,1)  =(sum_right_senk/n_senk)*100;
% Ergebnis_senk(w,2)  = Schwelle_Hohlfuss;
% Ergebnis_senk(w,3)  = Schwelle_Senkfuss;
% 
% % Ergebnisse Schwelle Andere
% right_andere = zeros(78,1);
% for iii = start_pic: end_pic
%     if ZUORDUNG_FLAECHENVERHAELTNIS(iii,1) == ndata(iii,nr)  && ndata(iii,nr) == 3
%         right_andere(iii,1) = 1;
%     else
%         right_andere(iii,2) = 0;
%     end
% end
% 
% sum_right_andere = sum( right_andere(:,1) );
% Ergebnis_andere(w,1)  =(sum_right_andere/n_andere)*100;
% Ergebnis_andere(w,2)  = Schwelle_Hohlfuss;
% Ergebnis_andere(w,3)  = Schwelle_Senkfuss;

% Ergebnisse Schwelle Gesamt

right_ges = zeros(78,1);
for iii = start_pic: end_pic
    if ZUORDUNG_FLAECHENVERHAELTNIS(iii,1) == ndata(iii,nr)
        right_ges(iii,1) = 1;
    else
        right_ges(iii,1) = 0;
    end
end

sum_right_ges = sum( right_ges(:,1) );
% Ergebnis_gesamt(w,1)  =(sum_right_ges/n_ges)*100
Ergebnis_schwellen =(sum_right_ges/n_ges)*100
Ergebnis_gesamt(w,2)  = Schwelle_Hohlfuss;
Ergebnis_gesamt(w,3)  = Schwelle_Senkfuss;

w = w+1;
% end
% end

[val_ges idx_ges]       = max(Ergebnis_gesamt(:,1));
[val_hohl idx_hohl]     = max(Ergebnis_hohl(:,1));
[val_senk idx_senk]     = max(Ergebnis_senk(:,1));
[val_andere idx_andere] = max(Ergebnis_andere(:,1));
end

%% Classify Pathologie

if klassifizierung == 1
load DATA.mat
nr = 2;
n_ges = length(find(isnan(ndata(:,1)) == 0) );
right_ges_classify = zeros(78,1);

for picture = start_pic:end_pic

ndata_new = ndata;
ZUORDUNG_FLAECHENVERHAELTNIS_cut    = ZUORDUNG_FLAECHENVERHAELTNIS;
% BEL_MITTELFUSS_cut                  = BEL_MITTELFUSS;
% BEL_RUECKFUSS_cut                   = BEL_RUECKFUSS;
% BEL_VORFUSS_cut                     = BEL_VORFUSS;
% RECHTS_LINKS_cut                    = RECHTS_LINKS;

xx(1,1) = ZUORDUNG_FLAECHENVERHAELTNIS(picture,2);
% xx(1,2) = ZUORDUNG_FLAECHENVERHAELTNIS(picture,3);
% xx(1,3) = ZUORDUNG_FLAECHENVERHAELTNIS(picture,4);
% xx(1,4) = BEL_MITTELFUSS(picture);
% xx(1,5) = BEL_RUECKFUSS(picture);
% xx(1,6) = BEL_VORFUSS(picture);
% xx(1,7) = RECHTS_LINKS(picture);

ndata_new(picture,:) = [];

ZUORDUNG_FLAECHENVERHAELTNIS_cut(picture,:) = [];
% BEL_MITTELFUSS_cut(picture,:)               = [];
% BEL_RUECKFUSS_cut(picture,:)                = [];
% BEL_VORFUSS_cut(picture,:)                  = [];  
% RECHTS_LINKS_cut(picture,:)                 = [];

xxgroup         = ndata_new(:,nr);

xxtraining(:,1) = ZUORDUNG_FLAECHENVERHAELTNIS_cut(:,2);
% xxtraining(:,2) = ZUORDUNG_FLAECHENVERHAELTNIS_cut(:,3);
% xxtraining(:,3) = ZUORDUNG_FLAECHENVERHAELTNIS_cut(:,4);
% xxtraining(:,4) = BEL_MITTELFUSS_cut(:);
% xxtraining(:,5) = BEL_RUECKFUSS_cut(:);
% xxtraining(:,6) = BEL_VORFUSS_cut(:);
% xxtraining(:,7) = RECHTS_LINKS_cut(:);

[xxclass(picture),xxerr(picture),xxPosterior(picture,1:3),xxlogp,xxcoeff] = classify(xx,xxtraining,xxgroup);

% for iii = start_pic: end_pic
    if xxclass(picture) == ndata(picture,nr)
        right_ges_classify(picture,1) = 1;
    else
        right_ges_classify(picture,1) = 0;
    end
% end



end
sum_right_ges_classify = sum( right_ges_classify(:,1) );
Ergebnis_gesamt_classify  =(sum_right_ges_classify/n_ges)*100

FINISH(:,1) = ndata(:,2);
FINISH(:,2) = xxclass(:);
FINISH(:,3) = right_ges_classify(:);
FINISH(:,4) = xxerr(:);
FINISH(:,5:7) = xxPosterior(:,1:3);
end

