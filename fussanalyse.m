clear all;
close all;
clc;

[ndata, text, alldata] = xlsread('list.xlsx');

%% Auswahl Plot
original                            = 0;
rgb_einzeln                         = 0;
histogramm_einzelfarben             = 0;
pic_sw_und_smoothed                 = 0;
histogramm_fuss                     = 0;
subplot_roh_belastet_klassifiziert  = 0;

test = zeros(78,1);

for ii = 1:78
    pic_org = imread(['FootPower/',text{ii,1}]);
%% read image
%  pic_org = imread('FootPower/10-259-L.jpg'); 

% Hohlfuss
% pic_org = imread('FootPower/928-324-L.jpg');      %+
% pic_org = imread('FootPower/928-324-R.jpg');      %+
% 
% pic_org = imread('FootPower/928-323-L.jpg');      %+
% pic_org = imread('FootPower/928-323-R.jpg');      %+
% 
% pic_org = imread('FootPower/10-272-L.jpg');       %+
% pic_org = imread('FootPower/10-272-R.jpg');       %+
% 
% pic_org = imread('FootPower/10-270-L.jpg');       %+
% pic_org = imread('FootPower/10-270-R.jpg');       %+
% 
% pic_org = imread('FootPower/10-266-L.jpg');       %+ platt/senk statt hohl      
% pic_org = imread('FootPower/10-266-R.jpg');       %+ platt/senk statt hohl 
% 
% pic_org = imread('FootPower/10-265-L.jpg');       %+ 
% pic_org = imread('FootPower/10-265-R.jpg');       %+        
% 
% pic_org = imread('FootPower/10-263-L.jpg');       %-
% pic_org = imread('FootPower/10-263-R.jpg');       %+-
% 
% pic_org = imread('FootPower/10-262-L.jpg');       %+
% pic_org = imread('FootPower/10-262-R.jpg');       %+
% 
% pic_org = imread('FootPower/10-261-L.jpg');       %+
% pic_org = imread('FootPower/10-261-R.jpg');       %+
% 
% pic_org = imread('FootPower/10-260-L.jpg');       %+
% pic_org = imread('FootPower/10-260-R.jpg');       %+-
% 
% pic_org = imread('FootPower/928-321-L.jpg');      %-
% pic_org = imread('FootPower/928-321-R.jpg');      %+-
% 
% pic_org = imread('FootPower/10-258-L.jpg');       %-
% pic_org = imread('FootPower/10-258-R.jpg');       %+
% 
% pic_org = imread('FootPower/10-257-L.jpg');       %+
% pic_org = imread('FootPower/10-257-R.jpg');
% 
% pic_org = imread('FootPower/10-254-L.jpg');       %+
% pic_org = imread('FootPower/10-254-R.jpg');       %+ platt/senk statt hohl    
% 
% pic_org = imread('FootPower/10-251-L.jpg');       %+-
% pic_org = imread('FootPower/10-251-R.jpg');       %+-
% 
% pic_org = imread('FootPower/10-250-L.jpg');       %-
% pic_org = imread('FootPower/10-250-R.jpg');       %-
% 
% pic_org = imread('FootPower/10-249-L.jpg');       %+
% pic_org = imread('FootPower/10-249-R.jpg');       %+-
% 
% pic_org = imread('FootPower/10-248-L.jpg');       %+
% pic_org = imread('FootPower/10-248-R.jpg');       %+
% 
% pic_org = imread('FootPower/10-246-L.jpg');       %+
% pic_org = imread('FootPower/10-246-R.jpg');       %+
% 
% pic_org = imread('FootPower/10-244-L.jpg');       %+
% pic_org = imread('FootPower/10-244-R.jpg');       %+
% 
% pic_org = imread('FootPower/928-320-L.jpg');      %+
% pic_org = imread('FootPower/928-320-R.jpg');      %+
% 
% pic_org = imread('FootPower/10-240-L.jpg');       %+
% pic_org = imread('FootPower/10-240-R.jpg');       %+

% Senkfuss
% pic_org = imread('FootPower/10-268-L.jpg');       %+
% pic_org = imread('FootPower/10-268-R.jpg');       %+-
% 
% pic_org = imread('FootPower/928-322-L.jpg');      %+
% pic_org = imread('FootPower/928-322-R.jpg');      %+
% 
% pic_org = imread('FootPower/10-255-L.jpg');       %+
% pic_org = imread('FootPower/10-255-R.jpg');       %+
% 
% pic_org = imread('FootPower/10-253-L.jpg');       %+
% pic_org = imread('FootPower/10-253-R.jpg');       %+
% 
% pic_org = imread('FootPower/10-245-L.jpg');       %+
% pic_org = imread('FootPower/10-245-R.jpg');       %+
% 
% pic_org = imread('FootPower/10-242-L.jpg');       %-
% pic_org = imread('FootPower/10-242-R.jpg');       %+-
% 
% pic_org = imread('FootPower/10-241-L.jpg');       %+
% pic_org = imread('FootPower/10-241-R.jpg');       %+

% Plattfuß
% pic_org = imread('FootPower/10-267-L.jpg');       %+
% pic_org = imread('FootPower/10-267-R.jpg');       %+
% 
% pic_org = imread('FootPower/10-256-L.jpg');       %+
% pic_org = imread('FootPower/10-256-R.jpg');       %+

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

%% Histogramm der Einzelnen Farben
% HIST_RGB = zeros(255,3);
% 
% for col = 1:D_size(2)
%     for row = 1:D_size(1)
%         for i = 1:254
%             if pic_org(row,col,1) == i
%             HIST_RGB(i+1,1) =   HIST_RGB(i+1,1)+1;
%             elseif pic_org(row,col,2) == i
%             HIST_RGB(i+1,2) =   HIST_RGB(i+1,2)+1;
%             elseif pic_org(row,col,3) == i
%             HIST_RGB(i+1,3) =   HIST_RGB(i+1,3)+1;
%             end
%         end
%     end
% end

%% Fuß in weiß darstellen

schwelle = 15;
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
            
%% Schwelle für belastete Flächen

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

%% Test

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


%% Classify

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

group = zeros(xtraining_size,1);

n = 1;
o = 1;

for m = 1:( D_size(1)*D_size(2) )
        if pic_smooth_vec(m,1) == 255
            
            if  pic_org_bel_vec(m,3) == 255
                
            xtraining(n,:) = pic_org_vec(m,:);
            group(n) = 1;
            
            n = n+1;
            
            elseif  pic_org_bel_vec(m) == pic_org_vec(m)
                
            xtraining(n,:) = pic_org_vec(m,:);
            group(n) = 2;
            
            n = n+1;
            
            elseif pic_org_bel_vec(m,1) == 255
            
            x(o,:) = pic_org_vec(m,:);
            o = o+1;    
            end
        end
end


[class,err,Posterior,logp,coeff] = classify(x,xtraining,group);

pic_finish_vec = zeros(D_size(1)*D_size(2),1,'uint8');

n = 1;
for m = 1:( D_size(1)*D_size(2) )
        if pic_smooth_vec(m,1) == 255
            
            if  pic_org_bel_vec(m,1) == 255
                
                if class(n) == 1
                    
                    pic_finish_vec(m,1) = 0;
                    pic_finish_vec(m,2) = 0;
                    pic_finish_vec(m,3) = 255;
                    
                    n = n+1;
                    
                elseif class(n) == 2
                
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

%% Pathologie

%% Achse Vorfuß
% Oben links
start_oben    = round(D_size(1)*0.25);
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

%% Achse Rückfuß
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

%% Mittelachse
m_z_l = round(mean((zeile_links_oben+zeile_links_unten)/2));
m_s_l = round(mean((spalte_links_oben+spalte_links_unten)/2));

m_z_r = round(mean((zeile_rechts_oben+zeile_rechts_unten)/2));
m_s_r = round(mean((spalte_rechts_oben+spalte_rechts_unten)/2));

BEL_VORFUSS = spalte_rechts_oben - spalte_links_oben;


[x1 x2] = find(pic_finish(m_z_l:m_z_r,1:end,3) == 255);       
[x3 x4] = find(pic_finish(m_z_r:m_z_l,1:end,3) == 255);

if isempty(x1) && isempty(x2)&& isempty(x3)&& isempty(x4)
% disp('Hohlfuss')
test(ii,1) = 1;
elseif isempty(x1)
BEL_MITTELFUSS = (max(x4) - min(x4));
elseif isempty(x3)
BEL_MITTELFUSS = (max(x2) - min(x2));
end

test(ii,2) = (BEL_MITTELFUSS / BEL_VORFUSS)*100;

if  test(ii,2) <= 40
%     disp('Hohlfuss')
test(ii,1) = 1;
elseif  test(ii,2) > 40
%     dips('Platt-/Senkfuss')
test(ii,1) = 2;
end

end

%% Plot
if original == 1
figure;
image(pic_org)
axis equal
end

if rgb_einzeln == 1
figure;
subplot(1,3,1)
image(pic_org)
hold on
axis equal
image(R)

subplot(1,3,2)
image(pic_org)
hold on
axis equal
image(G)

subplot(1,3,3)
image(pic_org)
hold on
axis equal
image(B)
end

if histogramm_einzelfarben == 1
figure
hold on
plot(HIST_RGB(:,1),'r')
plot(HIST_RGB(:,2),'g')
plot(HIST_RGB(:,3),'b')
end

if pic_sw_und_smoothed == 1
figure(1)
subplot(1,2,1)
image(pic_sw) 
axis equal

figure(1)
subplot(1,2,2)
image(pic_sw_smooth) 
axis equal
end

if histogramm_fuss == 1
figure
hold on
plot(HIST(:,1),'r')
plot(HIST(:,2),'g')
plot(HIST(:,3),'b')

plot(median,HIST(median,2),'or')

plot(val_max_first,HIST(val_max_first,2),'ob')
plot(val_max_second,HIST(val_max_second,2),'ob')
plot(val_min,HIST(val_min,2),'ok')
end

if subplot_roh_belastet_klassifiziert == 1
close all
figure(2)
subplot(1,4,1)
image(pic_org) 
axis equal

figure(2)
subplot(1,4,2)
image(pic_org_bel) 
axis equal

figure(2)
subplot(1,4,3)
image(pic_finish_temp)
axis equal

figure(2)
subplot(1,4,4)
image(pic_finish)
axis equal

hold on
plot(spalte_links_oben,zeile_links_oben,'og')
plot(spalte_rechts_oben,zeile_rechts_oben,'og')
plot([spalte_links_oben spalte_rechts_oben],[zeile_links_oben zeile_rechts_oben],'-r')

plot(spalte_links_unten,zeile_links_unten,'og')
plot(spalte_rechts_unten,zeile_rechts_unten,'og')
plot([spalte_links_unten spalte_rechts_unten],[zeile_links_unten zeile_rechts_unten],'-r')

plot([spalte_links_oben spalte_links_unten],[zeile_links_oben zeile_links_unten],'-r');
plot([spalte_rechts_oben spalte_rechts_unten],[zeile_rechts_oben zeile_rechts_unten],'-r');

plot([m_s_l m_s_r],[m_z_l m_z_r ],'-r');
end