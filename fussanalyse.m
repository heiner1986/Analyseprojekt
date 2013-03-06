clear all;
close all;
clc;

%% read image
% pic_org = imread('Knick-senkfuss_re.jpg');
% pic_org = imread('FootPower/10-270-L.jpg');
pic_org = imread('FootPower/10-270-R.jpg');

%% Parameter
D_size = size(pic_org);

%% Plot Original
figure;
image(pic_org)
axis equal

%% Fuﬂ in RGB Einzelfarben
R = pic_org;
G = pic_org;
B = pic_org;

R(:,:,2:3) = zeros;
G(:,:,1) = zeros;
G(:,:,3) = zeros;
B(:,:,1:2) = zeros;

% Plot
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

%% Histogramm der Einzelnen Farben
HIST_RGB = zeros(255,3);

for col = 1:D_size(2)
    for row = 1:D_size(1)
        for i = 1:255
            if pic_org(row,col,1) == i
            HIST_RGB(i+1,1) =   HIST_RGB(i+1,1)+1;
            elseif pic_org(row,col,2) == i
            HIST_RGB(i+1,2) =   HIST_RGB(i+1,2)+1;
            elseif pic_org(row,col,3) == i
            HIST_RGB(i+1,3) =   HIST_RGB(i+1,3)+1;
            end
        end
    end
end

figure
hold on
plot(HIST_RGB(:,1),'r')
plot(HIST_RGB(:,2),'g')
plot(HIST_RGB(:,3),'b')

%% Fuﬂ in schwarz darstellen

schwelle = find(HIST_RGB(:,3) < 50,1,'last')

schwelle = 240;
pic_sw = pic_org;

for i = 1:D_size(1)
    for j = 1:D_size(2)
     
        if pic_sw(i,j,3) > 15;

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

figure
subplot(1,2,1)
image(pic_sw)
axis equal

% R‰nder gl‰tten
pic_smooth = pic_sw;

for col = 2:( D_size(2)-1 )
    for row = 2: ( D_size(1)-1 )
        summe = 0;
        idx       = zeros(1,8);
        idx(1)    = pic_smooth( ( row-1),(col-1),1 );
        idx(2)    = pic_smooth( ( row  ),(col-1),1 );
        idx(3)    = pic_smooth( ( row+1),(col-1),1 );
        idx(4)    = pic_smooth( ( row-1),(col  ),1 );
        idx(5)    = pic_smooth( ( row+1),(col  ),1 );
        idx(6)    = pic_smooth( ( row-1),(col+1),1 );
        idx(7)    = pic_smooth( ( row  ),(col+1),1 );
        idx(8)    = pic_smooth( ( row+1),(col+1),1 );
        
        summe = sum(idx);
        
        if summe < 1020
            pic_smooth(row,col,1) = 0;
            pic_smooth(row,col,2) = 0;
            pic_smooth(row,col,3) = 0;
        end
    end
end

subplot(1,2,2)
image(pic_smooth)
axis equal

%% Histogramm der Einzelnen Farben nur auf Fuﬂ
HIST= zeros(255,3);

for col = 1:D_size(2)
    for row = 1:D_size(1)  
        for i = 1:255
            if pic_smooth(row,col,3) == 255 && pic_smooth(row,col,2) == 255 && pic_smooth(row,col,3) == 255 
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
            
figure
hold on
plot(HIST(:,1),'r')
plot(HIST(:,2),'g')
plot(HIST(:,3),'b')

%% Schwelle f¸r belastete Fl‰chen

% Median

median = 1;
summe =0;
while summe < (sum(HIST(:,2))/2)
    summe = summe+HIST(median,2);
    median = median+1;
end
plot(median,HIST(median,2),'or')

% Maximas
max_first = max( HIST(1:median,2) );
[val_max_first idx_max_second] = find(HIST == max_first);

max_second = max( HIST(median:end,2) );
[val_max_second idx_max_second] = find(HIST == max_second);

% Minimum zwischen Maximas
minmin = min( HIST(val_max_first:val_max_second,2) );
[val_min idx_min] = find(HIST == minmin);

% plot
plot(val_max_first,HIST(val_max_first,2),'ob')
plot(val_max_second,HIST(val_max_second,2),'ob')
plot(val_min,HIST(val_min,2),'ok')

%% Fuﬂ schwarz und Belastungsfl‰chen weiﬂ darstellen
schwelle_bel_vorn = val_min ;
schwelle_bel_hinten = val_min ;
pic_org_bel = pic_smooth;
for i = 1:D_size(1)
    for j = 1:D_size(2)
        if pic_smooth(i,j,1) == 255
            if pic_org(i,j,2) > val_min
                   
                pic_org_bel(i,j,1) = 0;
                pic_org_bel(i,j,2) = 255;
                pic_org_bel(i,j,3) =0 ;
                
            end
        end
    end
end

    
figure
image(pic_org_bel)
axis equal
hold on

%% Histogramm Farbst‰rke
% 
% gruen = G(:,:,2);
% size_gruen = size(gruen);
% gruen_vec = reshape(gruen,1,size_gruen(1)*size_gruen(2));
% gruen_vec_double = cast(gruen_vec,'double')
% 
% figure
% hist(gruen_vec_double,255)

%% Fuﬂ belastet verfeinert darstellen
% pic_org_sw_clean = pic_org_bel;
% 
% for i = 2:(D_size(1)-1)
%     for j = 1:D_size(2)
%         if pic_org_sw_clean(i-1,j,1) < pic_org_sw_clean(i,j,1) && pic_org_sw_clean(i+1,j,1)== pic_org_sw_clean(i,j,1) 
%             pic_org_sw_clean(i-1,j,1) = 255;
%             pic_org_sw_clean(i-1,j,2) = 255;
%             pic_org_sw_clean(i-1,j,2) = 255;
%         end
%     end
% end
% 
% figure
% image(pic_org_sw_clean)
% axis equal