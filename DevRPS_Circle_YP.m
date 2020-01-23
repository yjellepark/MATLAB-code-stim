close all
clear 
%% parameter setting
n1_dia = [236 35];
% pixel size - area 
d1_dia = [250 42];
fnames = { 'a' 'b'};
char_fnames=char(fnames);


%% main
f=figure;
number_of_samples=size(n1_dia,2);
for n = 1:number_of_samples

n1_rad = (n1_dia(n)/2)*2;
d1_rad = (d1_dia(n)/2)*2;
   
% %n1_area = 2000;
% d1_area = n1_area*(1/R1(n));
% n2_area = 2000;
% d2_area = n2_area*(1/R2(n));

%n1_rad = (n1_area(n))^(1/2);
%d1_rad = (d1_area(n))^(1/2);
% n2_rad = (n2_area/pi)^(1/2);
% d2_rad = (d2_area/pi)^(1/2);

gap_h = 10;
gap_p = 160;

xsize = 600;
ysize = 1048;

img = imread('baseim.jpg');
img = imresize(img, [ysize  xsize]);

imshow( img, 'border', 'tight' ); %//show your image
hold on;
color=0;
rectangle('Position', [0,0, xsize, ysize],'FaceColor',[color color color],'EdgeColor',[color color color]);


Y = 0;

%n1 circle
pos_n1 = [xsize*(1/2)-(n1_rad),   ysize*(50/100)-(gap_h/2+2*n1_rad)-Y,   2*n1_rad,  2*n1_rad]; 
col_n = [200/255,200/255,200/255];
rectangle('Position', pos_n1,'Curvature',[1 1],'FaceColor',col_n,'EdgeColor',col_n);



%d1 circle
pos_d1 = [xsize*(1/2)-(d1_rad),   ysize*(50/100)+(gap_h/2)-Y,   2*d1_rad,  2*d1_rad]; 
col_d= [200/255,200/255,200/255];  
rectangle('Position', pos_d1,'Curvature',[1 1],'FaceColor',col_d,'EdgeColor',col_d);


% 
% %n2 circle
% pos_n2 = [xsize*(3/4)-(n2_rad),   ysize/2-(gap_h/2+2*n2_rad),   2*n2_rad,  2*n2_rad]; 
% col_n = [1,1,1];
% rectangle('Position', pos_n2,'Curvature',[1 1],'FaceColor',col_n,'EdgeColor',col_n);
% 
% 
% 
% %d2 circle
% pos_d2 = [xsize*(3/4)-(d2_rad),   ysize/2+(gap_h/2),   2*d2_rad,  2*d2_rad];  %%wierd position
% col_d= [0,0,0];   
% rectangle('Position', pos_d2,'Curvature',[1 1],'FaceColor',col_d,'EdgeColor',col_d);

%rectangle('Position', [100 200 20 10],'Curvature',[1 1],'FaceColor',col_d,'EdgeColor',col_d);
hold off
drawnow;

frm = getframe( f );
imwrite( frm.cdata,['tryout/' fnames{n} '.png'] ); 



end
