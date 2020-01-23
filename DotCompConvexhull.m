close all
clear
%% parameter setting
top_dot_num=[12];
bot_dot_num=[60];

top_dot_area=[2400];
bot_dot_area=[2400];

fnames = {'2_1_P1_R3_12_60'};
%char_fnames=char(fnames);

xsize = 700; % box size
ysize = 500;

bot_grid_size=0.4; % grid size varies between 0.1~0.9
bot_bound_low=0.09;
bot_bound_high=0.91;
bot_rad_size=0.3; % range 0.1~0.6

top_grid_size=0.4; % grid size varies between 0.1~0.9
top_bound_low=0.13;
top_bound_high=0.87;




code_type=1; % 1) area_controlled condition; 2) size_controlled condition

%% presettings
number_of_samples=size(top_dot_num,2);
dot_max_num=120;
box_size=150;
box_area=box_size^2;

gap1_bot=150; %dot_grid
gap2_bot=150; %box_grid

gap1_top=150;
gap2_top=150; 

f=figure;
img = imread('baseim.jpg');
img = imresize(img, [ysize  xsize]);
imshow( img, 'border', 'tight' );

hold on;
color=192/255;
rectangle('Position', [0,0, xsize,ysize],'FaceColor',[color color color],'EdgeColor',[color color color]);

%rectangle('Position', [100,100, 200, 400],'FaceColor',[color color color],'EdgeColor',[color color color]);
%% main: all dots should be equally distributed in a given rectangle.

if code_type==1
    %% 1) area controlled condition: sum of the area of dots should be equal for each array and dot sizes vary (e.g. 80%~120%)
 
     
    
    for sample_i=1:number_of_samples
        % BOTTOM
        area_per_dot=bot_dot_area(sample_i)/bot_dot_num(sample_i);
        %rad_per_dot=sqrt(area_per_dot/pi);
        grid_num=ceil(sqrt(bot_dot_num(sample_i)));
        grid_delta=(bot_bound_high-bot_bound_low)*box_size/(grid_num-1);
        [bot_gridx, bot_gridy]=meshgrid( box_size*bot_bound_low: grid_delta:box_size*bot_bound_high,box_size*bot_bound_low: grid_delta:box_size*bot_bound_high);
        grid_tot=size(bot_gridx,1)*size(bot_gridx,1);
        
        ar=1:grid_tot;
        vert=[1,size(bot_gridx,2),size(bot_gridx,2)*size(bot_gridx,1)-size(bot_gridx,1)+1,size(bot_gridx,2)*size(bot_gridx,1)];
        
        grid_array=ar(randperm(grid_tot));
        for j=1:4
            a=sum((1:size(grid_array,2)).*(grid_array==vert(j)));
            b=grid_array(j);
            grid_array(j)=vert(j);
            grid_array(a)=b;
        end
        grid_array=grid_array(1:bot_dot_num(sample_i));
        
        bot_grid=[reshape(bot_gridx,grid_tot,1),reshape(bot_gridy,grid_tot,1)];
        bot_grid=bot_grid(grid_array,:);
        
        %bot_grid=bot_grid+bot_grid_size*(rand(bot_dot_num(sample_i),2)*grid_delta-grid_delta/2); %%dots
        
        % bot_grid=[bot_grid(:,1)+box_size/2, bot_grid(:,2)+100+box_size];
        
        
        bot_grid=[bot_grid(:,1)+xsize/2-box_size/2+gap1_bot, bot_grid(:,2)+ysize/2-box_size/2];%
        
        pair_dist=pdist2(bot_grid,bot_grid);
        nonzero_dist=pair_dist(pair_dist>0);
        max_rad=bot_rad_size*min(nonzero_dist);
        %random_rad=(rand(bot_dot_num(sample_i),1)*0.4+0.8)*max_rad;
        
        random_area=area_per_dot*(rand(bot_dot_num(sample_i),1)*0.4+0.8);
        random_rad=sqrt(random_area./pi);
        bot_grid=bot_grid+bot_grid_size*([zeros(4,2);rand(bot_dot_num(sample_i)-4,2)]*grid_delta-grid_delta/2);
        bot_grid(2,2)=bot_grid(2,2)+grid_delta/3;
        bot_grid(3,1)=bot_grid(3,1)+grid_delta/3;
        bot_grid(4,1)=bot_grid(4,1)+grid_delta/3;
        bot_grid(4,2)=bot_grid(4,2)+grid_delta/3;
        
        
        color=[38/255, 121/255, 178/255];
        %box_start1=xsize/2-box_size/2;
        rectangle('Position', [xsize/2-box_size/2+gap2_bot,ysize/2-box_size/2, box_size, box_size],'FaceColor',color,'EdgeColor',color);
        for i=1:bot_dot_num(sample_i)
            
            pos_bot_dot=[bot_grid(i,1)-random_rad(i), bot_grid(i,2)-random_rad(i), random_rad(i)*2,random_rad(i)*2];
            color_dot=[252/255, 127/255, 36/255];
            rectangle('Position', pos_bot_dot,'Curvature',[1 1],'FaceColor',color_dot,'EdgeColor',color_dot);
        end
        %bot_dot_area=(random_rad.^2)*pi;
        %bot_area=sum(bot_dot_area);
        % TOP
        area_per_dot=top_dot_area(sample_i)/top_dot_num(sample_i);
        grid_num=ceil(sqrt(top_dot_num(sample_i)));
        grid_delta=(top_bound_high-top_bound_low)*box_size/(grid_num-1);
        [top_gridx, top_gridy]=meshgrid( box_size*top_bound_low: grid_delta:box_size*top_bound_high,box_size*top_bound_low: grid_delta:box_size*top_bound_high);
        grid_tot=size(top_gridx,1)*size(top_gridx,1);
        ar=1:grid_tot;
        vert=[1,size(top_gridx,2),size(top_gridx,2)*size(top_gridx,1)-size(top_gridx,1)+1,size(top_gridx,2)*size(top_gridx,1)];
        
        grid_array=ar(randperm(grid_tot));
        for j=1:4
            a=sum((1:size(grid_array,2)).*(grid_array==vert(j)));
            b=grid_array(j);
            grid_array(j)=vert(j);
            grid_array(a)=b;
        end
        grid_array=grid_array(1:top_dot_num(sample_i));
        
        top_grid=[reshape(top_gridx,grid_tot,1),reshape(top_gridy,grid_tot,1)];
        top_grid=top_grid(grid_array,:);
        
        top_grid=[top_grid(:,1)+xsize/2-box_size/2-gap1_top, top_grid(:,2)+ysize/2-box_size/2];
        top_grid=top_grid+top_grid_size*([zeros(4,2);rand(top_dot_num(sample_i)-4,2)]*grid_delta-grid_delta/2);
        
        random_area=area_per_dot*(rand(top_dot_num(sample_i),1)*0.4+0.8);
        random_rad=sqrt(random_area./pi);
        
        top_grid(2,2)=top_grid(2,2)+grid_delta/3;
        top_grid(3,1)=top_grid(3,1)+grid_delta/3;
        top_grid(4,1)=top_grid(4,1)+grid_delta/3;
        top_grid(4,2)=top_grid(4,2)+grid_delta/3;
        
        color=[252/255, 127/255, 36/255]; % orange
        rectangle('Position', [xsize/2-box_size/2-gap2_top ,ysize/2-box_size/2, box_size, box_size],'FaceColor',color,'EdgeColor',color); %rectangle location
        for i=1:top_dot_num(sample_i)
            
            pos_top_dot=[top_grid(i,1)-random_rad(i), top_grid(i,2)-random_rad(i), random_rad(i)*2,random_rad(i)*2];
            color_dot=[38/255, 121/255, 178/255]; %blue
            rectangle('Position', pos_top_dot,'Curvature',[1 1],'FaceColor',color_dot,'EdgeColor',color_dot);
        end
        
        drawnow;
        char_len=size(fnames{1,sample_i});
        frm = getframe( f );
        
        imwrite( frm.cdata,['Pilot/'  fnames{sample_i} '.png'] );
      
        %imwritesize(frm.cdata,['R1_stimuli/' char_fnames(sample_i,:) '.tif'], 300);
        
        
    end % sample loop
    
    
    
elseif code_type==2
    %% 2) size controlled condition: equal dot size (but it is better if dot size can be manipulated within a certain range.
    for sample_i=1:number_of_samples
        grid_num=ceil(sqrt(bot_dot_num(sample_i)));
        grid_delta=(bot_bound_high-bot_bound_low)*box_size/(grid_num-1);
        [bot_gridx, bot_gridy]=meshgrid( box_size*bot_bound_low: grid_delta:box_size*bot_bound_high,box_size*bot_bound_low: grid_delta:box_size*bot_bound_high);
        grid_tot=size(bot_gridx,1)*size(bot_gridx,1);
        ar=1:grid_tot;
        
        
        grid_array=ar(randperm(grid_tot));
        grid_array=grid_array(1:bot_dot_num(sample_i));
        
        bot_grid=[reshape(bot_gridx,grid_tot,1),reshape(bot_gridy,grid_tot,1)];
        bot_grid=bot_grid(grid_array,:);
        
        bot_grid=bot_grid+bot_grid_size*(rand(bot_dot_num(sample_i),2)*grid_delta-grid_delta/2);
        %bot_grid=[bot_grid(:,1)+box_size/2, bot_grid(:,2)+100+box_size]; %same
        bot_grid=[bot_grid(:,1)+xsize/2-box_size/2, bot_grid(:,2)+ysize/2];
        pair_dist=pdist2(bot_grid,bot_grid);
        nonzero_dist=pair_dist(pair_dist>0);
        max_rad1=bot_rad_size*min(nonzero_dist);
        
        
        
        %TOP
        grid_num=ceil(sqrt(top_dot_num(sample_i)));
        grid_delta=(top_bound_high-top_bound_low)*box_size/(grid_num-1);
        [top_gridx, top_gridy]=meshgrid( box_size*top_bound_low: grid_delta:box_size*top_bound_high,box_size*top_bound_low: grid_delta:box_size*top_bound_high);
        grid_tot=size(top_gridx,1)*size(top_gridx,1);
        ar=1:grid_tot;
        
        grid_array=ar(randperm(grid_tot));
        grid_array=grid_array(1:top_dot_num(sample_i));
        
        top_grid=[reshape(top_gridx,grid_tot,1),reshape(top_gridy,grid_tot,1)];
        top_grid=top_grid(grid_array,:);
        
        %top_grid=[top_grid(:,1)+box_size/2, top_grid(:,2)+box_size/2];
        top_grid=[top_grid(:,1)+xsize/2-box_size/2, top_grid(:,2)+ysize/2-box_size];
        top_grid=top_grid+top_grid_size*(rand(top_dot_num(sample_i),2)*grid_delta-grid_delta/2); %same
        
        pair_dist=pdist2(top_grid,top_grid);
        nonzero_dist=pair_dist(pair_dist>0);
        max_rad2=bot_rad_size*min(nonzero_dist);
        max_rad=min(max_rad1,max_rad2);
        
        
        random_rad1=(rand(bot_dot_num(sample_i),1)*0.4+0.8)*max_rad;
        random_rad2=(rand(top_dot_num(sample_i),1)*0.4+0.8)*max_rad;
        %random_rad1=ones(bot_dot_num(sample_i),1)*max_rad;
        %random_rad2=ones(top_dot_num(sample_i),1)*max_rad;
        
        
        
        % plot
        color=[38/255, 121/255, 178/255];
        rectangle('Position', [xsize/2-box_size/2,ysize/2, box_size, box_size],'FaceColor',color,'EdgeColor',color);
        for i=1:bot_dot_num(sample_i)
            
            pos_bot_dot=[bot_grid(i,1)-random_rad1(i), bot_grid(i,2)-random_rad1(i), random_rad1(i)*2,random_rad1(i)*2];
            color_dot=[0 0 0];
            rectangle('Position', pos_bot_dot,'Curvature',[1 1],'FaceColor',color_dot,'EdgeColor',color_dot);
            
        end
        
        color=[252/255, 127/255, 36/255];
        rectangle('Position', [xsize/2-box_size/2,ysize/2-box_size, box_size, box_size],'FaceColor',color,'EdgeColor',color); %rectangle location
        for i=1:top_dot_num(sample_i)
            
            pos_top_dot=[top_grid(i,1)-random_rad2(i), top_grid(i,2)-random_rad2(i), random_rad2(i)*2,random_rad2(i)*2];
            color_dot=[1 1 1];
            rectangle('Position', pos_top_dot,'Curvature',[1 1],'FaceColor',color_dot,'EdgeColor',color_dot);
            
        end
        drawnow;
        
        frm = getframe( f );
        imwrite( frm.cdata,['Practice/' fnames{sample_i} '.png'] );
        
        bot_dot_area=(random_rad1.^2)*pi;
        area_bot=mean(bot_dot_area);
        top_dot_area=(random_rad2.^2)*pi;
        area_top=mean(top_dot_area);
    end % sample_i
    
    
end



