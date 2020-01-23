%%%%%%%%%%%%% Blob Maker modified my Yunji_Park on 8/12/2017 %%%%%%%%%%%%%%
clear
close all
%rng shuffle

%% set parameters

Area_1=[100]; % orange area
Area_2=[200]; % blue area

x_length = [13];
y_length = [5];


%Area_1=randi(42); % white area
%Area_2=randi(168); 

fnames = {'a'};
char_fnames=char(fnames);


%% main for loop
f=figure('units','normalized','outerposition',[0.2 0.2 0.4 0.6]);
for Ai=1:size(Area_1,2)
    %R=Area_1(Ai)/Area_2(Ai);
    
    
    %% create base ellipse
    %Features added by Mohammad Akanda 4/21/2017
    % Major and minor axis of basic ellipse
    a = x_length(Ai);
    b = y_length(Ai);
    % Basis points
    n = 18;dt = 2*pi/n;t = [0:dt:2*pi-dt];
    % Generates a perturbed ellipse
    pr = rand(1,n)*(0.2);
    pt = -dt/3 + 2*dt/3*rand(1,n)*(0.2);
    x = (2*a/3 + 2*a/3*pr).*cos(t + pt);
    y = (2*b/2 + 2*b/3*pr).*sin(t + pt);
    % Interpolate (using (r, t) pairs)
    t = t + pt;
    r = sqrt(x.^2 + y.^2);
    ll =800;%number of points consider
    % Smoothness conformity
    tq = linspace(t(1), 2*pi, ll);
    t = [t-2*pi t t+2*pi];
    r = [r r r];
    rq = interp1(t, r, tq, 'spline');
    xx = rq.*cos(tq); yy = rq.*sin(tq);
    ww = [xx' yy'];
    ww1  = ww(1:ll-10,:);
    ww2 = ww(ll-9:ll,:);
    for jj = 1:10
        if ww1(1,2)>ww2(jj,2)
            ww1 = [ww1;ww2(jj,:)];
        else
        end
    end
    ww1 = [ww1;ww1(1,:)];
    ww1(1:5,:);
    ww1(end-4:end,:);
    
    %% find all ratios
%     
%     all_ratios=zeros(size(ww1,1)-20,1);
%     c=0;
%     for i=2:size(all_ratios,1)
%         t1 = 1;
%         t2 = i;
%         a1 = [ww1(t1,1) ww1(t2,1)];
%         b1 = [ww1(t1,2) ww1(t2,2)];
%         aa = linspace(a1(1),a1(2),7);
%         bb = linspace(b1(1),b1(2),7);
%         %noise generation to line
%         ff = 4;
%         %ee = abs(bb(1)-bb(end))/ff;
%         ee=rand(1)+1;
%         r = 2*ee.*rand(1,3) -ee;
%         r = [0 0 r 0 0];
%         bb1 = bb +r;
%         aa1 = linspace(a1(1),a1(2),70);
%         bb1 = interp1(aa, bb1, aa1, 'spline');
%         %first segment of purturned ellipse
%         ee1 = ww1(t1:t2,:);
%         dd = [aa1' bb1'];
%         tt1 = dd(1:end-1,:);
%         tt2 = dd(2:end,:);
%         con1 = [ee1;flipud(tt1)];
%         A1 = polyarea(con1(:,1),con1(:,2));
%         %second segment of purturned ellipse
%         ee2 = [ww1(t2:end,:);ww1(2:t1,:)];
%         ef1 = ww1(t2+1:end,:);
%         con2 = [dd;ef1 ;ww1(2:t1,:)];
%         A2 = polyarea(con2(:,1),con2(:,2));
%         %tot = A1+A2;
%         c=c+1;
%         all_ratios(c,1) = A1/A2;
%         save(['dummy_file/m' num2str(i)], 'con1','con2','ee1' )
%         
%     end
    
%     %% finding the best ratio
%     all_ratios1=abs(all_ratios-R);
%     [q1,q2]=min(all_ratios1);
    
    %% plot best
%     load(['dummy_file/m' num2str(q2)])
%    
%     
%     A1_final = polyarea(con1(:,1),con1(:,2));
%     A2_final = polyarea(con2(:,1),con2(:,2));
%     Ratio_final=A1_final/A2_final;
%     
   

 %   fill(con1(:,1),con1(:,2),[252/255, 127/255, 36/255])
 %   hold on
 %   fill(con2(:,1),con2(:,2),[38/255, 121/255, 178/255])
 %   plot(ee1(:,1),ee1(:,2),'b')
  

figure
fill(ww1(:,1),ww1(:,2),'black')

hold on
plot(ww1(:,1),ww1(:,2),'black')

 set(gcf,'Color',[192/255 192/255 192/255]);
    
    axis([-30 30 -30 30]);
    axis off;
    hold off
    drawnow;
    frm = getframe( f );
    imwrite( frm.cdata,['Pilot/' char_fnames(Ai,:) '.png'] ); 
    
    
end % main forloop

