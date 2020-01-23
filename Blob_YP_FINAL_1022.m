%%%%%%%%%%%%% Blob Maker modified my Yunji_Park on 8/12/2017 %%%%%%%%%%%%%%
clear
close all
%rng shuffle

%% set parameters

Area_1=[33	22	66	60	68	16	20	24	48	16	84	19	23	25	39	60	92	14	21	27	21	38	42	23	19	17	24	20	24	23	18	22	37	64	40	17	12	26	30	18	54	22	15	29	26	72	62]; % orange area
Area_2=[110	40	120	75	85	60	75	45	90	20	105	79	96	48	75	75	115	70	105	90	70	95	105	115	95	73	103	75	90	115	90	60	101	120	75	85	60	65	75	30	90	110	75	67	60	108	93]; % blue area

x_length = [14	8	13	13	11	15	11	11	15	7	16	10	12	10	14	13	17	11	13	12	14	13	13	14	12	11	13	13	12	12	14	10	13	15	13	12	10	11	11	8	14	14	11	15	12	15	14];
y_length = [10.21428571	7.75	14.30769231	10.38461538	13.90909091	5.066666667	8.636363636	6.272727273	9.2	5.142857143	11.8125	9.8	9.916666667	7.3	8.142857143	10.38461538	12.17647059	7.636363636	9.692307692	9.75	6.5	10.23076923	11.30769231	9.857142857	9.5	8.181818182	9.769230769	7.307692308	9.5	11.5	7.714285714	8.2	10.61538462	12.26666667	8.846153846	8.5	7.2	8.272727273	9.545454545	6	10.28571429	9.428571429	8.181818182	6.4	7.166666667	12	11.07142857];


%Area_1=randi(42); % white area
%Area_2=randi(168); 

fnames = {'3_2_P1_33_110'	'3_2_P2_22_40'	'3_2_P2_66_120'	'3_2_P3_60_75'	'3_2_P3_68_85'	'4_3_P1_16_60'	'4_3_P1_20_75'	'4_3_P2_24_45'	'4_3_P2_48_90'	'4_3_P3_16_20'	'4_3_P3_84_105'	'6_5_P1_19_79'	'6_5_P1_23_96'	'6_5_P2_25_48'	'6_5_P2_39_75'	'6_5_P3_60_75'	'6_5_P3_92_115'	'2_1_P1_14_70'	'2_1_P1_21_105'	'2_1_P2_27_90'	'2_1_P2_21_70'	'2_1_P3_38_95'	'2_1_P3_42_105'	'3_1_P1_23_115'	'3_1_P1_19_95'	'3_1_P2_17_73'	'3_1_P2_24_103'	'3_1_P3_20_75'	'3_1_P3_24_90'	'3_2_P1_23_115'	'3_2_P1_18_90'	'3_2_P2_22_60'	'3_2_P2_37_101'	'3_2_P3_64_120'	'3_2_P3_40_75'	'4_3_P1_17_85'	'4_3_P1_12_60'	'4_3_P2_26_65'	'4_3_P2_30_75'	'4_3_P3_18_30'	'4_3_P3_54_90'	'6_5_P1_22_110'	'6_5_P1_15_75'	'6_5_P2_29_67'	'6_5_P2_26_60'	'6_5_P3_72_108'	'6_5_P3_62_93'};
char_fnames=char(fnames);

resize_x=400;
resize_y=400;


%% main for loop
f=figure('units','normalized','outerposition',[0.2 0.2 0.4 0.6]);
for Ai=1:size(Area_1,2)
    R=Area_1(Ai)/Area_2(Ai);
    
    
    %% create base ellipse
    %Features added by Mohammad Akanda 4/21/2017
    % Major and minor axis of basic ellipse
    a = x_length(Ai);
    b = y_length(Ai);
    % Basis points
    n = 18;dt = 2*pi/n;t = [0:dt:2*pi-dt];
    % Generates a perturbed ellipse
    pr = rand(1,n)*(0.5);
    pt = -dt/3 + 2*dt/3*rand(1,n)*(0.6);
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
    
    all_ratios=zeros(size(ww1,1)-20,1);
    c=0;
    for i=2:size(all_ratios,1)
        t1 = 1;
        t2 = i;
        a1 = [ww1(t1,1) ww1(t2,1)];
        b1 = [ww1(t1,2) ww1(t2,2)];
        aa = linspace(a1(1),a1(2),7);
        bb = linspace(b1(1),b1(2),7);
        %noise generation to line
        ff = 4;
        %ee = abs(bb(1)-bb(end))/ff;
        ee=rand(1)+1;
        r = 2*ee.*rand(1,3) -ee;
        r = [0 0 r 0 0];
        bb1 = bb +r;
        aa1 = linspace(a1(1),a1(2),70);
        bb1 = interp1(aa, bb1, aa1, 'spline');
        %first segment of purturned ellipse
        ee1 = ww1(t1:t2,:);
        dd = [aa1' bb1'];
        tt1 = dd(1:end-1,:);
        tt2 = dd(2:end,:);
        con1 = [ee1;flipud(tt1)];
        A1 = polyarea(con1(:,1),con1(:,2));
        %second segment of purturned ellipse
        ee2 = [ww1(t2:end,:);ww1(2:t1,:)];
        ef1 = ww1(t2+1:end,:);
        con2 = [dd;ef1 ;ww1(2:t1,:)];
        A2 = polyarea(con2(:,1),con2(:,2));
        %tot = A1+A2;
        c=c+1;
        all_ratios(c,1) = A1/A2;
        save(['dummy_file/m' num2str(i)], 'con1','con2','ee1' )
        
    end
    
    %% finding the best ratio
    all_ratios1=abs(all_ratios-R);
    [q1,q2]=min(all_ratios1);
    
    %% plot best
    load(['dummy_file/m' num2str(q2)])
   
    
    A1_final = polyarea(con1(:,1),con1(:,2));
    A2_final = polyarea(con2(:,1),con2(:,2));
    Ratio_final=A1_final/A2_final;
    
   

    fill(con1(:,1),con1(:,2),[1,1,1])
    hold on
    fill(con2(:,1),con2(:,2),[0,0,0])
    plot(ee1(:,1),ee1(:,2),'black')
    set(gcf,'Color',[192/255 192/255 192/255]);
    
    axis([-23 23 -23 23]);
    axis off;
    hold off
    drawnow;
    frm = getframe( f );
     A = imresize(frm.cdata,[resize_y resize_x]);

    imwrite(A,['Pilot/' char_fnames(Ai,:) '.png'] ); 
    
    
end % main forloop

