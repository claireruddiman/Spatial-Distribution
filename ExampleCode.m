% 
% This code analyzes spatial distribution of holes in the internal elastic
% lamina (HIEL) and myoendothelial junction (MEJ) localized lipids/proteins
% in the following sections:  
% 
% 
% % %  SECTION ONE
% 
% % %  Thresholds and identifies holes in the internal elastic lamina
% % %  (HIEL). These are objects void of fluorescent signal.
% 
% 
% % %  SECTION TWO
% 
% % %  Thresholds and identifies a fluorescent signal that
% % %  represents an protein of interest that will be used to evaluate the
% % %  spatial pattern of the HIEL. In this example, claudin 5 (cldn5), is
% % %  the fluoresecent signal and will be referred to throughout.
% 
% 
% % %  SECTION THREE
% 
% % %  This is a quality control section that removes objects identified as
% % %  HIEL that do not meet the circularity requirement. In this section,
% % %  there is also the option to manually add back in HIEL removed by 
% % %  circularity, and to remove additional nonspecific HIEL that may have
% % %  been detected.
% 
% 
% % %  SECTION FOUR
% 
% % %  This section runs simulations to evaluate the spatial pattern of
% % %  HIEL with respect to the fluorescent signal that was thresholded in
% % %  section two. This will refer to 3 separate Matlab function files.
% 
% 
% % %  SECTION FIVE
% 
% % %  Thresholds and identifies high intensity fluoresecent puncta that
% % %  overlaps in HIEL. This is useful to determine the extent of protein
% % %  localization in the HIEL. In this example phosphatidylserine, or PS,
% % %  may be referred to throughout.
% 
% 
% % %  SECTION SIX
% 
% % %  This section runs simulations to evaluate the spatial pattern of
% % %  fluorescent puncta in HIEL with respect to the fluorescent signal 
% % %  that was thresholded in section two. This will refer to 3 separate 
% % %  Matlab function files.
% 
% 
% % %  SECTION SEVEN
% 
% % %  Saves all matrices and variables that will be used for analysis.
% 
% 
% 
% 
% % 
% % % % % 
% % % % % % % % %  SECTION ONE
% % % % % 
% % 
% 
% 
% Number of iterations that the Monte Carlo, Positive, and Negative control
% simulations will be run. 
sim1=1; % For HIEL (section1) to fluorescent signal (section2)
sim2 =1; % For puncta-HIEL (section 5) to fluorescent signal (section2)
% Keep simulations at 1 iteration while troubleshooting code. 
% Recommended to increase to 1000 simulations for final data.
% Will need to run on a virtual supercomputing cluser.


% Define range of pixels wish to be included for analysis.
lowlim1_X = 0; %adjust all based on image/ROI
uplim1_X = 1024;
lowlim1_Y = 0; 
uplim1_Y = 1024; 
Pixel_min = 0; 
Pixel_max = 1024;
% In this example, all pixels of a 1024x1024 image are considered.


% Define RGB file. Must be an RGB image (can do that in ImageJ).
RGB = imread('M9V7Z1_hydr.png');

% Convert to grayscale.
I=rgb2gray(RGB); 

% Blur the image.
P = fspecial('motion',5,40); %adjust to optimize
% Step 2 of blurring.
I2 = imfilter(I,P,'replicate'); 

% Idenitfy puncta of interest by converting the blurred image to binary.
bw =imbinarize(I2,'adaptive','ForegroundPolarity','dark','Sensitivity',0.1); %adjust to optimize
% In this example, the 'puncta' of interest are actually areas void of fluorescent signal. 
% For this reason, 'dark' option is selected. For identificaiton of
% fluorescent puncta, you do not need the specifications of
% 'ForegroundPolarity' and 'dark'.


% Refine the identification of HIEL.
SE = strel('disk',5); % ajust to optimize
I3 = imerode(bw,SE);


% % FIGURE 1 
% % Show boundaries of ROI on top of the originial image. 
% % Visual confirmation that HIEL are being correctly identified.

Figure1 =figure; 
[B_M9V7Z1A1,L] = bwboundaries(I3); %detected outline of ROIs.
imshow(I3); hold on; %show detected outline
imshow(I); hold on;  %show original image
colors=('c');
m = size(B_M9V7Z1A1,1);
n = size (B_M9V7Z1A1,2);
area = zeros(m,n); 
stats = zeros(m,n);
stats = regionprops(L,'Area');

% Plot boundaries of ROIs. 
for k = 1:length(B_M9V7Z1A1)
  boundary = B_M9V7Z1A1{k};
    cidx = mod(k,length(colors))+1;
    plot(boundary(:,2), boundary(:,1),...
       colors(cidx),'LineWidth',2);

      area(k) = stats(k).Area;
% % %  actual number of pixels in the area    


% %   add text to show object numbers 
  rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
  col = boundary(rndRow,2); row = boundary(rndRow,1);
  h = text(col+1, row-1, num2str(L(row,col)));
  set(h,'Color',colors(cidx),'FontSize',12,'FontWeight','bold');
% %   randomize text position for better visibility

end

area = transpose(area);


% % Convert cell array of border pixels and plot centers of each object
Q = cellfun(@mean,B_M9V7Z1A1, 'UniformOutput', false); 
% get the average of the pixels along the border to calculate center of ROI
% NOTE that both Q and b are cell arrays


H = cell2mat(Q);
% % % H is a matrix, which is needed to plot the various coordinates
 
X = H(:,2:2:end);
Y = H(:,1:2:end);

% % % plots the X,Y coordinates of H (centers of the IEL holes)
plot(X,Y,'.','MarkerSize',12, 'Color', 'm');

% Draw vertical and horizontal lines on the image to see the image area you 
% are analyzing.  This is only relevant if you want to analyze a small area 
% of the original image.   
xline(lowlim1_X);hold on;
xline(uplim1_X);hold on;
yline(lowlim1_Y);hold on;
yline(uplim1_Y);hold on;
saveas(Figure1, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure1.png');
hold off;

% Round the final coordinates of HIEL center detection and store in matrix.
coord_M9V7Z1A1 = horzcat(round(X),round(Y));
% 
% 
% 
% % 
% % % % % 
% % % % % % % % %  SECTION TWO
% % % % % 
% % 
% 
% 
% 
% Threshold fluorescent signal claudin 5 (cldn5) EC borders
% (interendothelial junctions).


% Define RGB file. Must be an RGB image (can do that in ImageJ).
% This is fluorescent signal that will be used to define HIEL spatial
% pattern.
RGB_ = imread('M9V7Z1_cldn5_1.png'); 

% Convert to grayscale.
I_=rgb2gray(RGB_); 

% Blur the image.
P_ = fspecial('motion',5,30); 
% Step 2 of blurring.
I2_ = imfilter(I_,P_,'replicate');

% Idenitfy puncta of interest by converting the blurred image to binary.
bw_ =imbinarize(I2_,'adaptive','Sensitivity',0.1);%adjust to optimize

% Refine the identification of fluorescent signal.
SE_ = strel('rectangle',[4,1]); % % % pixel size of disk that helps to define the IEL holes
I4_ = imerode(bw_,SE_);
se_ = strel('rectangle',[10,1]); 
I5_ = imclose(I4_,se_);
I3_ = bwareaopen(I5_,100);

% Plot the result of thresholding the fluorescent signal.
figure 
imshow(I3_)
[rows,cols] = find(I3_); % Extract XY coordinates from I3_, which is the claudin5 final binary image
claudin5_full_M9V7Z1A1 = horzcat(cols,rows);

% Removes any fluoresecent signal from analysis that is outside of the
% desired range.
claudin5_full_M9V7Z1A1(claudin5_full_M9V7Z1A1(:, 1) < lowlim1_X & (claudin5_full_M9V7Z1A1(:, 1) >= Pixel_min), :) = [];
claudin5_full_M9V7Z1A1(claudin5_full_M9V7Z1A1(:, 1) <= Pixel_max & (claudin5_full_M9V7Z1A1(:, 1) > uplim1_X), :) = [];
claudin5_full_M9V7Z1A1(claudin5_full_M9V7Z1A1(:, 2) < lowlim1_Y & (claudin5_full_M9V7Z1A1(:, 2) >= Pixel_min), :) = [];
claudin5_full_M9V7Z1A1(claudin5_full_M9V7Z1A1(:, 2) <= Pixel_max & (claudin5_full_M9V7Z1A1(:, 2) > uplim1_Y), :) = [];

% Define sizes of HIEL and fluorescent signal matrices.
[mmm nnn] = size(coord_M9V7Z1A1);
nnn = mmm;
[mmmm nnnn] = size(claudin5_full_M9V7Z1A1);
nnnn = mmmm;
distance_ = zeros(nnn,mmmm);
% 
% 
% 
% % 
% % % % % 
% % % % % % % % %  SECTION THREE
% % % % % 
% % 
% 
% 
% Remove objects which do not correspond to IEL holes through a
% circularity prediction.


% % FIGURE 2 
colors2=('m');
Figure2 = figure;

imshow(I3); hold on;
imshow(I); hold on;

% Display vertical lines to indicate analyzed area.
xline(lowlim1_X);hold on;
xline(uplim1_X);hold on;
yline(lowlim1_Y);hold on;
yline(uplim1_Y);hold on;


% Begin loop to check circularity of detected HIEL
for k = 1:length(B_M9V7Z1A1)
    boundary = B_M9V7Z1A1{k};
    cidx = mod(k,length(colors))+1;
    plot(boundary(:,2), boundary(:,1),...
       colors(cidx),'LineWidth',2);
    
    % circumference of an HIEL, identified by the number of pixels. not as accurate as c_f2. see below.
    c(k)=size(B_M9V7Z1A1{k},1); 
    
    % HIEL area
    area(k) = stats(k).Area;
    
    % calculate circularity
    metric(k) = (4*pi*area(k))/c(k)^2 ; 
    
    % HIEL area
    AA(k) = stats(k).Area; 
    
    % HIEL radius
    r(k) = sqrt(AA(k)/pi); 
    
    % HIEL diameter
    d(k) = 2*r(k); 
    
    % HIEL circumference, as calculated from area
    c_f2(k) = 2*pi*r(k); 
    label(k) = k;
           for j = 1 : mmmm % size of fluorescent signal (cldn5) matrix
                scale_M9V7Z1A1 = 8.6719; % this is the number of pixels per micron
                distance_(k, j) = (sqrt((coord_M9V7Z1A1(k, 1) - claudin5_full_M9V7Z1A1(j, 1)) ^ 2 + ...
            (coord_M9V7Z1A1(k, 2) - claudin5_full_M9V7Z1A1(j, 2)) ^ 2))/scale_M9V7Z1A1;
           end     
           
   % if circularity is greater than 0.5 then keep the predicted object
    if metric(k) > 0.50 
        G{k} = B_M9V7Z1A1{k};
        
        % HIEL area
        AA(k) = stats(k).Area;
        
        % HIEL radius
        r(k) = sqrt(AA(k)/pi); 
        
        % HIEL diameter
        d(k) = 2*r(k);
        
        % HIEL circumference as calculated from area. this is the most
        % accurate. Use this value for downstream analysis
        c_f2(k) = 2*pi*r(k); 
        
        label(k) = k;
        X(k) = H(k,2:2:end);
        Y(k) = H(k,1:2:end);
            for j = 1 : mmmm % size of fluorescent signal (cldn5) matrix
                scale_M9V7Z1A1 = 8.6719;  % this is the number of pixels per micron
                distance_(k, j) = (sqrt((coord_M9V7Z1A1(k, 1) - claudin5_full_M9V7Z1A1(j, 1)) ^ 2 + ...
            (coord_M9V7Z1A1(k, 2) - claudin5_full_M9V7Z1A1(j, 2)) ^ 2))/scale_M9V7Z1A1;

            end
        
        metric_string = sprintf('%2.2f',metric(k)); % % % store circularity values in a string
        text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','b',...
       'FontSize',12,'FontWeight','bold') % % % % print circularity values on plot

   % if circularity is less than 0.5 remove the object    
    else
        G{k} = [];
        AA (k) = [];
        r(k)=[]; %
        d(k) = [];%
        c_f2(k) = [];%
        
            for j = 1 : mmmm % size of fluorescent signal (cldn5) matrix
                scale_M9V7Z1A1 = 8.6719;
                distance_(k, j) = 0;
            end 
            
        label(k) = [];
        X(k) = 0;
        Y(k) = 0;
        boundary = B_M9V7Z1A1{k};
        cidx = mod(k,length(colors2))+1;
        plot(boundary(:,2), boundary(:,1),...
        colors2(cidx),'LineWidth',2);
    end

end

saveas(Figure2, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure2.png');
hold off;
% 
% 
% 
% % 
% % % % % 
% % % % % % % % % MANUALLY ADD BACK OBJECTS REMOVED BY CIRCULARITY
% (OPTIONAL)
% % % % % 
% % 
% 
% 
% 
% 
% UNCOMMENT THIS SECTION if you wish to manually add back in any of the
% objects removed by circularity.

% % %%Keep some of the values ruled out by circularity 
%  for k=[add number of objects in here];
%         G{k} = B_M9V7Z1A1{k};
%         AA(k) = stats(k).Area;
%         r(k) = sqrt(AA(k)/pi)/scale_M9V7Z1A1; %
%         d(k) = 2*r(k); %
%         c_f2(k) = 2*pi*r(k);%
%         label(k) = k;
%         X(k) = H(k,2:2:end);
%         Y(k) = H(k,1:2:end);
% 
%             for j = 1 : mmmm % size of fluorescent signal (cldn5) matrix
%                 scale_M9V7Z1A1 = 8.6719;  % this is the number of pixels per micron
%                 distance_(k, j) = (sqrt((coord_M9V7Z1A1(k, 1) - claudin5_full_M9V7Z1A1(j, 1)) ^ 2 + ...
%             (coord_M9V7Z1A1(k, 2) - claudin5_full_M9V7Z1A1(j, 2)) ^ 2))/scale_M9V7Z1A1;
%             end
% 
%         metric_string = sprintf('%2.2f',metric(k)); % % % store circularity values in a string
%         text(boundary(1,2)-35,boundary(1,1)+13,metric_string,'Color','b',...
%        'FontSize',12,'FontWeight','bold') % % % % print circularity values on plot
% 
%  end 
%    
% % 

% 
% % 
% % % % % 
% % % % % % % % % MANUALLY REMOVE ADDITIONAL NONSPECIFIC OBJECTS (OPTIONAL)
% % % % % 
% % 
% 
% 
% 
% 
G{1}=[]; % OBJECT 1 IS NONSPECIFIC.
% G{2}=[]; % ADD A LINE FOR EACH OBJECT THAT YOU WISH TO REMOVE.


F2 = (G(~cellfun('isempty',G))); % F2 is the cell array that excludes circularity values not within range
        F = transpose(F2);

distance_(1,:) = 0; % OBJECT 1 IS NONSPECIFIC.
% distance_(2,:) = 0; % ADD A LINE FOR EACH OBJECT THAT YOU WISH TO REMOVE.

          distance_(~any(distance_,2),:) = []; %delete rows where all values are zero
          distance_sort_old = transpose(sort(distance_,2)); %sorting
          distance_avgMIN_hist_M9V7Z1A1 = transpose((distance_sort_old(1,:))); %%include this to be able to plot a histogram of the minimum distances for each IEL hole to claudin 5 
          distance_avgMIN_M9V7Z1A1= mean(distance_sort_old(:,1));
%         
          Total = 100; %flexible, could also do multiple of these?
          distance_sort_M9V7Z1A1 = distance_sort_old(1:Total, 1:size(distance_sort_old,2)); % take top 1000 values from this array 


X(1) = 0; % OBJECT 1 IS NONSPECIFIC.
% X(2) = 0; % ADD A LINE FOR EACH OBJECT THAT YOU WISH TO REMOVE.


        X=transpose(X);
        X(:,~any(X,1)) = []; %rows
        X=transpose(X);

        

Y(1) = 0; % OBJECT 1 IS NONSPECIFIC.
% Y(2) = 0; % ADD A LINE FOR EACH OBJECT THAT YOU WISH TO REMOVE.


        Y=transpose(Y);
        Y(:,~any(Y,1)) = []; %rows
        Y=transpose(Y);


AA(1) = 0; % OBJECT 1 IS NONSPECIFIC.
% AA(2) = 0; % ADD A LINE FOR EACH OBJECT THAT YOU WISH TO REMOVE.


        AA( :,~any(AA,1) ) = []; %rows
        AREA_final =transpose(AA);
        

c_f2(1) = 0; % OBJECT 1 IS NONSPECIFIC.
% c_f2 (2) = 0; % ADD A LINE FOR EACH OBJECT THAT YOU WISH TO REMOVE.

        c_f2( :,~any(c_f2,1) ) = []; %rows
        C_final =transpose(c_f2);

        

r(1) = 0; % OBJECT 1 IS NONSPECIFIC.
% r (2) = 0; % ADD A LINE FOR EACH OBJECT THAT YOU WISH TO REMOVE.

        r( :,~any(r,1) ) = []; %rows
        r_final =transpose(r);


d(1) = 0; % OBJECT 1 IS NONSPECIFIC.
% d (2) = 0; % ADD A LINE FOR EACH OBJECT THAT YOU WISH TO REMOVE.
        d( :,~any(d,1) ) = []; %rows
        d_final =transpose(d);


label(1) = 0; % OBJECT 1 IS NONSPECIFIC.
% label (2) = 0; % ADD A LINE FOR EACH OBJECT THAT YOU WISH TO REMOVE.

        label( :,~any(label,1) ) = []; %rows
        label_final =transpose(label);
        

coord_M9V7Z1A1 = horzcat(round(X),round(Y));


%%Export a list of the final hole numbers
Param_M9V7Z1A1 = horzcat(label_final,AREA_final,r_final,C_final,X,Y,distance_avgMIN_hist_M9V7Z1A1);

AvgRadius_M9V7Z1A1 = mean(Param_M9V7Z1A1(:,3));
StdDevRadius = AvgRadius_M9V7Z1A1/4.2;

% Keep values only within desired analysis range. 
Param_M9V7Z1A1(Param_M9V7Z1A1(:, 5) < lowlim1_X & (Param_M9V7Z1A1(:, 5) >= Pixel_min), :) = [];
Param_M9V7Z1A1(Param_M9V7Z1A1(:, 5) <= Pixel_max & (Param_M9V7Z1A1(:, 5) > uplim1_X), :) = [];
Param_M9V7Z1A1(Param_M9V7Z1A1(:, 6) < lowlim1_Y & (Param_M9V7Z1A1(:, 6) >= Pixel_min), :) = [];
Param_M9V7Z1A1(Param_M9V7Z1A1(:, 6) <= Pixel_max & (Param_M9V7Z1A1(:, 6) > uplim1_Y), :) = [];
    
% Calculate image area in microns, used for downstream analysis. 
ImageArea_M9V7Z1A1 = ((uplim1_X-lowlim1_X)/scale_M9V7Z1A1) * ((uplim1_Y-lowlim1_Y)/scale_M9V7Z1A1) ; 
    
PartialImage_Height_M9V7Z1A1 = (uplim1_Y-lowlim1_Y);
PartialImage_Width_M9V7Z1A1 = (uplim1_X-lowlim1_X);


% Define number of HIEL detected
IELHoles_M9V7Z1A1 = length(coord_M9V7Z1A1);

% Calculate number of HIEL per image area.
HolesPerArea_M9V7Z1A1 = IELHoles_M9V7Z1A1/ImageArea_M9V7Z1A1;    
   
coord_M9V7Z1A1 = horzcat(Param_M9V7Z1A1(:,5),Param_M9V7Z1A1(:,6));

% Distance 
distance_avgMIN_hist_M9V7Z1A1 = Param_M9V7Z1A1(:,7);
distance_avgMIN_M9V7Z1A1 = mean(Param_M9V7Z1A1(:,7));
distance_sort_M9V7Z1A1 = transpose(distance_sort_M9V7Z1A1);
distance_sort_XY_M9V7Z1A1 = horzcat(X,Y,distance_sort_M9V7Z1A1);


distance_sort_XY_M9V7Z1A1(distance_sort_XY_M9V7Z1A1(:, 1) < lowlim1_X & (distance_sort_XY_M9V7Z1A1(:, 1) >= Pixel_min), :) = [];
distance_sort_XY_M9V7Z1A1(distance_sort_XY_M9V7Z1A1(:, 1) <= Pixel_max & (distance_sort_XY_M9V7Z1A1(:, 1) > uplim1_X), :) = [];
distance_sort_XY_M9V7Z1A1(distance_sort_XY_M9V7Z1A1(:, 2) < lowlim1_Y & (distance_sort_XY_M9V7Z1A1(:, 2) >= Pixel_min), :) = [];
distance_sort_XY_M9V7Z1A1(distance_sort_XY_M9V7Z1A1(:, 2) <= Pixel_max & (distance_sort_XY_M9V7Z1A1(:, 2) > uplim1_Y), :) = [];


distance_sort_XY_M9V7Z1A1(:,1) = []; %remove first column (X)
distance_sort_XY_M9V7Z1A1(:,1) = []; %comment out if using full range for Y

distance_sort_XY_M9V7Z1A1 = transpose(distance_sort_XY_M9V7Z1A1);

distance_sort_M9V7Z1A1 = distance_sort_XY_M9V7Z1A1; 

% Transpose for exporting
distance_avgMIN_hist_M9V7Z1A1 = transpose(distance_avgMIN_hist_M9V7Z1A1);


% % FIGURE 3
% Plot the boundary of HIEL and center of HIEL over top of the image.

Figure3 = figure; 
imshow(I); hold on;

for k = 1:length(F)
  boundary = F{k};
    cidx = mod(k,length(colors))+1;
    plot(boundary(:,2), boundary(:,1),...
       colors(cidx),'LineWidth',2);

% %   add text to show object numbers 
  rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
  col = boundary(rndRow,2); row = boundary(rndRow,1);
  h = text(col+1, row-1, num2str(L(row,col)));
  set(h,'Color',colors(cidx),'FontSize',12,'FontWeight','bold');
%   randomize text position for better visibility

end

% % Convert cell array of border pixels and plot centers of each object
Q = cellfun(@mean,F, 'UniformOutput', false); 
% get the average of the pixels along the border to calculate center of ROI
% NOTE that both Q and b are cell arrays

H = cell2mat(Q);
% % % H is a matrix, which is needed to plot the various coordinates

X = H(:,2:2:end);
Y = H(:,1:2:end);

plot(X,Y,'.','MarkerSize',12, 'Color', 'm'); hold on;
xline(lowlim1_X);hold on;
xline(uplim1_X);hold on;
yline(lowlim1_Y);hold on;
yline(uplim1_Y);hold on;
saveas(Figure3, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure3.png');
hold off; 



% % FIGURE 4
% Plot HIEL centers on XY plot
Y1=Y;
Figure4 =figure; 
plot(X,Y1,'.','MarkerSize',10, 'Color', 'm'); hold on;
xlim([0 1024])
ylim([0 1024])
saveas(Figure4, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure4.png');
hold off;



% % FIGURE 5
% Plot HIEL centers with respect to other fluorescent signal of interest
Figure5 =figure;  
plot(claudin5_full_M9V7Z1A1(:,1), claudin5_full_M9V7Z1A1(:,2),'.','MarkerSize',2, 'Color',  'c'); hold on
plot(coord_M9V7Z1A1(:,1), coord_M9V7Z1A1(:,2),'.','MarkerSize',12, 'Color', 'm'); hold on;
xlim([0 1024]);
ylim([0 1024]);
saveas(Figure5, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure5.png');
hold off;



% % FIGURE 6
% Visual check to verify accuracy of thresholding
% Plot centers and secondary fluorescent signal of interest ontop of the
% original fluorescent image 
RGB_merge = imread('M9V7Z1_DAPI_cldn5_PS_hydr.png');
Figure6_copy =figure;
imshow(RGB_merge);
hold on;
plot(claudin5_full_M9V7Z1A1(:,1), claudin5_full_M9V7Z1A1(:,2),'.','MarkerSize',2, 'Color',  'w')
hold on
plot(coord_M9V7Z1A1(:,1), coord_M9V7Z1A1(:,2),'.','MarkerSize',12, 'Color', 'm')
hold on;

colors3=('k');
for k = 1:length(F)
  boundary_J = F{k};
    cidx = mod(k,length(colors))+1;
    plot(boundary_J(:,2), boundary_J(:,1),...
       colors3(cidx),'LineWidth',2);

% %   add text to show object numbers 
  rndRow = ceil(length(boundary_J)/(mod(rand*k,7)+1));
  col = boundary_J(rndRow,2); row = boundary_J(rndRow,1);
  h = text(col+1, row-1, num2str(L(row,col)));
  set(h,'Color',colors2(cidx),'FontSize',8,'FontWeight','bold');
%   randomize text position for better visibility

end
saveas(Figure6_copy, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure16_copy.png');
hold off;
% 
% 
% 
% 
% % 
% % % % % 
% % % % % % % % %  SECTION FOUR
% % % % % 
% % 
% 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % Begin Monte Carlo simulation
ite = 1;
D = cell(100, 1);
Min_d = cell(100,1);
Min_d_hist = cell(100,1);
distance_mc = zeros(mmmm,nnn);

% This loop commnicates with the ExampleMonteCarlo file.m
for jj=1:sim2 % sim2 is defined at the beginning of this code and indicates how many times the simulation will be run
    [distance_mc_sort_M9V7Z1A1,distance_mc_avgMIN_M9V7Z1A1,distance_mc_avgMIN_hist_M9V7Z1A1] = ExampleMonteCarlo(claudin5_full_M9V7Z1A1);
    D{jj} = distance_mc_sort_M9V7Z1A1;
    Min_d{jj} = distance_mc_avgMIN_M9V7Z1A1;
    Min_d_hist{jj} = distance_mc_avgMIN_hist_M9V7Z1A1;
end

% Analyze Monte Carlo simulation
Min_d; 

D;
%top 10000 is defined in the monte carlo simulation
Dtr = cellfun(@(x)transpose(x), D(:,1),'uniformOutput',false); % transpose

Dd_new_M9V7Z1A1 = cellfun(@(x)sort(x), D,'uniformOutput',false); %sort each cell from low to high
    %average each cell element by element to get 1 final matrix. then use this
    %for ddlow etc.
Dd_new_avg_M9V7Z1A1 = mean(cat(3,Dd_new_M9V7Z1A1{:}),3); %average output for MC simulations, sort high to low for each HIEL

% %sort by column (each random circle) low to high (distance to claudin 5)
Ddlow_M9V7Z1A1 = sort(Dd_new_avg_M9V7Z1A1); 


% lowest 100 distances to claudin 5 for each of the simulated IEL holes
Total = 100;
Dd_new_avg_top1000 = Dd_new_avg_M9V7Z1A1(1:Total); 
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%Begin simulated positive control 
%
ite = 1;
D_pc = cell(100, 1);
Min_p_pc = cell(100,1);
Min_d_hist_pc = cell(100,1);
distance_pc = zeros(mmmm,nnn);

% This loop commnicates with the ExamplePositiveControl file.m
for jj=1:sim2 % sim2 is defined at the beginning of this code and indicates how many times the simulation will be run.
    [distance_pc_sort_M9V7Z1A1,distance_pc_avgMIN_M9V7Z1A1,distance_pc_avgMIN_hist_M9V7Z1A1] = ExamplePositiveControl(claudin5_full_M9V7Z1A1,coord_M9V7Z1A1);
    D_pc{jj} = distance_pc_sort_M9V7Z1A1;
    Min_p_pc{jj} = distance_pc_avgMIN_M9V7Z1A1;
    Min_d_hist_pc{jj} = distance_pc_avgMIN_hist_M9V7Z1A1;
end

% Analyze positive control2 simulation
Min_p_pc; 

D_pc;
%top 100 is defined in the positive control simulation
Dtr_pc = cellfun(@(x)transpose(x), D_pc(:,1),'uniformOutput',false); % transpose

Dd_new_pc_M9V7Z1A1 = cellfun(@(x)sort(x), D_pc,'uniformOutput',false); %sort each cell from low to high
    %average each cell element by element to get 1 final matrix. then use this
    %for ddlow etc.
Dd_new_avg_pc_M9V7Z1A1 = mean(cat(3,Dd_new_pc_M9V7Z1A1{:}),3); %average output for MC simulations, sort high to low for each IEL hole

% %sort by column (each circle) low to high (distance to claudin 5)
Ddlow_pc_M9V7Z1A1 = sort(Dd_new_avg_pc_M9V7Z1A1,2); 

% lowest 100 distances to claudin 5 for each of the simulated IEL holes
Total = 100;
Dd_new_avg_top_pce_M9V7Z1A1 = Ddlow_pc_M9V7Z1A1(1:Total);
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%Begin simulated negative control 
%
ite = 1;
D_nc = cell(100, 1);
Min_d_nc = cell(100,1);
Min_d_hist_nc = cell(100,1);
distance_nc = zeros(mmmm,nnn);

% This loop commnicates with the ExampleNegativeControl file.m
for jj=1:sim2 % sim2 is defined at the beginning of this code and indicates how many times the simulation will be run.
    [distance_nc_sort_M9V7Z1A1,distance_nc_avgMIN_M9V7Z1A1,distance_nc_avgMIN_hist_M9V7Z1A1] = ExampleNegativeControl(claudin5_full_M9V7Z1A1,coord_M9V7Z1A1);
    D_nc{jj} = distance_nc_sort_M9V7Z1A1;
    Min_d_nc{jj} = distance_nc_avgMIN_M9V7Z1A1;
    Min_d_hist_nc{jj} = distance_nc_avgMIN_hist_M9V7Z1A1;
end

% Analyze neg control simulation
Min_d_nc; 

D_nc;
%top 100 is defined within the neg control file.
Dtr_nc_M9V7Z1A1 = cellfun(@(x)transpose(x), D_nc(:,1),'uniformOutput',false); % transpose
        
Dd_new_nc_M9V7Z1A1 = cellfun(@(x)sort(x), D_nc,'uniformOutput',false); %sort each cell from low to high
    %average each cell element by element to get 1 final matrix. then use this
    %for ddlow etc.
Dd_new_avg_nc_M9V7Z1A1 = mean(cat(3,Dd_new_nc_M9V7Z1A1{:}),3); %average output for MC simulations, sort high to low for each IEL hole

%sort by column (each circle) low to high (distance to claudin 5)
Ddlow_nc_M9V7Z1A1 = sort(Dd_new_avg_nc_M9V7Z1A1,2); 

% lowest 100 distances to claudin 5 for each of the simulated IEL holes
Total = 100;
Dd_new_avg_top_pce_M9V7Z1A1 = Ddlow_nc_M9V7Z1A1(1:Total);
%
%
%
% 
% % 
% % % % % 
% % % % % % % % %  SECTION FIVE
% % % % % 
% % 
% 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%PUNCTA DETECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%PUNCTA DETECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%PUNCTA DETECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%PUNCTA DETECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define RGB file. Must be an RGB image (can do that in ImageJ).
% This file is the puncta that are localized to HIEL.
RGB = imread('M9V7Z1_PS.png'); 

% convert to grayscale
I_p=rgb2gray(RGB); 
figure; imshow(I_p); % display grayscale

% Blur the image.
P = fspecial('motion',1); 
% Blur step 2.
I2_p = imfilter(I_p,P,'replicate'); 

figure; imshow(I2_p); %display after blurring

% Idenitfy puncta of interest by converting the blurred image to binary.
bw_p =imbinarize(I2_p,'adaptive','Sensitivity',0.1); %adjust to optimize

% Refine the identification of puncta.
SE = strel('disk',3); %adjust to optimize
I3_p = imerode(bw_p,SE);%adjust to optimize
I4_p = bwpropfilt(I3_p,'Area',[2 100]);

Figure2p = figure; hold on; imshow(I4_p); hold off;

% Read in more RGB images to use in figures throughout section
RGB_4 = imread('M9V7Z1_PS_hydr.png');
RGB_merge = imread('M9V7Z1_DAPI_cldn5_PS_hydr.png');


% Threshold puncta fluoresence images into 20 different levels
thresh = multithresh(I_p,20); % adjust levels bassed  on how strong or weak your puncta signal are

imshow(I_p);
seg_I = imquantize(I_p, thresh);
RGB_2 = label2rgb(seg_I);
figure;
imshow(RGB_2)

% Consider only the 20th level of threshold (highest) as the definition of
% puncta.
I4_p =  seg_I~=0 & seg_I~=1 & seg_I~=2 & seg_I~=3 &  seg_I~=4 &...
 seg_I~=5 & seg_I~=6 & seg_I~=7 & seg_I~=8 & seg_I~=9 & seg_I~=10 & seg_I~=11 &...
 seg_I~=12 & seg_I~=13 & seg_I~=14 & seg_I~=15 & seg_I~=16 & seg_I~=17 & seg_I~=18 & seg_I~=19;% & seg_I~=8

[B_p,L_p] = bwboundaries(I4_p); 


% % FIGURE 3P
% Plot HIEL centers with respect to other fluorescent signal of interest

Figure3p_test = figure; 
RGB_3 = label2rgb(I4_p);

imshow(RGB_4); hold on;

colors=('k');
for k = 1:length(B_p)
    boundary = B_p{k};
    cidx = mod(k,length(colors))+1;
    plot(boundary(:,2), boundary(:,1),...
       colors(cidx),'LineWidth',2);

% %   add text to show object numbers 
  rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
  col = boundary(rndRow,2); row = boundary(rndRow,1);
  h = text(col+1, row-1, num2str(L_p(row,col)));
  set(h,'Color',colors(cidx),'FontSize',8);
% %   randomize text position for better visibility

end
hold on;

saveas(Figure3p_test, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure3p_test.png');
hold off;


G_p_check = cell2mat(B_p);
X_p_check = G_p_check(:,2:2:end);
Y_p_check = G_p_check(:,1:2:end);
Figure_verifypuncta =figure; imshow(RGB_4); hold on; plot(X_p_check, Y_p_check, '.', 'Color', 'y', 'MarkerSize',5);
saveas(Figure_verifypuncta, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure_verifypuncta.png');

hold off; 
    
    
% Calculate centerpoints of puncta    

Figure_p_centers =figure; 
hold on;
Q_p = cellfun(@mean,B_p, 'UniformOutput', false); 
% get the average of the pixels along the border to calculate center

Q_p(cellfun(@(Q_p) any(isnan(Q_p)),Q_p)) = [];
H_p = cell2mat(Q_p);
% NOTE both Q and b are cell arrays

% plots the X_p,Y_p coordinates of H (centers of puncta)
X_p = H_p(:,2:2:end);
Y_p = H_p(:,1:2:end);
plot(X_p,Y_p,'.','MarkerSize',5, 'Color', 'm');

saveas(Figure_p_centers, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure_p_centers.png');
hold off;
coord_p_centers_M9V7Z1A1 = horzcat(round(X_p),round(Y_p));


% % UNCOMMENT if want to analyze within certain bounds.
% % Only consider puncta within the specified bounds of analysis
% coord_p_centers_M9V7Z1A1(coord_p_centers_M9V7Z1A1(:, 1) < lowlim1_X & (coord_p_centers_M9V7Z1A1(:, 1) >= Pixel_min), :) = [];
% coord_p_centers_M9V7Z1A1(coord_p_centers_M9V7Z1A1(:, 1) <= Pixel_max & (coord_p_centers_M9V7Z1A1(:, 1) > uplim1_X), :) = [];
% coord_p_centers_M9V7Z1A1(coord_p_centers_M9V7Z1A1(:, 2) < lowlim1_Y & (coord_p_centers_M9V7Z1A1(:, 2) >= Pixel_min), :) = [];
% coord_p_centers_M9V7Z1A1(coord_p_centers_M9V7Z1A1(:, 2) <= Pixel_max & (coord_p_centers_M9V7Z1A1(:, 2) > uplim1_Y), :) = [];


% Uncomment this this code if want to consider only a specific area of image

% % EXAMPLE: if only want to consider ones that are between 200<X_p<700
% % where 200 and 700 are pixel bounds
% lowlim1 = 200;
% uplim1 = 700;
% puncta_final(puncta_final(:,1) < lowlim1, :) = []; 
% puncta_final(puncta_final(:,1) > uplim1, :) = []; 

puncta_final_M9V7Z1A1 = coord_p_centers_M9V7Z1A1;



% % FIGURE PUNCTA CHECK
% Plot detected puncta centers on top of IEL image.

FigurePunctaCheck = figure;
RGb_merge2 = imread('M9V7Z1_hydr.png');
imshow(RGb_merge2); hold on;  

        
 for i =1:length(coord_p_centers_M9V7Z1A1)
   plot(coord_p_centers_M9V7Z1A1(i,1),coord_p_centers_M9V7Z1A1(i,2), '.','MarkerSize',12, 'Color',  'b')
   labels2(i) = cellstr(num2str(([0+i])));
   text(coord_p_centers_M9V7Z1A1(i,1), coord_p_centers_M9V7Z1A1(i,2), labels2(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'w','FontSize',6);             
 end 

xline(lowlim1_X, 'Color','w');hold on;
xline(uplim1_X,'Color','w');hold on;
yline(lowlim1_Y,'Color','w');hold on;
yline(uplim1_Y,'Color','w');hold on;
FigurePunctaOverlapCheck_2.InvertHardcopy = 'off';


 saveas(FigurePunctaCheck, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/FigurePunctaCheck.png');
 hold off;


 
% % FIGURE PUNCTA CHECK 2
% Plot detected puncta centers on top of IEL+puncta image

FigurePunctaCheck2 = figure;
imshow(RGB_4); hold on;  

        
 for i =1:length(coord_p_centers_M9V7Z1A1)
   plot(coord_p_centers_M9V7Z1A1(i,1),coord_p_centers_M9V7Z1A1(i,2), '.','MarkerSize',12, 'Color',  'b')
   labels2(i) = cellstr(num2str(([0+i])));
   text(coord_p_centers_M9V7Z1A1(i,1), coord_p_centers_M9V7Z1A1(i,2), labels2(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'w','FontSize',6);             
 end 

 
xline(lowlim1_X, 'Color','w');hold on;
xline(uplim1_X,'Color','w');hold on;
yline(lowlim1_Y,'Color','w');hold on;
yline(uplim1_Y,'Color','w');hold on;
FigurePunctaOverlapCheck_2.InvertHardcopy = 'off';


 saveas(FigurePunctaCheck2, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/FigurePunctaCheck2.png');
 hold off;


% % UNCOMMENT if want to analyze within certain bounds.

% puncta_final_M9V7Z1A1(puncta_final_M9V7Z1A1(:, 1) < lowlim1_X & (puncta_final_M9V7Z1A1(:, 1) > Pixel_min), :) = [];
% puncta_final_M9V7Z1A1(puncta_final_M9V7Z1A1(:, 1) < Pixel_max & (puncta_final_M9V7Z1A1(:, 1) > uplim1_X), :) = [];
% puncta_final_M9V7Z1A1(puncta_final_M9V7Z1A1(:, 2) < lowlim1_Y & (puncta_final_M9V7Z1A1(:, 2) > Pixel_min), :) = [];
% puncta_final_M9V7Z1A1(puncta_final_M9V7Z1A1(:, 2) < Pixel_max & (puncta_final_M9V7Z1A1(:, 2) > uplim1_Y), :) = [];

%puncta_final(puncta_final(:, 1) < 400 & (puncta_final(:, 1) > 780) & (puncta_final(:, 2) < 200) & (puncta_final(:,2)>210), :) = [];
%Remove when  X_p 358-382 and Y_p 709 - 724 
% puncta_final(puncta_final(:, 1) < 395 & (puncta_final(:, 1) > 350) & (puncta_final(:, 2) < 745) & (puncta_final(:,2)>700), :) = [];




% % % 
% % % 
% % % 
% % % The next several figures are visual checks to confirm that the code
% % % is thresholding and detecting the correct objects.
% % % 
% % % 
% % % 


Figure9p = figure;
imshow(RGB); hold on;
plot(puncta_final_M9V7Z1A1(:,1), puncta_final_M9V7Z1A1(:,2), '.','MarkerSize',2, 'Color',  'c'); hold on;
saveas(Figure9p, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure9p.png');
xline(lowlim1_X, 'Color','w');hold on;
xline(uplim1_X,'Color','w');hold on;
yline(lowlim1_Y,'Color','w');hold on;
yline(uplim1_Y,'Color','w');hold on;
hold off;

RGB_merge = imread('M9V7Z1_PS.png');
Figure9p = figure;
imshow(RGB_merge); hold on;
plot(puncta_final_M9V7Z1A1(:,1), puncta_final_M9V7Z1A1(:,2), '.','MarkerSize',2, 'Color',  'c'); hold on;
saveas(Figure9p, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure9p.png');
xline(lowlim1_X, 'Color','w');hold on;
xline(uplim1_X,'Color','w');hold on;
yline(lowlim1_Y,'Color','w');hold on;
yline(uplim1_Y,'Color','w');hold on;
hold off;

RGB_merge = imread('M9V7Z1_PS.png');
Figure_puncta_centers = figure;
imshow(RGB_merge); hold on;
plot(coord_p_centers_M9V7Z1A1(:,1), coord_p_centers_M9V7Z1A1(:,2), '.','MarkerSize',2, 'Color',  'w'); hold on;
saveas(Figure_puncta_centers, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure_puncta_centers.png');
xline(lowlim1_X, 'Color','w');hold on;
xline(uplim1_X,'Color','w');hold on;
yline(lowlim1_Y,'Color','w');hold on;
yline(uplim1_Y,'Color','w');hold on;
hold off;


% 
% % 
% % % % % 
% % % % % % % % %  PUNCTA CALCULATIONS
% % % % % 
% % 
% 

% Calculate distance of puncta centers to HIEL centers.

for k = 1:length(coord_p_centers_M9V7Z1A1)
   label_p(k) = k;
end

label_p = transpose(label_p);

% Define distance of PS puncta to HIEL center that defines it as
% being in an HIEL.
Threshold_p = 0.75; % in microns

% Calculate the distances of ALL puncta to ALL HIEL
distance_p_M9V7Z1A1 = pdist2(coord_M9V7Z1A1,coord_p_centers_M9V7Z1A1);

% Sort distances
distance_p_M9V7Z1A1 = sort(distance_p_M9V7Z1A1);
distance_p_ALLPUNCTA_M9V7Z1A1=(transpose(distance_p_M9V7Z1A1(1,:))); % distance of all puncta to nearest IEL hole

% Convert distances to microns
distance_p_M9V7Z1A1 = distance_p_M9V7Z1A1/scale_M9V7Z1A1;
distance_p_thresh = distance_p_M9V7Z1A1;

% Remove distances that are above the threshold.
distance_p_thresh(distance_p_thresh > Threshold_p) = NaN;
distance_p_min_M9V7Z1A1 = transpose(distance_p_thresh(1,:));

distance_p_min_final_M9V7Z1A1 = distance_p_min_M9V7Z1A1(~isnan(distance_p_min_M9V7Z1A1));

% Add in puncta labels to the remaining distances in order to
% identify the puncta that overlap with HIEL.
distance_p_min_M9V7Z1A1 = horzcat(label_p, distance_p_min_M9V7Z1A1);
distance_p_min_index_M9V7Z1A1  = distance_p_min_M9V7Z1A1(~isnan(distance_p_min_M9V7Z1A1(:,2))); 

% Report puncta label and its distance to HIEL in a matrix.
distance_p_min_final_M9V7Z1A1 = horzcat(distance_p_min_index_M9V7Z1A1, distance_p_min_final_M9V7Z1A1); 


% This matrix contains the puncta that were identified as
% overlapping with HIEL, with their distances to HIEL centers
% reported (within the thresholded limit of Threshold_p =0.75
% microns).
distance_p_min_final_M9V7Z1A1_copy_orig = distance_p_min_final_M9V7Z1A1; 
% 
% 
% 
% % % % 
% % % % 
% % % % MANUALLY REMOVE PUNCTA (if they were incorrectly identified as overlapping with HIEL)
% % % % using row numbers of distance_p_min_final_M9V7Z1A1 matrix
% % % % and the "TF" system below.
% % % % 
% % % % 
% 
% 
% 
TF1 =distance_p_min_final_M9V7Z1A1(:,1) ==30;
TF2 =distance_p_min_final_M9V7Z1A1(:,1) ==31;
TF3 =distance_p_min_final_M9V7Z1A1(:,1) ==60;
TF4 =distance_p_min_final_M9V7Z1A1(:,1) ==59;
TF5 =distance_p_min_final_M9V7Z1A1(:,1) ==66;
%       TF6 =distance_p_min_final_M9V7Z1A1(:,1) ==NEWNUMBER;
%         
TFall = TF1 | TF2 | TF3 | TF4 | TF5 ; % | TF6 ;
distance_p_min_final_M9V7Z1A1(TFall, :)=[];
%        
%     


% %        
% % FIGURE
% % Visual confirmation that threshold of 0.75um correctly identifies a
% % puncta overlapping with HIEL.
% %    


RGb_merge2 = imread('M9V7Z1_hydr.png'); %%this is the figure that has the puncta and the IEL holes only
FigurePunctaOverlapCheck = figure;
PunctaIndex_M9V7Z1A1 = horzcat(label_p, coord_p_centers_M9V7Z1A1);
imshow(RGb_merge2); hold on;  
  
% These are the labels that correspond to the puncta identifier. Must MANUALLY copy
% and paste these from distance_p_min_final_M9V7Z1A1. 
labels = cellstr(num2str(([[[[[[3;12;13;22;23;24;32;37;38;39;41;43;55;56;61;63;64;65;67;68]]]]]])));

        
 for i = distance_p_min_final_M9V7Z1A1(:,1)
   plot(PunctaIndex_M9V7Z1A1(i,2), PunctaIndex_M9V7Z1A1(i,3), '.','MarkerSize',12, 'Color',  'b')
   text(PunctaIndex_M9V7Z1A1(i,2), PunctaIndex_M9V7Z1A1(i,3), labels, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','Color','m','FontSize',18);             
 end
 
xline(lowlim1_X, 'Color','w');hold on;
xline(uplim1_X,'Color','w');hold on;
yline(lowlim1_Y,'Color','w');hold on;
yline(uplim1_Y,'Color','w');hold on;
FigurePunctaOverlapCheck_2.InvertHardcopy = 'off';

 saveas(FigurePunctaOverlapCheck, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/FigurePunctaOverlapCheck.png');
 hold off;

% %        
% % FIGURE
% % Visual confirmation that threshold of 0.75um correctly identifies a
% % puncta overlapping with HIEL.
% %    
       
FigurePunctaOverlapCheck_2 = figure;
PunctaIndex_M9V7Z1A1 = horzcat(label_p, coord_p_centers_M9V7Z1A1);
imshow(RGB_merge); hold on;  
  

        
 for i = distance_p_min_final_M9V7Z1A1(:,1)
   plot(PunctaIndex_M9V7Z1A1(i,2), PunctaIndex_M9V7Z1A1(i,3), '.','MarkerSize',12, 'Color',  'b')
   text(PunctaIndex_M9V7Z1A1(i,2), PunctaIndex_M9V7Z1A1(i,3), labels, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'm','FontSize',18);             
 end
 
xline(lowlim1_X, 'Color','w');hold on;
xline(uplim1_X,'Color','w');hold on;
yline(lowlim1_Y,'Color','w');hold on;
yline(uplim1_Y,'Color','w');hold on;
FigurePunctaOverlapCheck_2.InvertHardcopy = 'off';

saveas(FigurePunctaOverlapCheck_2, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/FigurePunctaOverlapCheck_2.png');
hold off;

        
% Calculate the distance of puncta to claudin5       
distance_p_cldn5_M9V7Z1A1 = pdist2(claudin5_full_M9V7Z1A1,coord_p_centers_M9V7Z1A1);
distance_p_cldn5_M9V7Z1A1 = sort(distance_p_cldn5_M9V7Z1A1);
distance_p_cldn5_M9V7Z1A1 = distance_p_cldn5_M9V7Z1A1/scale_M9V7Z1A1;
distance_p_cldn5_M9V7Z1A1 = transpose(distance_p_cldn5_M9V7Z1A1(1,:));

distance_p_cldn5_M9V7Z1A1_Label = horzcat(label_p, distance_p_cldn5_M9V7Z1A1);
        
  
% These are the labels that correspond to the puncta identifier. Must MANUALLY copy
% and paste these from distance_p_min_final_M9V7Z1A1.       
distance_p_cldn5_overlap_M9V7Z1A1 = distance_p_cldn5_M9V7Z1A1([[[[[[[3 12 13 22 23 24 32 37 38 39 41 43 55 56 61 63 64 65 67 68]]]]]]]); 
distance_p_min_final_overlap_M9V7Z1A1 = horzcat(distance_p_min_final_M9V7Z1A1,distance_p_cldn5_overlap_M9V7Z1A1); %  #puncta, dist to nearest HIEL, dist to nearest cldn5
       

% Final matirx that includes puncta number, distance to closest HIEL, and minimum distance to claudin5.        
distance_ALLPUNCTA_final_M9V7Z1A1= horzcat(label_p,  distance_p_ALLPUNCTA_M9V7Z1A1, distance_p_cldn5_M9V7Z1A1); % #puncta, dist to nearest HIEL, dist to nearest cldn5
               
        
% Calculate proportion of HIEL that have a puncta
% number of IEL holes
Total_HIEL_M9V7Z1A1 = length(coord_M9V7Z1A1);
PS_pos_M9V7Z1A1 = length(distance_p_min_index_M9V7Z1A1); %length of 
PercentPosHIEL_M9V7Z1A1 = PS_pos_M9V7Z1A1/Total_HIEL_M9V7Z1A1;
        

% Calculate proportion of puncta overlying HIEL.
Total_Puncta_M9V7Z1A1 = length(coord_p_centers_M9V7Z1A1);
PercentPunctaWithHIELOverlap_M9V7Z1A1 = PS_pos_M9V7Z1A1/Total_Puncta_M9V7Z1A1;
% 
% 
% 
% 
% % 
% % % % % 
% % % % % % % % %  SECTION SIX
% % % % % 
% % 
% 
%         
%
%
% % begin Monte Carlo simulation for puncta to claudin5
ite = 1;
D = cell(100, 1);
Min_d_nc = cell(100,1);
Min_d_hist_nc = cell(100,1);
distance_mc_p = zeros(mmmm,nnn);
N_mc_p = length(coord_M9V7Z1A1); 


% This loop communicates with ExampleMonteCarloPuncta.m
for jj=1:sim1 % sim1 is defined at the beginning of this code and indicates how many times the simulation will be run.
    [distance_mc_p_sort_M9V7Z1A1,distance_mc_p_avgMIN_M9V7Z1A1,distance_mc_p_avgMIN_M9V7Z1A1_hist] = ExampleMonteCarloPuncta(claudin5_full_M9V7Z1A1,N_mc_p);
    D{jj} = distance_mc_p_sort_M9V7Z1A1;
    Min_d_nc{jj} = distance_mc_p_avgMIN_M9V7Z1A1;
    Min_d_hist_nc{jj} = distance_mc_p_avgMIN_M9V7Z1A1_hist;
end


% Analyze Monte Carlo simulation
Min_d_nc; 

D;
%top 100 is defined in the monte carlo simulation
Dtr = cellfun(@(x)transpose(x), D(:,1),'uniformOutput',false); % transpose

        
Dd_new = cellfun(@(x)sort(x), D,'uniformOutput',false); %sort each cell from low to high
    %average each cell element by element to get 1 final matrix. then use this
    %for Ddlow_p_M9V7Z1A1 etc.
Dd_new_p_avg_M9V7Z1A1 = mean(cat(3,Dd_new{:}),3); %average output for mc_p simulations, sort high to low for each IEL hole

Ddlow_p_M9V7Z1A1 = sort(Dd_new_p_avg_M9V7Z1A1); %sort by column (each random circle) low to high (distance to claudin 5)
Ddlow_p_M9V7Z1A12 = sort(Dd_new_p_avg_M9V7Z1A1);
Ddlow_p_M9V7Z1A13 = sort(Dd_new_p_avg_M9V7Z1A1);
Ddlow_p_M9V7Z1A14 = sort(Dd_new_p_avg_M9V7Z1A1);
Ddlow_p_M9V7Z1A15 = sort(Dd_new_p_avg_M9V7Z1A1);

% % lowest 100 distances to claudin 5 for each of the simulated IEL holes
Total = 100;
Dd_new_p_avg_M9V7Z1A1_top1000 = Dd_new_p_avg_M9V7Z1A1(1:Total); 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%begin simulated positive control 

ite = 100;
D_pc = cell(100, 1);
Min_p_pc = cell(100,1);
Min_d_hist_pc = cell(100,1);
distance_pc = zeros(mmmm,nnn);

% This code communicates with ExamplePositiveControlPuncta.m
for jj=1:sim1 % sim1 is defined at the beginning of this code and indicates how many times the simulation will be run.
    [distance_pc_p_sort_M9V7Z1A1,distance_pc_p_avgMIN_M9V7Z1A1,distance_pc_p_avgMIN_M9V7Z1A1_hist] = ExamplePositiveControlPuncta(claudin5_full_M9V7Z1A1,coord_p_centers_M9V7Z1A1,PunctaIndex_M9V7Z1A1);
    D_pc{jj} = distance_pc_p_sort_M9V7Z1A1;
    Min_p_pc{jj} = distance_pc_p_avgMIN_M9V7Z1A1;
    Min_d_hist_pc{jj} = distance_pc_p_avgMIN_M9V7Z1A1_hist;
end

% Analyze positive control simulation
Min_p_pc; 

D_pc;
%top 100 is defined in the positive control simulation
Dtr_pc = cellfun(@(x)transpose(x), D_pc(:,1),'uniformOutput',false); % transpose

Dd_new_p_pc = cellfun(@(x)sort(x), D_pc,'uniformOutput',false); %sort each cell from low to high

Dd_new_p_avg_M9V7Z1A1_pc = mean(cat(3,Dd_new_p_pc{:}),3); %average output for mc_p simulations, sort high to low for each IEL hole

Ddlow_p_M9V7Z1A1_pc = sort(Dd_new_p_avg_M9V7Z1A1_pc,2); %sort by column (each random circle) low to high (distance to claudin 5)
Ddlow_p_M9V7Z1A12_pc = sort(Dd_new_p_avg_M9V7Z1A1_pc,2);
Ddlow_p_M9V7Z1A13_pc = sort(Dd_new_p_avg_M9V7Z1A1_pc,2);
Ddlow_p_M9V7Z1A14_pc = sort(Dd_new_p_avg_M9V7Z1A1_pc,2);
Ddlow_p_M9V7Z1A15_pc = sort(Dd_new_p_avg_M9V7Z1A1_pc,2);

% % lowest 100 distances to claudin 5 for each of the simulated IEL holes
Total = 100;
Dd_new_p_avg_M9V7Z1A1_top_pce = Ddlow_p_M9V7Z1A1_pc(1:Total);
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%begin simulated negative control 

ite = 1;
D_nc = cell(100, 1);
Min_d_nc = cell(100,1);
Min_d_hist_nc = cell(100,1);
distance_nc = zeros(mmmm,nnn);

% This code communicates with ExampleNegativeControlPuncta.m
for jj=1:sim1 %sim1 is defined at the beginning of this code and indicates how many times the simulation will be run.
    [distance_nc_p_sort_M9V7Z1A1,distance_nc_p_avgMIN_M9V7Z1A1,distance_nc_p_avgMIN_M9V7Z1A1_hist] = ExampleNegativeControlPuncta(claudin5_full_M9V7Z1A1,coord_p_centers_M9V7Z1A1, PunctaIndex_M9V7Z1A1);
    D_nc{jj} = distance_nc_p_sort_M9V7Z1A1;
    Min_d_nc{jj} = distance_nc_p_avgMIN_M9V7Z1A1;
    Min_d_hist_nc{jj} = distance_nc_p_avgMIN_M9V7Z1A1_hist;
end

% Analyze neg control simulation
Min_d_nc; 

D_nc;
%top 100 is defined within the neg control file.
Dtr_nc = cellfun(@(x)transpose(x), D_nc(:,1),'uniformOutput',false); % transpose
        
Dd_new_p_nc = cellfun(@(x)sort(x), D_nc,'uniformOutput',false); %sort each cell from low to high
Dd_new_p_avg_M9V7Z1A1_nc = mean(cat(3,Dd_new_p_nc{:}),3); %average output for mc_p simulations, sort high to low for each IEL hole
Ddlow_p_M9V7Z1A1_nc = sort(Dd_new_p_avg_M9V7Z1A1_nc,2); %sort by column (each random circle) low to high (distance to claudin 5)
Ddlow_p_M9V7Z1A12_nc = sort(Dd_new_p_avg_M9V7Z1A1_nc,2);
Ddlow_p_M9V7Z1A13_nc = sort(Dd_new_p_avg_M9V7Z1A1_nc,2);
Ddlow_p_M9V7Z1A14_nc = sort(Dd_new_p_avg_M9V7Z1A1_nc,2);
Ddlow_p_M9V7Z1A15_nc = sort(Dd_new_p_avg_M9V7Z1A1_nc,2);


% % lowest 100 distances to claudin 5 for each of the simulated IEL holes
Total = 100;
Dd_new_p_avg_M9V7Z1A1_top_pce = Ddlow_p_M9V7Z1A1_nc(1:Total);
%
%
%
% 
% % 
% % % % % 
% % % % % % % % %  SECTION SEVEN
% % % % % 
% % 
% 
%
%
%
% Save all the worthwhile variables into a file. 
save('Puncta_Cldn5_Partial_M9V7Z1A1.mat', 'claudin5_full_M9V7Z1A1', 'coord_M9V7Z1A1', 'B_M9V7Z1A1','distance_sort_M9V7Z1A1','distance_avgMIN_M9V7Z1A1', 'distance_avgMIN_hist_M9V7Z1A1',...
    'distance_nc_sort_M9V7Z1A1','distance_nc_avgMIN_M9V7Z1A1', 'distance_nc_avgMIN_hist_M9V7Z1A1','Dd_new_avg_M9V7Z1A1','Dd_new_avg_pc_M9V7Z1A1','Dd_new_avg_nc_M9V7Z1A1', 'Ddlow_M9V7Z1A1','Ddlow_nc_M9V7Z1A1',...
    'Ddlow_pc_M9V7Z1A1','distance_mc_avgMIN_M9V7Z1A1', 'distance_mc_avgMIN_hist_M9V7Z1A1','distance_pc_sort_M9V7Z1A1',...
    'distance_pc_avgMIN_M9V7Z1A1', 'distance_pc_avgMIN_hist_M9V7Z1A1','Param_M9V7Z1A1', 'scale_M9V7Z1A1', 'PartialImage_Height_M9V7Z1A1', 'PartialImage_Width_M9V7Z1A1','ImageArea_M9V7Z1A1',...
    'IELHoles_M9V7Z1A1', 'HolesPerArea_M9V7Z1A1', 'coord_p_centers_M9V7Z1A1', 'puncta_final_M9V7Z1A1', 'distance_p_M9V7Z1A1', 'distance_p_thresh',  'distance_p_min_final_M9V7Z1A1', 'distance_p_min_M9V7Z1A1', ...
     'distance_p_min_index_M9V7Z1A1','Total_HIEL_M9V7Z1A1','PS_pos_M9V7Z1A1', 'PercentPosHIEL_M9V7Z1A1', 'Total_Puncta_M9V7Z1A1', 'PercentPunctaWithHIELOverlap_M9V7Z1A1',  'distance_nc_p_sort_M9V7Z1A1','distance_nc_p_avgMIN_M9V7Z1A1', ...
     'distance_nc_p_avgMIN_M9V7Z1A1_hist','Dd_new_p_avg_M9V7Z1A1','Dd_new_p_avg_M9V7Z1A1_pc','Dd_new_p_avg_M9V7Z1A1_nc', 'Ddlow_p_M9V7Z1A1','Ddlow_p_M9V7Z1A1_nc',...
    'Ddlow_p_M9V7Z1A1_pc','distance_mc_p_avgMIN_M9V7Z1A1', 'distance_mc_p_avgMIN_M9V7Z1A1_hist',...
    'distance_pc_p_sort_M9V7Z1A1','distance_pc_p_avgMIN_M9V7Z1A1', 'distance_pc_p_avgMIN_M9V7Z1A1_hist',  'distance_ALLPUNCTA_final_M9V7Z1A1', 'distance_p_min_final_overlap_M9V7Z1A1','distance_p_cldn5_overlap_M9V7Z1A1');

 
