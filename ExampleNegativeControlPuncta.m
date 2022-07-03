function [distance_nc_sort_M9V7Z1A1, distance_nc_avgMIN_M9V7Z1A1, distance_nc_avgMIN_hist_M9V7Z1A1] = ExampleNegativeControlPuncta(claudin5_full_M9V7Z1A1,coord_M9V7Z1A1,PunctaIndex_M9V7Z1A1)
N_nc = 64;  % #HIEL identified in ExampleCode
scale =8.6719; %scale for how many pixels equal 1 micron
AvgRadius = 9.66;%average HIEL radius in pixels identfied in ExampleCode
totalPoints = size(claudin5_full_M9V7Z1A1,1);

lowlimX = 0;
upperlimX = 1024;
newX_nc_2 = (lowlimX-upperlimX).*rand(10000,1) + upperlimX;

lowlimY = 0;
upperlimY = 1024;
newY_nc_2 = (lowlimY-upperlimY).*rand(10000,1) + upperlimY;

NC_xycoord_2 = [round(newX_nc_2) round(newY_nc_2)];
x_nc = NC_xycoord_2(:,1);
y_nc = NC_xycoord_2(:,2); 


minAllowableDistance = 1;
numberOfPoints = N_nc+1; %make this equal to 1 more than the number of IEL holes
nMax = 10000;

% Initialize first point.
keeperX = x_nc(1);
keeperY = y_nc(1);
    
% Try dropping down more points.
counter = 2;
k =1 ; %

while counter < numberOfPoints & k < totalPoints
    k=k+1; %%
  % Get a trial point.
  thisX = x_nc(k);
  thisY = y_nc(k);
  thisXY = horzcat(thisX,thisY);
   
  % See how far is is away from existing keeper points.
  distances = sqrt((thisX-keeperX).^2 + (thisY - keeperY).^2);
  distances2 = pdist2(thisXY, claudin5_full_M9V7Z1A1);
  minDistance2 = min(distances2);
  minDistance = min(distances);

    if  length(thisX)< numberOfPoints && minDistance2 >scale*1.5 && minDistance2 <scale*4.5 && minDistance > AvgRadius &&...
            thisX <upperlimX && thisX > lowlimX && thisY <upperlimY && thisY >lowlimY  %% CHANGE FOR PARTIAL
    keeperX(counter) = thisX;
    keeperY(counter) = thisY;
    counter = counter + 1;
    
    else length(keeperX) < numberOfPoints
        k= k+1; 
    end
end


    keeperX = transpose(keeperX);
    keeperY = transpose(keeperY);
    xycoord_NC = horzcat(keeperX,keeperY);
    
    
    N_puncta=20; %MANUALLY ADJUST THIS. Number of HIEL that are filled with punctate in the image 
    xycoord_NC = datasample(xycoord_NC, N_puncta);

    
    Missing_values = numberOfPoints - length(xycoord_NC);


[mmm_nc nnn_nc] = size(xycoord_NC);
nnn_nc = mmm_nc;

[mmmm nnnn] = size(claudin5_full_M9V7Z1A1);
distance_mc = zeros(nnn_nc,mmm_nc);

for i = 1 : nnn_nc % shorter, "coord_M9V7Z1A1"
    for j = 1 : mmmm %longer, "claudin5"
        scale = 8.6719;
        distance_nc_M9V7Z1A1(i, j) = (sqrt((xycoord_NC(i, 1) - claudin5_full_M9V7Z1A1(j,1)) ^ 2 + ...
            (xycoord_NC(i, 2) - claudin5_full_M9V7Z1A1(j,2)) ^ 2))/scale;
        %     slope(i,j) = (coord_M9V7Z1A1(i, 2) - coord_M9V7Z1A1(j, 2))/(coord_M9V7Z1A1(i, 1) - coord_M9V7Z1A1(j, 1));
    end
end

% Calculate the average min distance of simulated puncta to claudin 5 staining.
distance_nc_M9V7Z1A1 = transpose(distance_nc_M9V7Z1A1);
distance_nc_sort_old = sort(distance_nc_M9V7Z1A1, 1); % sort distances within each column
distance_nc_sort_old = sortrows(distance_nc_sort_old.',1).'; %sort by first row
[~,inx]=sort(distance_nc_sort_old(1,:));  %sort by first row
distance_nc_sort_old = distance_nc_sort_old(:,inx);  %sort by first row

% This is the matrix to be used for downstream analysis!!!
% This is the minimum distance of simulated puncta to claudin 5 
distance_nc_avgMIN_hist_M9V7Z1A1 = (distance_nc_sort_old(1,:));

% This is the average value of the matrix above.
distance_nc_avgMIN_M9V7Z1A1= mean(distance_nc_sort_old(1,:));

% Keep top 100 closest distances of HIEL to claudin5 in case needed for
% future analyses.
Total = 100;
distance_nc_sort_M9V7Z1A1 = distance_nc_sort_old(1:Total, 1:N_puncta); % take top 1000 values from this array 



    howmanyzeros_M9V7Z1A1 = nnz(~distance_nc_sort_M9V7Z1A1);

%%Plot the result of one simulation! 
%COMMENT OFF all of this plot when you are ready to run >1 simulation, or
%you will get a plot for every simulation and that takes up a lot of
%memory.
 
%FIGURE 5
Figure_neg_ctrl_partial_cldn5 =figure;  
plot(claudin5_full_M9V7Z1A1(:,1), claudin5_full_M9V7Z1A1(:,2),'.','MarkerSize',2, 'Color',  'c'); hold on
hold on;
 for i = N_puncta
   plot(PunctaIndex_M9V7Z1A1(i,2), PunctaIndex_M9V7Z1A1(i,3), '.','MarkerSize',12, 'Color',  'b')
 end
 hold on;

plot(xycoord_NC(:,1), xycoord_NC(:,2),'.','MarkerSize',12, 'Color', 'k'); hold on;
xlim([0 1024]);
ylim([0 1024]);
legend('Claudin5','1st order artery', 'Negative Control','FontSize',12);
saveas(Figure_neg_ctrl_partial_cldn5, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure_neg_ctrl_partial_cldn5.png');
hold off;

end
