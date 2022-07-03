function [distance_pc_sort_M9V7Z1A1, distance_pc_avgMIN_M9V7Z1A1, distance_pc_avgMIN_hist_M9V7Z1A1] = ExamplePositiveControl(claudin5_full_M9V7Z1A1, coord_M9V7Z1A1)
N_pc = 64; % #HIEL identified in ExampleCode
scale =8.6719; %scale for how many pixels equal 1 micron
AvgRadius = 9.66; %average HIEL radius in pixels identfied in ExampleCode
totalPoints = size(claudin5_full_M9V7Z1A1,1);

lowlimX = 0;
upperlimX = 1024;
newX_pc_2 = (lowlimX-upperlimX).*rand(10000,1) + upperlimX;

lowlimY = 0;
upperlimY = 1024;
newY_pc_2 = (lowlimY-upperlimY).*rand(10000,1) + upperlimY;

PC_xycoord_2 = [round(newX_pc_2) round(newY_pc_2)];


x_pc = PC_xycoord_2(:,1);
y_pc = PC_xycoord_2(:,2); 


minAllowableDistance = 1;
numberOfPoints = N_pc+1; % make this equal to 1 more than the number of IEL holes
nMax = 10000;

% Initialize first point.
keeperX = x_pc(1);
keeperY = y_pc(1);

% Try dropping down more points.
counter = 2;
k =1 ; %

while counter < numberOfPoints & k < totalPoints
    k=k+1; %%
  % Get a trial point.
  thisX = x_pc(k);
  thisY = y_pc(k);
  thisXY = horzcat(thisX,thisY);

  % See how far is is away from existing keeper points.
  distances = sqrt((thisX-keeperX).^2 + (thisY - keeperY).^2);
  distances2 = pdist2(thisXY, claudin5_full_M9V7Z1A1);
  minDistance2 = min(distances2);
  minDistance = min(distances);

%   Create the new point if it is within 1.5*scale of fluorescent signal
%   (claudin5 in this example), and if it is within the
%   specified analyzed area.
    if  length(thisX)< numberOfPoints && minDistance2 >=0 && minDistance2 <scale*1.5 && minDistance > AvgRadius &&...
            thisX <upperlimX && thisX > lowlimX && thisY <upperlimY && thisY >lowlimY
        % must alter the above line if only analyzing a partial image
    keeperX(counter) = thisX;
    keeperY(counter) = thisY;
    counter = counter + 1;
    
    else length(keeperX) < numberOfPoints
        k= k+1;
   
    end
end

keeperX = transpose(keeperX);
keeperY = transpose(keeperY);
xycoord_PC = horzcat(keeperX,keeperY);

Missing_values = numberOfPoints - length(xycoord_PC);


[mmm_pc nnn_pc] = size(xycoord_PC);
nnn_pc = mmm_pc;

[mmmm nnnn] = size(claudin5_full_M9V7Z1A1);
distance_pc_M9V7Z1A1 = zeros(nnn_pc,mmm_pc);

for i = 1 : nnn_pc % shorter, "coord_M9V7Z1A1"
    for j = 1 : mmmm %longer, "claudin5"
        scale = 8.6719;
        distance_pc_M9V7Z1A1(i, j) = (sqrt((xycoord_PC(i, 1) - claudin5_full_M9V7Z1A1(j,1)) ^ 2 + ...
            (xycoord_PC(i, 2) - claudin5_full_M9V7Z1A1(j,2)) ^ 2))/scale;
    end
end

% Calculate the average min distance of simulated HIEL to claudin 5 staining.
distance_pc_M9V7Z1A1 = transpose(distance_pc_M9V7Z1A1);
distance_pc_sort_old = sort(distance_pc_M9V7Z1A1, 1); % sort distances within each column
distance_pc_sort_old = sortrows(distance_pc_sort_old.',1).'; %sort by first row
[~,inx]=sort(distance_pc_sort_old(1,:));  %sort by first row
distance_pc_sort_old = distance_pc_sort_old(:,inx);  %sort by first row


% This is the matrix to be used for downstream analysis!!!
% This is the minimum distance of simulated HIEL to claudin 5 
distance_pc_avgMIN_hist_M9V7Z1A1 = (distance_pc_sort_old(1,:));

 % This is the average value of the matrix above.
distance_pc_avgMIN_M9V7Z1A1= mean(distance_pc_sort_old(1,:));


% Keep top 100 closest distances of HIEL to claudin5 in case needed for
% future analyses.
Total = 100;
distance_pc_sort_M9V7Z1A1 = distance_pc_sort_old(1:Total, 1:N_pc); % take top 100 values from this array 

howmanyzeros = nnz(~distance_pc_sort_M9V7Z1A1);

% % COMMENT OUT THIS SECTION WHEN INCREASING NUMBER OF SIMULATIONS >1
% % % %FIGURE 5
% % confirm that the positive control is working with this visual check
Figure_pos_ctrl_partial_cldn5 =figure;  
plot(claudin5_full_M9V7Z1A1(:,1), claudin5_full_M9V7Z1A1(:,2),'.','MarkerSize',2, 'Color',  'c'); hold on
plot(coord_M9V7Z1A1(:,1), coord_M9V7Z1A1(:,2),'.','MarkerSize',12, 'Color', 'm'); hold on;
plot(xycoord_PC(:,1), xycoord_PC(:,2),'.','MarkerSize',12, 'Color', 'b'); hold on;
xlim([0 1024]);
ylim([0 1024]);
legend('Claudin5','1st order artery', 'Positive Control','FontSize',12);
saveas(Figure_pos_ctrl_partial_cldn5, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure_pos_ctrl_partial_cldn5.png');
hold off;
% 

end
