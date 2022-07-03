%distance_mc - nearest cldn5 coordinate to any give hole
function [distance_mc_p_sort_M9V7Z1A1, distance_mc_p_avgMIN_M9V7Z1A1, distance_mc_p_avgMIN_hist_M9V7Z1A1] = ExampleMonteCarloPuncta(claudin5_full_M9V7Z1A1, N_mc_p)

% pixel bounds of analyzed area
m_mc_p = size(N_mc_p,1);
n_mc_p = size (N_mc_p,2);
r_mc_p = zeros(m_mc_p,n_mc_p);
xy_mc_p = zeros(m_mc_p,n_mc_p);

xlow_p = 0;
xhigh_p = 1024;
ylow_p = 0;
yhigh_p =1024;

xhigh_p = xhigh_p - xlow_p;
yhigh_p = yhigh_p - ylow_p;
% % preallocate storage for r for speed
% for k = 1:length(hi_mc)
    for  i=1:N_mc_p
           
            % choose a random radius for the simulated HIEL
            r_mc_p(i) = round(normrnd(9.32,1.99)); 
            % First variable is mean of HIEL size detected in ExampleCode
            % Second is standard deviation of HIEL size detected in ExampleCod
           
            % generate a new HIEL simulated center that is within the
            % bounds of the analyzed area
            newpt_mc_p = @()r_mc_p(i)+ round([((xhigh_p-2*r_mc_p(i))* rand((1))+xlow_p) ((yhigh_p-2*r_mc_p(i))* rand((1))+ylow_p)]);
            
            xy_mc_p = newpt_mc_p();  % matrix to store XY coordinates
            fails = 0;     % to avoid looping forever

            while size(xy_mc_p,1) < N_mc_p
        
           % generate new point and test distance
            pt = newpt_mc_p();
            % test distance to all other points. the goal is to not
                % have two overlapping puncta
                if all(pdist2(xy_mc_p, pt) > 2*r_mc_p(i)) %
                    xy_mc_p = [xy_mc_p; pt];  % add it
                    fails = 0;      % reset failure counter
                else
        % increase failure counter,
                    fails = fails + 1;
        % give up if exceeded some threshold
                        if fails > 5000
                            error('this is taking too long...');
                        end  
                 end
             end 

    end

hold off

%%Plot the result of one simulation! 
%COMMENT OFF all of this plot when you are ready to run >1 simulation, or
%you will get a plot for every simulation and that takes up a lot of
%memory.

N_mc_p_puncta = 20; %MANUALLY ADJUST THIS. Number of HIEL that are filled with punctate in the image 
Simulated_Puncta = datasample(xy_mc_p, N_mc_p_puncta); 
    
% % 
Figure_mc_p_example = figure; hold on;
plot(claudin5_full_M9V7Z1A1(:,1), claudin5_full_M9V7Z1A1(:,2),'.','MarkerSize',8, 'Color',  'c'); hold on;
centers_mc_p = zeros (m_mc_p,n_mc_p);
    for j=1:size(Simulated_Puncta,1)

          centers_mc_p = [xy_mc_p(j,1), xy_mc_p(j,2)];
          viscircles(centers_mc_p, r_mc_p(j), 'color', 'k');
    end
    hold on;
    xlim([0 1024]);
    ylim([0 1024]);
    saveas(Figure_mc_p_example, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure_mc_example.png');
    hold off;

    
% % % % Calculate distances between center points

% Define sizes of matrices.
[mmm_mc_p nnn_mc_p] = size(Simulated_Puncta);

% distance_mc - nearest cldn5 coordinate to any give puncta
distance_mc_p_M9V7Z1A1 = zeros(nnn_mc_p,mmm_mc_p);


for i = 1 : N_mc_p_puncta % shorter, "coord_M9V7Z1A1"
    for j = 1 : length(claudin5_full_M9V7Z1A1) %longer, "claudin5"
        scale = 8.6719;
        distance_mc_p_M9V7Z1A1(i, j) = (sqrt((Simulated_Puncta(i, 1) - claudin5_full_M9V7Z1A1(j, 1)) ^ 2 + ...
            (Simulated_Puncta(i, 2) - claudin5_full_M9V7Z1A1(j, 2)) ^ 2))/scale;
        %     slope(i,j) = (coord_M9V7Z1A1(i, 2) - coord_M9V7Z1A1(j, 2))/(coord_M9V7Z1A1(i, 1) - coord_M9V7Z1A1(j, 1));
    end
end

% Calculate the average min distance of HIEL to claudin 5 staining.
distance_mc_p_M9V7Z1A1 = transpose(distance_mc_p_M9V7Z1A1);
distance_mc_p_sort_old = sort(distance_mc_p_M9V7Z1A1, 1); % sort distances within each column
distance_mc_p_sort_old = sortrows(distance_mc_p_sort_old.',1).'; %sort by first row
[~,inx]=sort(distance_mc_p_sort_old(1,:));  %sort by first row
distance_mc_p_sort_old = distance_mc_p_sort_old(:,inx);  %sort by first row

% This is the matrix to be used for downstream analysis!!!
% This is the minimum distance of simulated HIEL to claudin 5 
distance_mc_p_avgMIN_hist_M9V7Z1A1 = (distance_mc_p_sort_old(1,:));

% This is the average value of the matrix above.
distance_mc_p_avgMIN_M9V7Z1A1= mean(distance_mc_p_sort_old(1,:));


% take top 100 distances measured from this array. the first row indicates 
% the minimum distance measured for each HIEL to claudin5.
% Save 100 just in case it can be used in downstream analysis.
Total = 100;
distance_mc_p_sort_M9V7Z1A1 = distance_mc_p_sort_old(1:Total, 1:N_mc_p_puncta); 

clear distance_mc_p_sort_old; % delete original array
clear distance_mc_p_M9V7Z1A1;


end