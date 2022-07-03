function [distance_mc_sort_M9V7Z1A1, distance_mc_avgMIN_M9V7Z1A1, distance_mc_avgMIN_hist_M9V7Z1A1] = ExampleMonteCarlo(claudin5_full_M9V7Z1A1)

% Manually adjust. This is the number of simulated HIEL centers (should match # detected in ExampleCode).
N_mc = 64; 

m_mc = size(N_mc,1);
n_mc = size (N_mc,2);
r_mc = zeros(m_mc,n_mc);
xy_mc = zeros(m_mc,n_mc);

% pixel bounds of analyzed area
xlow = 0;
xhigh = 1024;
ylow = 0;
yhigh =1024;

    xhigh = xhigh - xlow;
    yhigh = yhigh - ylow;

    for  i=1:N_mc

            % choose a random radius for the simulated HIEL
            r_mc(i) = round(normrnd(9.66,1.99)); 
            % First variable is mean of HIEL size detected in ExampleCode
            % Second is standard deviation of HIEL size detected in ExampleCode
           
            % generate a new HIEL simulated center that is within the
            % bounds of the analyzed area
            newpt_mc = @()r_mc(i)+ round([((xhigh-2*r_mc(i))* rand((1))+xlow) ((yhigh-2*r_mc(i))* rand((1))+ylow)]);
            
            xy_mc = newpt_mc();  % matrix to store XY coordinates
            fails = 0;     % to avoid looping forever
           
            % generate new HIEl center
            while size(xy_mc,1) < N_mc
            pt = newpt_mc();
                % test distance to all other points. the goal is to not
                % have two overlapping HIEL.
                if all(pdist2(xy_mc, pt) > 2*r_mc(i)) 
                    xy_mc = [xy_mc; pt];  % add it if it passes the test
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



% % % % 
Figure_mc_example = figure; hold on;
plot(claudin5_full_M9V7Z1A1(:,1), claudin5_full_M9V7Z1A1(:,2),'.','MarkerSize',8, 'Color',  'c'); hold on;

%%preallocate space
centers_mc = zeros (m_mc,n_mc);
    for j=1:size(xy_mc,1)
          centers_mc = [xy_mc(j,1), xy_mc(j,2)];
          viscircles(centers_mc, r_mc(j), 'color', 'k');
%         end
    end
    hold on;
    xlim([0 1024]);
    ylim([0 1024]);
    saveas(Figure_mc_example, '/Users/claireruddiman/Dropbox/_MEJ Paper/_Figures/Github/M9V7Z1A1_Analysis/Figure_mc_example.png');
    hold off;
% %end   


% % % % Calculate distances between center points

% Define sizes of matrices.
[mmm_mc nnn_mc] = size(xy_mc);
nnn_mc = mmm_mc;
[mmmm nnnn] = size(claudin5_full_M9V7Z1A1);

% distance_mc - nearest cldn5 coordinate to any give HIEL
distance_mc_M9V7Z1A1 = zeros(nnn_mc,mmm_mc);


for i = 1 : nnn_mc % shorter, "coord_M9V7Z1A1"
    for j = 1 : mmmm %longer, "claudin5"
        scale = 8.6719; % this is the number of pixels per micron
        distance_mc_M9V7Z1A1(i, j) = (sqrt((xy_mc(i, 1) - claudin5_full_M9V7Z1A1(j, 1)) ^ 2 + ...
            (xy_mc(i, 2) - claudin5_full_M9V7Z1A1(j, 2)) ^ 2))/scale;
        %     slope(i,j) = (coord_M9V7Z1A1(i, 2) - coord_M9V7Z1A1(j, 2))/(coord_M9V7Z1A1(i, 1) - coord_M9V7Z1A1(j, 1));
    end
end

% Calculate the average min distance of HIEL to claudin 5 staining.
    distance_mc_M9V7Z1A1 = transpose(distance_mc_M9V7Z1A1);
    distance_mc_sort_old = sort(distance_mc_M9V7Z1A1, 1); % sort distances within each column
    distance_mc_sort_old = sortrows(distance_mc_sort_old.',1).'; %sort by first row
    [~,inx]=sort(distance_mc_sort_old(1,:));  %sort by first row
    distance_mc_sort_old = distance_mc_sort_old(:,inx);  %sort by first row

    % This is the matrix to be used for downstream analysis!!!
    % This is the minimum distance of simulated HIEL to claudin 5 
    distance_mc_avgMIN_hist_M9V7Z1A1 = (distance_mc_sort_old(1,:));

    % This is the average value of the matrix above.
    distance_mc_avgMIN_M9V7Z1A1= mean(distance_mc_sort_old(1,:));

% take top 100 distances measured from this array. the first row indicates 
% the minimum distance measured for each HIEL to claudin5.
% Save 100 just in case it can be used in downstream analysis.
Total = 100;
distance_mc_sort_M9V7Z1A1 = distance_mc_sort_old(1:Total, 1:N_mc); 
clear distance_mc_sort_old; % delete original array
clear distance_mc_M9V7Z1A1;

end