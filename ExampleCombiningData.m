%%The purpose of this script is to combine final data sets analyzed
%%obtained across several images that were analyzed with the "ExampleCode"
%%analysis template.


%Names of images included:
    %M2V1Z1A1
    %M2V2Z1A1
    %M2V2Z1A1
    %M2V2Z2A1
    %M2V2Z2A1
    %M2V2Z3A1
    %M2V2Z3A2
    %M3V3Z2
    %M3V3Z3A1
    %M3V3Z3A2
    %M3V3Z3A3
    %M3V3Z1A1
    %M3V3Z1A2
    %M3V3Z1A3
    %M3V1Z2A1
    %M4V1Z1A1
    %M4V1Z1A2
    %M4V2Z1A1
    

%SECTION 1 - BASIC INFORMATION
%How many images
    scale1 = 6.7448; %number of pixels per micron
    scale2 = 8.679;  %number of pixels per micron
    ImageNumber_3 = 22; 
   
%Number of HIEL (use length of Param)
IELHoles_3 = length(Param_M2V1Z1A1) + length(Param_M2V1Z1A2) +length(Param_M2V2Z1A1)+length(Param_M2V2Z1A2)+...
    length(Param_M2V2Z2A1)+ length(Param_M2V2Z3A1)+ length(Param_M2V2Z3A2)+ length(Param_M3V3Z2)+ length(Param_M3V3Z3A1)+...
    length(Param_M3V3Z3A2)+length(Param_M3V3Z3A3)+length(Param_M3V3Z1A1)+length(Param_M3V3Z1A2)+length(Param_M3V3Z1A3)+...
    length(Param_M3V1Z2A1)+ length(Param_M4V1Z1A1)+length(Param_M4V1Z1A2)+length(Param_M4V2Z1A1)+length(Param_M9V3Z1A1)+length(Param_M9V6Z1A1)+...
    length(Param_M12V2Z1A1)+length(Param_M16V1Z3A1);
IELHolesAvg_3 = IELHoles_3/ImageNumber_3;
    
%Average Radius per image
AvgRadius_3 = mean((mean(Param_M2V1Z1A1(:,3))/scale1) + (mean(Param_M2V1Z1A2(:,3))/scale1)+ (mean(Param_M2V2Z1A1(:,3))/scale1)+...
                    (mean(Param_M2V2Z1A2(:,3))/scale1) +(mean(Param_M2V2Z2A1(:,3))/scale1)+(mean(Param_M2V2Z3A1(:,3))/scale1)+...
                    (mean(Param_M2V2Z3A2(:,3))/scale1)+ (mean(Param_M3V3Z2(:,3))/scale2)+ (mean(Param_M3V3Z3A1(:,3))/scale2)+...
                    (mean(Param_M3V3Z3A2(:,3))/scale2)+(mean(Param_M3V3Z3A3(:,3))/scale2)+(mean(Param_M3V3Z1A1(:,3))/scale2)+...
                    (mean(Param_M3V3Z1A2(:,3))/scale2)+(mean(Param_M3V3Z1A3(:,3))/scale2)+(mean(Param_M3V1Z2A1(:,3))/scale2)+...
                    (mean(Param_M4V1Z1A1(:,3))/scale2)+(mean(Param_M4V1Z1A2(:,3))/scale2)+(mean(Param_M4V2Z1A1(:,3))/scale2)+...
                    (mean(Param_M9V6Z1A1(:,3))/scale2)+(mean(Param_M9V3Z1A1(:,3))/scale2)+(mean(Param_M12V2Z1A1(:,3))/scale2)+(mean(Param_M16V1Z3A1(:,3))/scale2));

%All radii for every HIEL detected
RadiiDist_3 = transpose (horzcat(transpose((Param_M2V1Z1A1(:,3)/scale1)),transpose((Param_M2V1Z1A2(:,3)/scale1)),transpose((Param_M2V2Z1A1(:,3)/scale1)),...
                    transpose((Param_M2V2Z1A2(:,3)/scale1)) , transpose((Param_M2V2Z2A1(:,3)/scale1)),transpose((Param_M2V2Z3A1(:,3)/scale1)),...
                    transpose((Param_M2V2Z3A2(:,3)/scale1)), transpose((Param_M3V3Z2(:,3)/scale2)), transpose((Param_M3V3Z3A1(:,3)/scale2)),...
                    transpose((Param_M3V3Z3A2(:,3)/scale2)),transpose((Param_M3V3Z3A3(:,3)/scale2)),transpose((Param_M3V3Z1A1(:,3)/scale2)),...
                    transpose((Param_M3V3Z1A2(:,3)/scale2)),transpose((Param_M3V3Z1A3(:,3)/scale2)),transpose((Param_M3V1Z2A1(:,3)/scale2)),...
                    transpose((Param_M4V1Z1A1(:,3)/scale2)),transpose((Param_M4V1Z1A2(:,3)/scale2)),transpose((Param_M4V2Z1A1(:,3)/scale2)),...
                    transpose((Param_M9V3Z1A1(:,3)/scale2)),transpose((Param_M9V6Z1A1(:,3)/scale2)),transpose((Param_M12V2Z1A1(:,3)/scale2)),...
                    transpose((Param_M16V1Z3A1(:,3)/scale2))));  
                
                
%Total Area Analyzed in microns
TotalArea_3 = ImageArea_M2V1Z1A1 + ImageArea_M2V1Z1A2 + ImageArea_M2V2Z1A1 + ImageArea_M2V2Z1A2+...
ImageArea_M2V2Z2A1+ImageArea_M2V2Z3A1+ImageArea_M2V2Z3A2+FullImage_Area_M3V3Z2+ImageArea_M3V3Z3A1+ImageArea_M3V3Z3A2+...
ImageArea_M3V3Z3A3+ImageArea_M3V3Z1A1+ImageArea_M3V3Z1A2+ImageArea_M3V3Z1A3+ImageArea_M3V1Z2A1+ImageArea_M4V1Z1A1+ImageArea_M4V1Z1A2+ImageArea_M4V2Z1A1+...
+ImageArea_M9V6Z1A1+ImageArea_M9V3Z1A1+ImageArea_M12V2Z1A1+ImageArea_M16V1Z3A1;

%Total Area by mouse
TotalArea_3_M2 = ImageArea_M2V1Z1A1 + ImageArea_M2V1Z1A2 +...
    ImageArea_M2V2Z1A1 + ImageArea_M2V2Z1A2 +ImageArea_M2V2Z2A1+ImageArea_M2V2Z3A1+ImageArea_M2V2Z3A2;
TotalVessels_M2 = 2;

TotalArea_3_M3 = FullImage_Area_M3V3Z2+ImageArea_M3V3Z3A1+ImageArea_M3V3Z3A2+...
    ImageArea_M3V3Z3A3+ImageArea_M3V3Z1A1+ImageArea_M3V3Z1A2+ImageArea_M3V3Z1A3 +ImageArea_M3V1Z2A1;
TotalVessels_M3 = 2;

TotalArea_3_M4 = ImageArea_M4V1Z1A1+ ImageArea_M4V1Z1A2 +ImageArea_M4V2Z1A1;
TotalVessels_M4 = 2;


TotalArea_3_M9 = ImageArea_M9V3Z1A1+ ImageArea_M9V6Z1A1;
TotalVessels_M9 = 2;

TotalArea_3_M12 = ImageArea_M12V2Z1A1;
TotalVessels_M12 = 1;


TotalArea_3_M16 = ImageArea_M16V1Z3A1;
TotalVessels_M16 = 1;

  
%Density of HIEL


Correct_HolesPerArea_M2V1Z1A1 = length(coord_M2V1Z1A1)/ImageArea_M2V1Z1A1;
Correct_HolesPerArea_M2V1Z1A2 = length(coord_M2V1Z1A2)/ImageArea_M2V1Z1A2;
Correct_HolesPerArea_M2V2Z1A1 = length(coord_M2V2Z1A1)/ImageArea_M2V2Z1A1;
Correct_HolesPerArea_M2V2Z1A2 = length(coord_M2V2Z1A2)/ImageArea_M2V2Z1A2;
Correct_HolesPerArea_M2V2Z2A1 = length(coord_M2V2Z2A1)/ImageArea_M2V2Z2A1;
Correct_HolesPerArea_M2V2Z3A1 = length(coord_M2V2Z3A1)/ImageArea_M2V2Z3A1;
Correct_HolesPerArea_M2V2Z3A2 = length(coord_M2V2Z3A2)/ImageArea_M2V2Z3A2;
Correct_HolesPerArea_M3V3Z2 = length(coord_M3V3Z2)/FullImage_Area_M3V3Z2;
Correct_HolesPerArea_M3V3Z3A1 = length(coord_M3V3Z3A1)/ImageArea_M3V3Z3A1;
Correct_HolesPerArea_M3V3Z3A2 = length(coord_M3V3Z3A2)/ImageArea_M3V3Z3A2;
Correct_HolesPerArea_M3V3Z3A3 = length(coord_M3V3Z3A3)/ImageArea_M3V3Z3A3;
Correct_HolesPerArea_M3V3Z1A1 = length(coord_M3V3Z1A1)/ImageArea_M3V3Z1A1;
Correct_HolesPerArea_M3V3Z1A2 = length(coord_M3V3Z1A2)/ImageArea_M3V3Z1A2;
Correct_HolesPerArea_M3V3Z1A3 = length(coord_M3V3Z1A3)/ImageArea_M3V3Z1A3;
Correct_HolesPerArea_M3V1Z2A1 = length(coord_M3V1Z2A1)/ImageArea_M3V1Z2A1;
Correct_HolesPerArea_M4V1Z1A1 = length(coord_M4V1Z1A1)/ImageArea_M4V1Z1A1;
Correct_HolesPerArea_M4V1Z1A2 = length(coord_M4V1Z1A2)/ImageArea_M4V1Z1A2;
Correct_HolesPerArea_M4V2Z1A1 = length(coord_M4V2Z1A1)/ImageArea_M4V2Z1A1;
Correct_HolesPerArea_M9V3Z1A1 = length(coord_M9V3Z1A1)/ImageArea_M9V3Z1A1;
Correct_HolesPerArea_M9V6Z1A1 = length(coord_M9V6Z1A1)/ImageArea_M9V6Z1A1;
Correct_HolesPerArea_M12V2Z1A1 = length(coord_M12V2Z1A1)/ImageArea_M12V2Z1A1;
Correct_HolesPerArea_M16V1Z3A1 = length(coord_M16V1Z3A1)/ImageArea_M16V1Z3A1;  


                
Correct_HolesPerArea_3_dist = transpose(horzcat(Correct_HolesPerArea_M2V1Z1A1,Correct_HolesPerArea_M2V1Z1A2,Correct_HolesPerArea_M2V2Z1A1,...
Correct_HolesPerArea_M2V2Z1A2,Correct_HolesPerArea_M2V2Z2A1,Correct_HolesPerArea_M2V2Z3A1,Correct_HolesPerArea_M2V2Z3A2,...
Correct_HolesPerArea_M3V3Z2,Correct_HolesPerArea_M3V3Z3A1,Correct_HolesPerArea_M3V3Z3A2,Correct_HolesPerArea_M3V3Z3A3,...
Correct_HolesPerArea_M3V3Z1A1,Correct_HolesPerArea_M3V3Z1A2,Correct_HolesPerArea_M3V3Z1A3,Correct_HolesPerArea_M3V1Z2A1,Correct_HolesPerArea_M4V1Z1A1,...
Correct_HolesPerArea_M4V1Z1A2,Correct_HolesPerArea_M4V2Z1A1,Correct_HolesPerArea_M9V3Z1A1,Correct_HolesPerArea_M9V6Z1A1,Correct_HolesPerArea_M12V2Z1A1,Correct_HolesPerArea_M16V1Z3A1)) ;   



HIEL_radi = (horzcat(transpose(Param_M2V1Z1A1(:,3)/scale1),...
transpose(Param_M2V1Z1A2(:,3)/scale1),...
transpose(Param_M2V2Z1A1(:,3)/scale1),...
transpose(Param_M2V2Z1A2(:,3)/scale1),...
transpose(Param_M2V2Z2A1(:,3)/scale1),...
transpose(Param_M2V2Z3A1(:,3)/scale1),...
transpose(Param_M2V2Z3A2(:,3)/scale1),...
transpose(Param_M3V3Z2(:,3)/scale2),...
transpose(Param_M3V3Z3A1(:,3)/scale2),...
transpose(Param_M3V3Z3A2(:,3)/scale2),...
transpose(Param_M3V3Z3A3(:,3)/scale2),...
transpose(Param_M3V3Z1A1(:,3)/scale2),...
transpose(Param_M3V3Z1A2(:,3)/scale2),...
transpose(Param_M3V3Z1A3(:,3)/scale2),...
transpose(Param_M3V1Z2A1(:,3)/scale2),...
transpose(Param_M4V1Z1A1(:,3)/scale2),...
transpose(Param_M4V1Z1A2(:,3)/scale2),...
transpose(Param_M4V2Z1A1(:,3)/scale2),...
transpose(Param_M9V3Z1A1(:,3)/scale2),...
transpose(Param_M9V6Z1A1(:,3)/scale2),...
transpose(Param_M12V2Z1A1(:,3)/scale2),...
transpose(Param_M16V1Z3A1(:,3)/scale2)));



HIEL_radi = transpose(HIEL_radi);
HIEL_radi_mean = mean(HIEL_radi);



%SECTION 2 - AVERAGE MINIMUM DISTANCE

%Average Minimum Distance of HIEL to Claudin5
AvgMinIEL_3 = mean((distance_avgMIN_M2V1Z1A1 + distance_avgMIN_M2V1Z1A2 + distance_avgMIN_M2V2Z1A1+...
    distance_avgMIN_M2V2Z1A2+ distance_avgMIN_M2V2Z2A1+ distance_avgMIN_M2V2Z3A1+ distance_avgMIN_M2V2Z3A2 + distance_avgMIN_M3V3Z2+...
    distance_avgMIN_M3V3Z3A1 + distance_avgMIN_M3V3Z3A2+ distance_avgMIN_M3V3Z3A3+ distance_avgMIN_M3V3Z1A1+ distance_avgMIN_M3V3Z1A2+...
    distance_avgMIN_M3V3Z1A3+ distance_avgMIN_M3V1Z2A1+distance_avgMIN_M4V1Z1A1+distance_avgMIN_M4V1Z1A2+distance_avgMIN_M4V2Z1A1+...
    +distance_avgMIN_M9V3Z1A1+distance_avgMIN_M9V6Z1A1+distance_avgMIN_M12V2Z1A1+distance_avgMIN_M16V1Z3A1));

%Average Minimum Distance Monte Carlo Simulations
AvgMinMC_3 = mean((distance_mc_avgMIN_M2V1Z1A1 + distance_mc_avgMIN_M2V1Z1A2 + distance_mc_avgMIN_M2V2Z1A1 +...
    distance_mc_avgMIN_M2V2Z1A2 +distance_mc_avgMIN_M2V2Z2A1+distance_mc_avgMIN_M2V2Z3A1+distance_mc_avgMIN_M2V2Z3A2+distance_mc_avgMIN_M3V3Z2+...
    distance_mc_avgMIN_M3V3Z3A1+distance_mc_avgMIN_M3V3Z3A2+distance_mc_avgMIN_M3V3Z3A3+distance_mc_avgMIN_M3V3Z1A1+distance_mc_avgMIN_M3V3Z1A2+...
    distance_mc_avgMIN_M3V3Z1A3+distance_mc_avgMIN_M3V1Z2A1+distance_mc_avgMIN_M4V1Z1A1+distance_mc_avgMIN_M4V1Z1A2+distance_mc_avgMIN_M4V2Z1A1+...
    distance_mc_avgMIN_M9V6Z1A1+distance_mc_avgMIN_M9V3Z1A1+distance_mc_avgMIN_M12V2Z1A1+distance_mc_avgMIN_M16V1Z3A1));

       
%Average Minimum Distance Positive Control Simulations
AvgMinPC_3 = mean((distance_pc_1um_avgMIN_M2V1Z1A1 + distance_pc_1um_avgMIN_M2V1Z1A2+ distance_pc_1um_avgMIN_M2V2Z1A1+...
    distance_pc_1um_avgMIN_M2V2Z1A2+distance_pc_1um_avgMIN_M2V2Z2A1+distance_pc_1um_avgMIN_M2V2Z3A1+distance_pc_1um_avgMIN_M2V2Z3A2+distance_pc_1um_avgMIN_M3V3Z2+...
    distance_pc_1um_avgMIN_M3V3Z3A1+distance_pc_1um_avgMIN_M3V3Z3A2+distance_pc_1um_avgMIN_M3V3Z3A3+distance_pc_1um_avgMIN_M3V3Z1A1+distance_pc_1um_avgMIN_M3V3Z1A2+...
    distance_pc_1um_avgMIN_M3V3Z1A3+distance_pc_1um_avgMIN_M3V1Z2A1+distance_pc_1um_avgMIN_M4V1Z1A1+distance_pc_1um_avgMIN_M4V1Z1A2+distance_pc_1um_avgMIN_M4V2Z1A1+...
    distance_pc_1um_avgMIN_M9V6Z1A1+distance_pc_1um_avgMIN_M9V3Z1A1+distance_pc_1um_avgMIN_M12V2Z1A1+distance_pc_1um_avgMIN_M16V1Z3A1));


%Average Minimum Distance Negative Control Simulations
AvgMinNC_3 = mean((distance_nc_avgMIN_M2V1Z1A1 + distance_nc_avgMIN_M2V1Z1A2 + distance_nc_avgMIN_M2V2Z1A1+...
    distance_nc_avgMIN_M2V2Z1A2+distance_nc_avgMIN_M2V2Z2A1+distance_nc_avgMIN_M2V2Z3A1+distance_nc_avgMIN_M2V2Z3A2+distance_nc_avgMIN_M3V3Z2+...
    distance_nc_avgMIN_M3V3Z3A1+distance_nc_avgMIN_M3V3Z3A2+distance_nc_avgMIN_M3V3Z3A3+distance_nc_avgMIN_M3V3Z1A1+distance_nc_avgMIN_M3V3Z1A2+...
    distance_nc_avgMIN_M3V3Z1A3+distance_nc_avgMIN_M3V1Z2A1+distance_nc_avgMIN_M4V1Z1A1+distance_nc_avgMIN_M4V1Z1A2+distance_nc_avgMIN_M4V2Z1A1+...
    distance_nc_avgMIN_M9V3Z1A1+distance_nc_avgMIN_M9V6Z1A1+distance_nc_avgMIN_M12V2Z1A1+distance_nc_avgMIN_M16V1Z3A1));



%SECTION 3 - SHORTEST DIST TO CLDN5 FOR EACH IEL HOLE (DISTRIBUTION)

%Minimum Distance HIEL to Cldn5 DISTRIBUTION
MinIELhist_3 = horzcat(distance_avgMIN_hist_M2V1Z1A1,distance_avgMIN_hist_M2V1Z1A2,distance_avgMIN_hist_M2V1Z1A2,...
    distance_avgMIN_hist_M2V2Z1A2,distance_avgMIN_hist_M2V2Z2A1,distance_avgMIN_hist_M2V2Z3A1,distance_avgMIN_hist_M2V2Z3A2,...
    distance_avgMIN_M3V3Z2_hist,distance_avgMIN_hist_M3V3Z3A1,distance_avgMIN_hist_M3V3Z3A2,distance_avgMIN_hist_M3V3Z3A3,...
    distance_avgMIN_hist_M3V3Z1A1,distance_avgMIN_hist_M3V3Z1A2,distance_avgMIN_hist_M3V3Z1A3,distance_avgMIN_hist_M3V1Z2A1,distance_avgMIN_hist_M4V1Z1A1,...
    distance_avgMIN_hist_M4V1Z1A2,distance_avgMIN_hist_M4V2Z1A1,distance_avgMIN_hist_M9V6Z1A1,distance_avgMIN_hist_M9V3Z1A1,distance_avgMIN_hist_M12V2Z1A1,...
    distance_avgMIN_hist_M16V1Z3A1);  
MinIELhist_3_length = length(MinIELhist_3);

%Minimum Distance HIEL to Cldn5 Monte Carlo DISTRIBUTION
MinMChist_3 = horzcat(distance_mc_avgMIN_hist_M2V1Z1A1,distance_mc_avgMIN_hist_M2V1Z1A2,distance_mc_avgMIN_hist_M2V1Z1A2,....
    distance_mc_avgMIN_hist_M2V2Z1A2,distance_mc_avgMIN_hist_M2V2Z2A1,distance_mc_avgMIN_hist_M2V2Z3A1,distance_mc_avgMIN_hist_M2V2Z3A2,...
    distance_mc_avgMIN_M3V3Z2_hist,distance_mc_avgMIN_hist_M3V3Z3A1,distance_mc_avgMIN_hist_M3V3Z3A2,distance_mc_avgMIN_hist_M3V3Z3A3,...
    distance_mc_avgMIN_hist_M3V3Z1A1,distance_mc_avgMIN_hist_M3V3Z1A2,distance_mc_avgMIN_hist_M3V3Z1A3,distance_mc_avgMIN_hist_M3V1Z2A1,distance_mc_avgMIN_hist_M4V1Z1A1,...
    distance_mc_avgMIN_hist_M4V1Z1A2,distance_mc_avgMIN_hist_M4V2Z1A1,distance_mc_avgMIN_hist_M9V6Z1A1,distance_mc_avgMIN_hist_M9V3Z1A1,...
    distance_mc_avgMIN_hist_M12V2Z1A1,distance_mc_avgMIN_hist_M16V1Z3A1) ;  
MinMChist_3_length = length(MinMChist_3);


%Minimum Distance HIEL to Cldn5 Positive Control DISTRIBUTION
MinPChist_3 = horzcat(distance_pc_1um_avgMIN_hist_M2V1Z1A1,distance_pc_1um_avgMIN_hist_M2V1Z1A2,distance_pc_1um_avgMIN_hist_M2V1Z1A2,...
    distance_pc_1um_avgMIN_hist_M2V2Z1A2,distance_pc_1um_avgMIN_hist_M2V2Z2A1,distance_pc_1um_avgMIN_hist_M2V2Z3A1,distance_pc_1um_avgMIN_hist_M2V2Z3A2,...
    distance_pc_1um_avgMIN_M3V3Z2_hist,distance_pc_1um_avgMIN_hist_M3V3Z3A1,distance_pc_1um_avgMIN_hist_M3V3Z3A2,distance_pc_1um_avgMIN_hist_M3V3Z3A3,...
    distance_pc_1um_avgMIN_hist_M3V3Z1A1,distance_pc_1um_avgMIN_hist_M3V3Z1A2,distance_pc_1um_avgMIN_hist_M3V3Z1A3,distance_pc_1um_avgMIN_hist_M3V1Z2A1,distance_pc_1um_avgMIN_hist_M4V1Z1A1,...
    distance_pc_1um_avgMIN_hist_M4V1Z1A2,distance_pc_1um_avgMIN_hist_M4V2Z1A1,distance_pc_1um_avgMIN_hist_M9V6Z1A1,distance_pc_1um_avgMIN_hist_M9V3Z1A1,...
    distance_pc_1um_avgMIN_hist_M12V2Z1A1,distance_pc_1um_avgMIN_hist_M16V1Z3A1) ;  
MinPChist_3_length = length(MinPChist_3);



%Minimum HIEL to Cldn5 Distance Negative Control DISTRIBUTION
MinNChist_3 = horzcat(distance_nc_avgMIN_hist_M2V1Z1A1, distance_nc_avgMIN_hist_M2V1Z1A2,distance_nc_avgMIN_hist_M2V1Z1A2,...
    distance_nc_avgMIN_hist_M2V2Z1A2,distance_nc_avgMIN_hist_M2V2Z2A1,distance_nc_avgMIN_hist_M2V2Z3A1,distance_nc_avgMIN_hist_M2V2Z3A2,...
    distance_nc_avgMIN_M3V3Z2_hist,distance_nc_avgMIN_hist_M3V3Z3A1,distance_nc_avgMIN_hist_M3V3Z3A2,distance_nc_avgMIN_hist_M3V3Z3A3,...
    distance_nc_avgMIN_hist_M3V3Z1A1,distance_nc_avgMIN_hist_M3V3Z1A2,distance_nc_avgMIN_hist_M3V3Z1A3,distance_nc_avgMIN_hist_M3V1Z2A1,distance_nc_avgMIN_hist_M4V1Z1A1,...
    distance_nc_avgMIN_hist_M4V1Z1A2,distance_nc_avgMIN_hist_M4V2Z1A1,distance_nc_avgMIN_hist_M9V6Z1A1,distance_nc_avgMIN_hist_M9V3Z1A1,...
    distance_nc_avgMIN_hist_M12V2Z1A1,distance_nc_avgMIN_hist_M16V1Z3A1) ;  
MinNChist_3_length = length(MinNChist_3);



%double check that all of the datasets are the same length.
%     if MinIELhist_3_length ~= MinMChist_3_length ~= MinPChist_3_length ~=MinPCEhist_3_length ~=MinNChist_3_length
%     
%         error('does not match')
%     else
%         print('good to go')
% 
%     end
% 
% 
%MinDist_3 = horzcat(transpose(MinPChist_3),transpose(MinMChist_3),transpose(MinNChist_3),transpose(MinIELhist_3), transpose(MinPCEhist_3));

MinDist_3 = horzcat(transpose(MinMChist_3),transpose(MinNChist_3),transpose(MinIELhist_3), transpose(MinPCEhist_3));


%%FIGURE Average minimum distance. 
AvgMinDistance = figure; 
edges = linspace(0,5,21);%beginning of range, end of range, # of bins+1
histogram(MinPChist_3, 'BinEdges', edges, 'FaceColor','b', 'Normalization', 'probability');
hold on;
histogram(MinNChist_3, 'BinEdges', edges, 'FaceColor','w','EdgeColor', 'k', 'Normalization', 'probability');
hold on;
histogram(MinMChist_3, 'BinEdges', edges, 'FaceColor','k', 'Normalization', 'probability');
hold on;
histogram(MinIELhist_3, 'BinEdges', edges, 'FaceColor','m', 'Normalization', 'probability');
hold on;
xlabel('Distance to Claudin 5 stain (um)', 'FontSize', 10);
ylabel('% Centers', 'FontSize', 10);
legend('Positive Control','Negative Control','Monte Carlo Simulations','3rd order artery','FontSize',12);
title('Minimum distance to interendothelial junction');
AvgMinDistance.InvertHardcopy = 'off'; %removes the default option of a white background so that the white claudin5_ coordinates show up correctly on the saved image
saveas(AvgMinDistance, '/Users/claireruddiman/AnalysisFiles/CombinedAnalysis3rd/AvgMinDistance.png');
hold off;

    
    

%SECTION 4 - DISTRIBUTION OF SHORTEST 100 DISTANCES (IEL TO CLDN5)

% Avg Top 100, averaged across all analyses IEL 
Top100IEL_3 = horzcat(distance_sort_M2V1Z1A1, distance_sort_M2V1Z1A2,distance_sort_M2V2Z1A1, distance_sort_M2V2Z1A2, distance_sort_M2V2Z2A1,...
    distance_sort_M2V2Z3A1,distance_sort_M2V2Z3A2,distance_sort_M3V3Z2,distance_sort_M3V3Z3A1,distance_sort_M3V3Z3A2,distance_sort_M3V3Z3A3,...
    distance_sort_M3V3Z1A1,distance_sort_M3V3Z1A2,distance_sort_M3V3Z1A3,distance_sort_M3V1Z2A1,distance_sort_M4V1Z1A1,distance_sort_M4V1Z1A2,distance_sort_M4V2Z1A1,...
    distance_sort_M9V6Z1A1,distance_sort_M9V3Z1A1,distance_sort_M12V2Z1A1,distance_sort_M16V1Z3A1);

% Avg Top 100, averaged across all analyses Monte Carlo
Top100MC_3 = horzcat(Dd_new_avg_M2V1Z1A1, Dd_new_avg_M2V1Z1A2, Dd_new_avg_M2V2Z1A1, Dd_new_avg_M2V2Z1A2, Dd_new_avg_M2V2Z2A1,...
    Dd_new_avg_M2V2Z3A1,Dd_new_avg_M2V2Z3A2,Dd_new_avg_M3V3Z2,Dd_new_avg_M3V3Z3A1,Dd_new_avg_M3V3Z3A2,Dd_new_avg_M3V3Z3A3,...
    Dd_new_avg_M3V3Z1A1,Dd_new_avg_M3V3Z1A2,Dd_new_avg_M3V3Z1A3,Dd_new_avg_M3V1Z2A1,Dd_new_avg_M4V1Z1A1,Dd_new_avg_M4V1Z1A2,Dd_new_avg_M4V2Z1A1,...
    Dd_new_avg_M9V6Z1A1,Dd_new_avg_M9V3Z1A1,Dd_new_avg_M12V2Z1A1,Dd_new_avg_M16V1Z3A1);

% Avg Top 100, averaged across all analyses Positive Control
Top100PC_3 =  horzcat(Dd_new_avg_pc_M2V1Z1A1,Dd_new_avg_pc_M2V1Z1A2,Dd_new_avg_pc_M2V2Z1A1, Dd_new_avg_pc_M2V2Z1A2,Dd_new_avg_pc_M2V2Z2A1,...
    Dd_new_avg_pc_M2V2Z3A1,Dd_new_avg_pc_M2V2Z3A2,Dd_new_avg_M3V3Z2_pc,Dd_new_avg_pc_M3V3Z3A1,Dd_new_avg_pc_M3V3Z3A2,Dd_new_avg_pc_M3V3Z3A3,...
    Dd_new_avg_pc_M3V3Z1A1,Dd_new_avg_pc_M3V3Z1A2,Dd_new_avg_pc_M3V3Z1A3,Dd_new_avg_pc_M3V1Z2A1,Dd_new_avg_pc_M4V1Z1A1,Dd_new_avg_pc_M4V1Z1A2,Dd_new_avg_pc_M4V2Z1A1,...
    Dd_new_avg_pc_M9V6Z1A1,Dd_new_avg_pc_M9V3Z1A1,Dd_new_avg_pc_M12V2Z1A1,Dd_new_avg_pc_M16V1Z3A1);

% Avg Top 100, averaged across all analyses negative control
Top100NC_3 =  horzcat(Dd_new_avg_nc_M2V1Z1A1,Dd_new_avg_nc_M2V1Z1A2,Dd_new_avg_nc_M2V2Z1A1, Dd_new_avg_nc_M2V2Z1A2,Dd_new_avg_nc_M2V2Z2A1,...
    Dd_new_avg_nc_M2V2Z3A1,Dd_new_avg_nc_M2V2Z3A2,Dd_new_avg_M3V3Z2_nc,Dd_new_avg_nc_M3V3Z3A1,Dd_new_avg_nc_M3V3Z3A2,Dd_new_avg_nc_M3V3Z3A3,...
    Dd_new_avg_nc_M3V3Z1A1,Dd_new_avg_nc_M3V3Z1A2,Dd_new_avg_nc_M3V3Z1A3,Dd_new_avg_nc_M3V1Z2A1,Dd_new_avg_nc_M4V1Z1A1,Dd_new_avg_nc_M4V1Z1A2,Dd_new_avg_nc_M4V2Z1A1,...
    Dd_new_avg_nc_M9V6Z1A1,Dd_new_avg_nc_M9V3Z1A1,Dd_new_avg_nc_M12V2Z1A1,Dd_new_avg_nc_M16V1Z3A1);


Top100Distances = figure; 
edges = linspace(0,6,31);%beginning of range, end of range, # of bins+1
histogram(Top100PC_3, 'BinEdges', edges, 'FaceColor','b', 'Normalization', 'probability');
hold on;
histogram(Top100NC_3, 'BinEdges', edges, 'FaceColor','w','EdgeColor', 'k', 'Normalization', 'probability');
hold on;
histogram(Top100MC_3, 'BinEdges', edges, 'FaceColor','k', 'Normalization', 'probability');
hold on;
histogram(Top100IEL_3, 'BinEdges', edges, 'FaceColor','m', 'Normalization', 'probability');
hold on;
xlabel('Distance to Claudin 5 stain (um)', 'FontSize', 10);
ylabel('% Centers', 'FontSize', 10);
legend('Positive Control','Negative Control','Monte Carlo Simulations','3rd order artery','FontSize',12);
title('Shortest 100 distances per IEL hole to claudin 5');
Top100Distances.InvertHardcopy = 'off'; %removes the default option of a white background so that the white claudin5_ coordinates show up correctly on the saved image
saveas(Top100Distances, '/Users/claireruddiman/AnalysisFiles/CombinedAnalysis3rd/Top100Distances.png');
hold off;



save('3rdorderRivanna.mat', 'IELHoles_3', 'IELHolesAvg_3','AvgRadius_3', 'AvgMinIEL_3','AvgMinMC_3','AvgMinPC_3','AvgMinPCE_3',...
'AvgMinNC_3','MinIELhist_3','MinMChist_3','MinPChist_3','MinPCEhist_3','MinNChist_3');

