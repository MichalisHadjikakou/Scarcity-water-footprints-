%% NATIONAL SCARCITY-ADJUSTED FOOTPRINTS FOR AUSTRALIA (RIDOUTT, HADJIKAKOU ET AL.)

% Author: Michalis Hadjikakou, UNSW
% Last updated: 16 March 2018
% Purpose of script: Matching of scarcity-adjusted water footprints and IElab IO
% classification, addition of RoW from Eora26, and calculation of total water use and scarcity-adjusted water multipliers for
% all agricultural and non-agricultural sectors
% Input dataset: 'SM1_All_input_datasets_revised'.xlsx' (last updated
% 16/03/2018), IO tables from IELab, 2013 Global IO tables (in basic prices) and water use extensions
% from Eora26 (http://www.worldmrio.com/simplified/)- full link: http://www.worldmrio.com/ComputationsM/Phase199/Loop082/simplified/Eora26_2013_bp.zip?auth=
% Output datasets: Total multipliers and footprint results for each
% agricultural commodity, uncertainty estimates and impact decomposition
% Used functions: AUS_RoW_disaggregation (allows augmenting AUS national table with Eora26 RoW, Calculate_multipliers_modified (performs IO calculus),TIM_decomposition (carried out decomposition of total impact multipliers) 

% Additional notes:  
%%
% 
% * Inflation-adjustment added using producer price indices
% * Figures/graphs - moved to seperate standalone scripts

%% 1. Reading in all input files and concordances and specifying initial parameters

% Add local folder/path here if required


% 1.1 Initialising parameters and loading IO table

quotient = 'AFLQ'; % Two different location quotients (SLQ or AFLQ)- Since this is a national table the differences are negligible
year = 2013; % Set MRIO year - can use any year in timeseries depending on input data
RoW = 'RoW_yes'; % RoW enabled or no
sec = 101; % Number of sectors
agsec = 26; % Number of agricultural sectors
RoWsec = 26; % Eora 26 RoW 
reg = 1; % Number or regions
diag_supply = 'diag_no'; % Diagonalize supply matrix - no is preferable since some industry supply more than one product; very small overall impact on results
uncertainty_analysis = 'yes'; % Option to carry out Monte Carlo analysis of multipliers based on historical range in water use intensities (from SM1)
write_results = 'no'; % Controls whether results are written as .csv files or not - change to yes if you want all results to be written as .csv files
chosen_indicators = [1:3,6]; % Eliminates AWARE-irri and AWARE non-irri where necessary since AWARE average encompasses both

if strcmp(quotient,'SLQ')

    T = csvread(['20170911_Mother_AllCountries_001_T-Results_',num2str(year),'_238_Markup001(full).csv']); % Transaction matrix
    v = csvread(['20170911_Mother_AllCountries_001_V-Results_',num2str(year),'_238_Markup001(full).csv']); % Value added matrix 

else 

    T = csvread(['20170911_Mother_AllCountries_001_T-Results_',num2str(year),'_249_Markup001(full).csv']); % Transaction matrix
    v = csvread(['20170911_Mother_AllCountries_001_V-Results_',num2str(year),'_249_Markup001(full).csv']); % Value added matrix
    
end

% 1.2 Loading water use and scarcity extensions

Ext_H2O = xlsread('SM1_All_input_datasets_revised.xlsx','Water_intensities_aggregated','E2:E102');
Ext_scarcity = xlsread('SM1_All_input_datasets_revised.xlsx','Water_intensities_aggregated','K2:O102');

All_ext = [Ext_H2O,Ext_scarcity]; % Adding all extensions together
All_ext = [All_ext;zeros(sec,size(All_ext,2))]; % Adding zeros for commodities

disp('All data inputs loaded successfully');

%% 2. Modifying IO table depending on user-specified functions (e.g RoW or diagonalising Supply matrix) and calculating TIMs

% 2.1 Fixing Supply matrix - see MRIO format (not necessary if
% co-production is of interest)
if strcmp(diag_supply,'diag_yes')
 
    Supply = T(1:sec,sec+1:2*sec); % 'Cutting out' NSW supply matrix for diagonalisation
    sum_NSW = sum(Supply,2); % Rowsums of supply matrix
    T(1:sec,sec+1:2*sec)= diag(sum_NSW); % Diagonalising original NSW supply matrix
    
end

if strcmp(RoW,'RoW_yes')
    
    %Aggregating RoW in IElab T-matrix to one single row
    RoW_row_sums = sum(T(2*sec+1:end,:)); % Calculating the row sum of the whole RoW block
    T_new = [T(1:2*sec,1:2*sec+1); RoW_row_sums]; % RoW consolidated as one extra row and column 
    
    %Importing text files containing all necessary files from Eora MRIO -
    %source: http://www.worldmrio.com/simplified/ - download link: http://www.worldmrio.com/ComputationsM/Phase199/Loop082/simplified/Eora26_2013_bp.zip?auth=
    T_Eora = importdata(['Eora26_' num2str(year) '_bp\' 'Eora26_' num2str(year) '_bp_T.txt']); % 4915 by 4915 MRIO transaction matrix
    FD_Eora = importdata(['Eora26_' num2str(year) '_bp\' 'Eora26_' num2str(year) '_bp_FD.txt']); % Final demand matrix
    VA_Eora = importdata(['Eora26_' num2str(year) '_bp\' 'Eora26_' num2str(year) '_bp_VA.txt']); % Value added matrix
    Q_Eora = importdata(['Eora26_' num2str(year) '_bp\' 'Eora26_' num2str(year) '_bp_Q.txt']); % Environmental extensions matrix
    labels_VA_Eora = importdata(['Eora26_' num2str(year) '_bp\' 'labels_VA.txt']); % Environmental extensions - labels
    labels_FD_Eora = importdata(['Eora26_' num2str(year) '_bp\' 'labels_FD.txt']); % Final demand labels
    labels_T_Eora = importdata(['Eora26_' num2str(year) '_bp\' 'labels_T.txt']); % Transaction matrix sector labels
    labels_Q_Eora = importdata(['Eora26_' num2str(year) '_bp\' 'labels_Q.txt']); % Transaction matrix sector labels
    
    % Running custom-function to obtain T-matrix with 26-sector RoW 
    
    [Dis_T,ext_RoW,Total_RoW] = AUS_RoW_disaggregation(T_Eora,Q_Eora,VA_Eora,labels_T_Eora,T_new);
    
    Sum_IE_Lab = sum(T_new)+sum(v); % Total sum of IE_Lab matrix plus value added
    Total_output = [Sum_IE_Lab(1:end-1) Total_RoW]; % Adding RoW totals to vector

    Index_Q = strcmp(strtrim(labels_Q_Eora),'Water use, blue water (m3)	Total'); % Selecting blue water in this case but can modify accordingly
    RoW_water = ext_RoW(Index_Q,:); % Selecting blue water extension
    RoW_water = repmat(RoW_water,6,1); % Creating 4 metrics (Blue water, WSI HH, WSI World and AWARE)
    AWARE_ag = 46; % See below
    AWARE_non_ag = 20; % See below

    RoW_water(4:6,1) = RoW_water(4:6,1)*AWARE_ag; % Multiplying by ag AWARE factor
    RoW_water(4:6,2:end)= RoW_water(4:6,2:end)*AWARE_non_ag; % Multiplying by non-ag factor 
    RoW_water(4,2:end)= 0; % AWARE irrigation
    RoW_water(5,1)=0;% AWARE non-irrigation

    % Regarding ROW we should use global average WSI value
    % 
    % For WSI WORLD EQ, the global average is 1
    % For AWARE: global average = 43
    % For WSI(HH_EQ): global average = 1
    % 
    % 
    % See more detail about AWARE
    % Interpretation – Relative to the world average 
    % It should be noted that a factor value of 1 is not equivalent to the factor for the average water consumption in the world, i.e. the world average factor to use when the location is not known. This value is calculated as the consumption-weighted average of the factors, which are based on 1/AMD and not AMD, hence the world consumption-based average has a value of 43 for unknown use and 20 and 46 respectively for non-agricultural and agricultural water consumption respectively.
    % 
    % In other words, if the location is unknown and the water use is non-agricultural then relevant AWAREW factor = 20. If the water use in unknown location is agricultural the n  factor = 46. The average for generic water use is 43

    Water_all = [All_ext;RoW_water']; % Combining AUS and RoW water extensions 

    % Modifying RoW extensions to ensure scarcity adjustments 
    Water_TIMs = Calculate_multipliers_modified(Total_output,Water_all,Dis_T);
    
    disp('Eora table and extensions loaded successfully');

else    
      %T = T(1:(reg+1)*sec,1:(reg+1)*sec); % Cutting out RoW
      %v = v(1:(reg+1)*sec); % Cutting out RoW
      T = [T(1:2*sec,:);sum(T(2*sec+1:end,:))];
      Water_all = [All_ext;zeros(1,size(All_ext,2))];
      Total_output = sum(T)+sum(v);  
      Water_TIMs = Calculate_multipliers_modified(Total_output,Water_all,T);
end


disp('All table modifications completed successfully')

% 2.2 Writing .csv file with final multiplier (L/$) results
if strcmp(write_results,'yes')

csvwrite(['TIMs_all_',date,'_',num2str(year),'_',quotient,'_',diag_supply,'_',diag_supply,RoW,'.csv'],Water_TIMs');

disp('All calculations completed successfully');

end

%% 3. Uncertainty analysis of average multipliers based on historical range (see SM1)

if strcmp(uncertainty_analysis,'yes')
   
   % 3.1 Initialising
   
   N = 1000; % Number of simulations 
   Unc_est_ag = xlsread('SM1_All_input_datasets_revised.xlsx','Sector details','F2:F27');
   Unc_est_other = 0.13; % See SM1 spreadsheet tab named Water_intensity_timeseries - 'Industry (2008-2015)'
   Direct_intensities = Water_all(:,chosen_indicators); % Selecting only AWARE AVERAGE to minimise compuational requirements
   Indicators = {'Water_use';'WSI_HH';'WSI_World';'AWARE_average'};
   
   Inputs_array = zeros(size(Direct_intensities,1),size(Direct_intensities,2),N); % Creating iputs array
   Results_array = zeros(size(Direct_intensities,1),size(Direct_intensities,2),N); % Creating results array
   Sect_unc = [Unc_est_ag;repmat(Unc_est_other,sec-agsec,1)]; % Uncertainty vector with margin for each sector
   b = Sect_unc+1; % Upper-bound of uncertainty
   a = 1-Sect_unc; % Lower bound of uncertainty
   %R = rand(sec,size(Direct_intensities,2),N); % Random 3D matrix to save results 
   
   Sect_original = Direct_intensities(sec+1:end,:); % These sectors remain unchanged - i.e no uncertainty assumption 
   
   % 3.2 Producing random distributions for direct water intensities -
   % this is carried out seperately for each indicator
   
   for i = 1:N % Number of iterations
   
       for j = 1:size(Direct_intensities,2)% Number of sectors
       
        Inputs_array(1:sec,j,i) = a.*Direct_intensities(1:sec,j) + (b.*Direct_intensities(1:sec,j)-a.*Direct_intensities(1:sec,j)).* rand(sec,1);
       
       end
       
       Inputs_array(sec+1:end,:,i)=Sect_original; %These sectors always have the same direct intensity    
   
   end
   
   % 3.3 Calculating TIMs with random distributions
   
   for i = 1:N % Number of iterations
   
       for j = 1:size(Direct_intensities,2)% Number of sectors
           
        Results_array(:,:,i) = Calculate_multipliers_modified(Total_output,Inputs_array(:,:,i),Dis_T)';
   
       end
       
   end

   % 3.4 Processing uncertainty results - this provides 3D matrices wtih
   
   stats = [5 25 50 75 95]; % Percentiles needed - in this case 5th,25th, 50th,75th and 95th  percentile
   Ag_results = zeros(agsec,N, size(Direct_intensities,2)); % Creating new empty array
   Ag_stats = zeros(agsec,length(stats),size(Direct_intensities,2)); % Creating another empty array to save stats
   
   for i = 1:size(Direct_intensities,2)
   
    Ag_results(:,:,i) = squeeze(Results_array(1:agsec,i,:)); % Restructuring matrix to consolidate results for each indicator
    Ag_stats(:,:,i) = prctile(Ag_results(:,:,i),stats,2);
   
    if strcmp(write_results,'yes')
    
       csvwrite(['TIMs_uncertainty_',date,'_',num2str(year),'_',quotient,'_',diag_supply,'_',diag_supply,RoW,'_indicator_',num2str(i),'.csv'],Ag_stats(:,:,i));
        
    end
    
   end
   
end


%% 4. Calculating final results for figures (Fig.2 and Fig.3 in Ridoutt, Hadjikakou et al.)

% 4.1 Collecting all necessary data and rearranging
% Order: WF, WSI HH, WSI world, AWARE-irrigation, AWARE-non-irrigation, AWARE average (last column)
Direct = Water_all; % All sectors direct
Total = Water_TIMs';
Indirect = Total-Direct; % All sectors indirect
Direct_perc = Direct./Total; % Percentage direct water use
Indirect_perc = Indirect./Total; % Percentage indirect water use

Direct_AG = Direct(1:agsec,:); % Agricultural sectors direct
Indirect_AG = Indirect(1:agsec,:); % Agricultural sectors indirect
Total_AG = Total(1:agsec,:); % Agricultural sectors total
Direct_perc_AG = Direct_perc(1:agsec,:); % Direct % agriculture
Indirect_perc_AG = Indirect_perc(1:agsec,:); % Indirect % agriculture

% Adjustments to account for head to kg conversions
Dollar_to_kg_factors = xlsread('SM1_All_input_datasets_revised.xlsx','Water volumes','K3:K28');
Dollar_to_kg_factors(1:20)= Dollar_to_kg_factors(1:20)/1000; % Assumes 1000 kg to a ton
Dollar_to_kg_factors(21)= Dollar_to_kg_factors(21)/5555; % Assumes 5555L of milk per hear - For 2013/2014, average milk production for aust was 5555 litres per cow. This is from Dairy Australia’s “Dairy in Focus” report.
Dollar_to_kg_factors(22)= Dollar_to_kg_factors(22)/mean([591,556]); % Assumes 585kg per hear (beef) - An MLA publication reports a projection for 2017 of 7.1 million head to be slaughtered and 2.1 million t carcase weight. If the dressing percentage was 50% (it can range a lot also), this would give a value of 591 kg/head. The same MLA publication refers to an average carcase wt in 2013 of 278kg, which with a similar 50% dressing percentage would yield 556 kg/head. 
Dollar_to_kg_factors(23)= Dollar_to_kg_factors(23)/(21.66/0.47); % Assumes this per head (sheep)- According to ABS7218, Australia produced 21,899.3 thousand lambs in 2013/2014 and lamb meat of 474,267 t, which equates to 21.66 kg per lamb. However, this will be dressed weight (after removal of skin, head, organs, etc.) during slaughter. For Australian lambs, the dressed % is usually about 47% (Ridoutt et al 2012 JCLP 28:127-133). As such, a best estimate of the weight of a lamb is 21.66/0.47 = 46 kg
Dollar_to_kg_factors(24:26)= Dollar_to_kg_factors(24:26)/1000; % Assumes 1000 kg to a ton
Dollar_to_kg_factors(24)= Dollar_to_kg_factors(24)*0.70; % 2009/2010 Chickens slaughtered:  465676.6 thousand head (ABS 7215) Meat (carcase) produced: 834409.0  t (ABS 7215)= 1.79 kg carcase/head - Carcase yield (dressing %) = 0.70 (RIRDC 2012)
Dollar_to_kg_factors(26) = Dollar_to_kg_factors(26)*0.76; % Pigs - carcase yield (dressing percentage) = 0.76 (RIRDC 2010)

% Writing .csv files with product footprint - optional
if strcmp(write_results,'yes')
    
    csvwrite(['Footprints_direct',date,'_',num2str(year),'_',quotient,'_',diag_supply,'_',diag_supply,RoW,'.csv'],Direct_AG.*Dollar_to_kg_factors);% Can also change to only AG
    csvwrite(['Footprints_total',date,'_',num2str(year),'_',quotient,'_',diag_supply,'_',diag_supply,RoW,'.csv'],Total_AG.*Dollar_to_kg_factors);% Can also change to only AG
end

disp('Product footprints successfully calculated');

% 4.2 Plotting each barplot - see outsourced scripts/functions

%% 5. Performing TIMs decomposition - See Wiedmann (2017) Journal of Economic Structures for theory (https://journalofeconomicstructures.springeropen.com/articles/10.1186/s40008-017-0072-0)

% 5.1 Calculating Leontief inverse and diagonalising multipliers for each
% extension
Decomp_sectors = 4; % Sector itself, other agriculture, other AUS industry/service, RoW 

[TIMs_decomposed,ag_TIMs_decomposed,ag_TIM_perc] = TIM_decomposition(Direct_perc_AG,Decomp_sectors,agsec,sec,Total_output,Water_all,Dis_T);

% 5.2 Writing .csv file with final multiplier results - this corresponds to
% the results in SM3
if strcmp(write_results,'yes')
   
   for i= 1:size(Water_all,2)
    
       csvwrite(['TIMs_decomposition_',date,'_',num2str(year),'_',quotient,'_',diag_supply,'_',diag_supply,RoW,'_indicator_',num2str(i),'.csv'],ag_TIM_perc(:,:,i));
   
   end
   
end

disp('All calculations and data export completed successfully');
















































