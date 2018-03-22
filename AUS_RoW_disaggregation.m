function [Dis_T,All_DIMs_Eora,Total_RoW] = AUS_RoW_disaggregation(T_Eora,Q,VA,Eora_labels,T_IElab)
%ADDS an RoW block to form a two-region AUS-RoW MRIO - can add more regions
%if desired by adding more countries
%   %   *AIM:* Adds RoW region 
%   *Authors:* Michalis Hadjikakou (Deakin University, Australia)
%   *Last updated:* 09/11/2017
%   *Data inputs:*  
% * n : number of sectors for RoW
% * labels : Country name labels from Eora 26-sector full MRIO
% * T_Eora : Eora full transaction matrix with all countries
% * T_IElab : IElab full transaction matrix for Australia
% *Function outputs: *Transaction matrix with disaggregated RoW region
% (with n sectors)

% Initialising 

currency_rate = 0.967915; % US$/AU$ for 2013 - modify depending on year of interest
n = 26; % number of sectors for each country - this is used for indexing purposes (code remains flexible to changed MRIO sector resolution)
AUS_RoW_conc = zeros(2*n+1,length(T_Eora)); % Dimensions are compatible with 2 region plus residual rest of RoW column

% Converting from 1000 USD to million AUD

T_Eora = (T_Eora/currency_rate)/1000;
Q = (Q/currency_rate)/1000;
VA = (VA/currency_rate)/1000;

% Locating countries

Index_AUS = find(not(cellfun('isempty', strfind(Eora_labels, 'Australia')))); % Finding positions of Australian sectors in original matrix
%Index_USA = find(not(cellfun('isempty', strfind(Eora_labels, 'USA')))); % Finding positions of Australian sectors in original matrix
%Index_CHN = find(not(cellfun('isempty', strfind(Eora_labels, 'China')))); % Finding positions of Australian sectors in original matrix
AUS_RoW_conc(1:n,Index_AUS) = eye(n); % Filling in Australian sector positions with ones
%AUS_USA_CHN_RoW_conc(n+1:2*n,Index_USA) = eye(n); % Filling in USA sector positions with ones
%AUS_USA_CHN_RoW_conc(2*n+1:3*n,Index_CHN) = eye(n); % Filling in China sector positions with ones

RoW_before_Australia = (1:Index_AUS(1)-1)'; % Index of all countries before Australia 
RoW_after_Australia = ((Index_AUS(end)+1):length(T_Eora)-1)'; % Index of all countries after Australia 

AUS_RoW_conc(1:n,Index_AUS(1):Index_AUS(end)) = eye(n); % Australia sector aggregation
AUS_RoW_conc(n+1:2*n,1:Index_AUS(1)-1) = repmat(eye(n),1,(Index_AUS(1)-1)/n); % Aggregating all RoW countries before Australia (alphabetically)
AUS_RoW_conc(n+1:2*n,Index_AUS(end)+1:length(T_Eora)-1) = repmat(eye(n),1,length(RoW_after_Australia)/n); % Aggregating all RoW countries between USA and RoW 
AUS_RoW_conc(end,end) = 1; % RoW kept as one row in concordance but this can later be removed as rest of RoW is negligible 

Eora_26_simplified = AUS_RoW_conc*T_Eora*AUS_RoW_conc'; % Creating aggregated MRIO
Eora_26_simplified(end,:) = []; % Deleting last row as all inputs have now been absorbed into a single RoW region
Eora_26_simplified(:,end) = []; % Deleting last column as all exports to rest of RoW have now been absorbed into a single RoW region

%RoW_input_shares = Eora_26_simplified(3*n+1:4*n,1:4*n)./repmat(sum(Eora_26_simplified(3*n+1:4*n,1:4*n)),n,1); % Determining RoW input shares from rest of RoW
%RoW_inputs_rest_of_RoW = repmat(Eora_26_simplified(end,1:end-1),n,1).*RoW_input_shares; % Calculating rest of RoW inputs to be added to RoW inputs
%RoW_corrected = Eora_26_simplified(3*n+1:4*n,1:4*n)+RoW_inputs_rest_of_RoW; % Correcting RoW inputs to include all rest of RoW
%Eora_26_simplified(3*n+1:4*n,1:4*n)= RoW_corrected; % Indexing original matrix and correcting inputs
%Eora_26_simplified(end,:) = []; % Deleting last row as all inputs have now been absorbed into a single RoW region

%RoW_output_shares = Eora_26_simplified(1:4*n,3*n+1:4*n)./repmat(sum(Eora_26_simplified(1:4*n,3*n+1:4*n),2),1,n); % Determining RoW input shares from rest of RoW
%RoW_outputs_rest_of_RoW = repmat(Eora_26_simplified(1:end,end),1,n).*RoW_output_shares;% Calculating rest of RoW exports to be added to RoW inputs
%RoW_exports_corrected = Eora_26_simplified(1:4*n,3*n+1:4*n)+RoW_outputs_rest_of_RoW;% Correcting RoW exports to include all rest of RoW 
%Eora_26_simplified(1:4*n,3*n+1:4*n)= RoW_exports_corrected;% Indexing original matrix and correcting 
%Eora_26_simplified(:,end) = []; % Deleting last column as all exports to rest of RoW have now been absorbed into a single RoW region

Output_totals_RoW = sum(Eora_26_simplified(n+1:2*n,:),2); % Calculating row sums of all RoW exports and domestic production - these will act as weights for initial disaggregation
Output_weights_RoW = Output_totals_RoW'/sum(Output_totals_RoW); % Calculating weights 
Dis_T = disaggregate_USE(T_IElab,length(T_IElab),Output_weights_RoW); % Disaggregating using Hadjikakou's disaggregate function (see working directory)
Dis_T(length(Dis_T)-25:length(Dis_T),length(Dis_T)-25:length(Dis_T)) = Eora_26_simplified(n+1:2*n,n+1:2*n);% Setting RoW to be equal to RoW in Eora 

%New_T_label = [Eora_MRIO.Tlabels xlsread([datapath 'Eora_26/' 'labels_T.xlsx'],'labels_T','D1:D16')]; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% 2.3 Disaggregating RoW environmental extensions, VA and calculating DIMs
Ext_RoW = zeros(size(Q,1),size(AUS_RoW_conc,1)); % Empty matrix where aggregated extensions are to be saved

for ext = 1:size(Q,1)
    Ext_RoW(ext,:) = AUS_RoW_conc*Q(ext,:)'; % Looping through each extension in Q matrix and multiplying by concordance matrix
end

VA_agg = AUS_RoW_conc*VA'; % Aggregating original final demand
%VA_weights_RoW = VA_agg(3*n+1:end-1,:)./repmat(sum(VA_agg(3*n+1:end-1,:)),size(VA_agg(3*n+1:end-1,:),1),1); % Calculating VA weights for RoW

Sum_VA = sum(VA_agg,2); % Calculating sum of Value Added

Total_RoW = sum(Eora_26_simplified(n+1:2*n,n+1:2*n),1)+Sum_VA(end-n+1:end)'; % Sum of RoW T matrix plus RoW VA

All_DIMs_Eora = Ext_RoW(:,end-n:end-1)./repmat(Total_RoW,length(Ext_RoW),1); % Calculating direct impact multipliers for RoW that can be used later on in the disaggregated model

% Y_agg = AUS_USA_CHN_RoW_conc*FD;

end

