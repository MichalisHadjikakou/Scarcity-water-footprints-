function [All_TIMs] = Calculate_multipliers_modified(Total_output,Q,T)
%CALCULATES TOTAL IMPACT MULTIPLIERS FOR ANY NUMBER OF INDICATORS 
%   %   *AIM:* Standard total impact multipliers function - allows simplification of script 
%   *Authors:* Michalis Hadjikakou (Deakin University, Australia),
%   m.hadjikakou@unsw.edu.au
%   *Last updated:* 11/09/2017
%   *Data inputs: Total output, Q, T
% * Total output : Total output as row vector
% * Q : DIMs as column vector 
% * T : Transaction matrix  
% *Function outputs: *Total impact multipliers for any number of indicators
% (All_TIMs)

    Total_output = Total_output .* (Total_output>1);  % set small and negative values to zero;
    % this is necessary to prevent negative and huge DIMs (check with and without).    
    All_indicators = Q'; % Ensures that script can be used with many different indicators - Basically DIMs have already been calculated

    DIMs =  All_indicators; % can be done for several rows of DIMs - Columns become rows to ensure compatibility 

    DIMs(isinf(DIMs))=0 ;  % Set inf in DIMs to zero
    DIMs(isnan(DIMs))=0 ;  % Set nan in DIMs to zero
    DIMs = DIMs .* (DIMs>0);  % set negative values to zero!

    A = T ./ repmat(Total_output,length(Total_output),1); % This is identical to A = T * inv(diag(TotalIn));
    A(isinf(A))=0 ;  % Set inf in matrix A to zero
    A(isnan(A))=0 ;  % Set nan in matrix A to zero

    I = eye(size(A));
    L = inv(I-A);

    for indnum = 1:size(All_indicators,1)% Loop through indicator rows

        TIMs = DIMs(indnum,:) * L; % Calculating total impact multipliers for each indicator
        TIMs = TIMs .* (TIMs>0);  % set negative values to zero!
        All_TIMs(indnum,:) = TIMs; % Saving TIMs in matrix

    end


end

