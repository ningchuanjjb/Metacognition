% decodingDataSimplified.extraForMerged.tempMappingCellIndex(:,1)


decodingDataSimplified;
currentSessionIndex_B;


cellIndex_suite2p = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,1);


temp_range_B = FOVAllCellRange_multiFOV(currentSessionIndex_B,1):FOVAllCellRange_multiFOV(currentSessionIndex_B,2);
cellIndex_suite2p_B_raw = selevtivity_multiFOV.cellIndex_suite2p_multiFOV(temp_range_B);
cellIndex_suite2p_B = decodingDataSimplified.extraForMerged.tempMappingCellIndex_suite2p(:,end);
temp_range_AB = nan(1,length(cellIndex_suite2p_B));
for tempi=1:length(cellIndex_suite2p_B)
    temp_range_AB(tempi) = find(cellIndex_suite2p_B_raw==cellIndex_suite2p_B(tempi));
end
temp_range_AB = temp_range_AB + temp_range_B(1) - 1;

temp_cellIndexMapping_AB = temp_range_AB';