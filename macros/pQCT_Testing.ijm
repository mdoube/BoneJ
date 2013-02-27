
setBatchMode(true);
roiSelection = newArray("Bigger","Smaller","Left","Right","Top","Bottom","Central","Peripheral","SecondLargest","TwoLargestLeft","TwoLargestRight");
rotationSelection = newArray("According_to_Imax/Imin","Furthest_point","All_Bones_Imax/Imin","Not_selected_to_right","Selected_to_right");

for (i = 0; i < roiSelection.length; i++){
	roiSelectionOption = roiSelection[i];
	
	for (j = 0; j < roiSelection.length; j++){
		softRoiOption = roiSelection[j];

		for (k = 0; k < rotationSelection.length; k++){
			rotationOption = rotationSelection[k];
			print(roiSelectionOption, softRoiOption, rotationOption);
			run("Stratec pQCT", "select=C:/MyTemp/oma/Timon/tyo/Sulin2012/girls84m/I0020082.M01");
			run("Distribution Analysis", "  air_threshold=-40.0000 fat=40.0000 muscle_threshold=40.0000 marrow_threshold=80.0000 soft_tissue_threshold=150.0000 rotation_threshold=150.0000 area=280.0000 bmd=280.0000 roi_selection="+roiSelectionOption+" soft_tissue_roi_selection="+softRoiOption+" rotation_selection="+rotationOption+" analyse_mass_distribution analyse_concentric_density_distribution analyse_density_distribution analyse_soft_tissues manual_rotation_[+-_180_deg]=0.0000 save_visual_result_image_on_disk image_save_path=U:/programming/BoneJ/analysisTestCheck/"+roiSelectionOption+"_"+softRoiOption+"_"+substring(rotationOption,0,14)+"_");
//			wait(500);
			close();
			close();
		}
	}
}
print("pQCT test macro finished, all OK");
setBatchMode(false);
