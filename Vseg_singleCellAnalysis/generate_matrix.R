generate_matrix <- function(feature, region, algorithm, final_file, cut, work_dir){
	tmp_cmd = paste("perl ~/ShapeSegVnet/R/shapeImageSegVnet2.pl -e", feature, "-r",region, "-a",algorithm, "-f", final_file, "-c", cut, "-w", work_dir)
	system(tmp_cmd)
}