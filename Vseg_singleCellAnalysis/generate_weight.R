generate_weight <- function(distance_file, coord_file, final_file, out_file){
	tmp_cmd = paste("perl ~/ShapeSegVnet/R/edgeWeight.pl -d", distance_file, "-c",coord_file, "-f", final_file, "-o", out_file)
	system(tmp_cmd)
}