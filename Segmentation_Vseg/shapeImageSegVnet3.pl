#!/usr/bin/perl -w
use POSIX;
use strict;
use Getopt::Long;
use List::Util qw(min);

my $usage = <<USAGE;

perl shapeImageSegVnet.pl <parameters>
    -feature|w 	  <str>	  the edge and feature file (xCoord,yCoord,feature1,feature2,feature3)
    -region|r 	  <str>	  Xlower,Xupper,ylower,yupper
    -algorithm|a  <str>	  fastGreedy or edgeBetweenness
    -final|f      <str>	  the raw final file
    -cut|c        <num>   cutoff of minimum barcodes per group

USAGE



my $featurefile = "";
my $region = "";
my $search = "";
my $lowercut = 30;
my $finalF = "";
GetOptions ("feature|w=s" => \$featurefile,
            "region|r=s" => \$region,
            "algorithm|a=s" => \$search,
            "finalF|f=s" => \$finalF,
            "cut|c=i" => \$lowercut
            );

if ($featurefile eq "" or $region eq "" or $search eq "" or $lowercut eq "" or $finalF eq "") {
    print $usage,"\n";
    exit;
}




####################################################assign barcode to cell segmentation

my $gridR = 800; #the subregion row size
my $gridC = 800; #the subregion col size
my %cellSeg;
my %cellBcount;
my %cell2blist;
my %featureUMI;
my %coord2UMI;
my $weightTitle;
vnetSeg($featurefile,$region);

cellIDsort();
outputSeg();


######################################################generate files for seurat umap analysis

my $gThresh = 0;
my $binThresh = 0;
my $filter = 1;
my $proteinOnly = 0;
my @segID2seq;
my $IDcounter = 0;

my %segCoord;

my %segCount;
my %segRecord;

my %geneInfor;
my %geneCount;

my %segBcount;
my %segID2Coord;
my %uniqueGeneMt;
my %uniqueGeneAl;

open(FF,"$finalF") || die;
my %uniqueBarcode;

while(my $line=<FF>){
	chomp $line;
	my @dat = split(/\t/,$line);
	
	if($dat[13]=~/rRNA/ or  $dat[13]=~/tRNA/) {next;}
	
	#protein filter
	if($proteinOnly == 1 and $dat[13] ne 'protein_coding'){
		next;
	}

	
	my $geneid = $dat[11];
	my $umi = $dat[4];
	my $spt = $dat[5];
	
	#convert the id to barcode as a formal input 
	if(not exists $uniqueBarcode{$spt}){
		$uniqueBarcode{$spt} = $IDcounter;
		$segID2seq[$IDcounter] = $spt;
		$IDcounter++;
	}
	
	$geneCount{$geneid}++;
	$geneInfor{$geneid} = $dat[12];
	
	#if($dat[2]<$xlimLeft or $dat[2]>$xlimRight or $dat[3]<$ylimLeft or $dat[3]>$ylimRight){
	#	next;
	#}
	
	my $xadj  = $dat[2];
	my $yadj  = 12000-$dat[3];
	if(not exists $cellSeg{$yadj}{$xadj}){
		next;
	}
	
	#bin total count stat
	$segCount{$cellSeg{$yadj}{$xadj}}++;
	
	#bin total barcode Count;
	my $keyCoord = $yadj.'_'.$xadj;
	
	if(not exists $segID2Coord{$cellSeg{$yadj}{$xadj}}{$keyCoord}){
		$segID2Coord{$cellSeg{$yadj}{$xadj}}{$keyCoord}=1;
		$segBcount{$cellSeg{$yadj}{$xadj}}++;
		#bin center calculate
		$segCoord{$cellSeg{$yadj}{$xadj}}{'xsum'} +=$dat[2];
		$segCoord{$cellSeg{$yadj}{$xadj}}{'ysum'} +=$dat[3];
	}
	
	

	
	#bin gene Information
	if(not exists $segRecord{$cellSeg{$yadj}{$xadj}}{$geneid}){
		if($dat[11]=~/^mt-/){
			$uniqueGeneMt{$cellSeg{$yadj}{$xadj}}++;
		}
		$uniqueGeneAl{$cellSeg{$yadj}{$xadj}}++;
	}
	$segRecord{$geneid}{$cellSeg{$yadj}{$xadj}}++;
}
close(FF);


print $IDcounter,"\n";


#filter lower expressed bin
#print bin coordinate

my $coordName = 'spt_3U_3_Mask'.'Coord.txt';
open(CD,">$coordName") || die;
print CD ',','xcoord',',','ycoord',',','binID',',','total_counts',',','BarcodeCount',"\n";
foreach my $fbinID(sort {$a<=>$b} keys %segCount){
	if($segCount{$fbinID}<$binThresh and $filter==1){next;}
	my $cy = $segCoord{$fbinID}{'ysum'}/$segBcount{$fbinID};
	my $cx = $segCoord{$fbinID}{'xsum'}/$segBcount{$fbinID};
	print CD $fbinID,',',$cx,',',$cy,',',$segID2seq[$fbinID],',',$segCount{$fbinID},',',$segBcount{$fbinID},"\n";
}
close(CD);

my $cgename = 'spt_3U_3_Mask'.'CellGeneExp.txt';

open(EF,">$cgename") || die;
#print header line
print EF 'bin';
foreach my $geneID(sort keys %geneInfor){
	if($geneCount{$geneID}<$gThresh){next;}
	print EF ',',$geneInfor{$geneID};
}
print EF "\n";

	
foreach my $efbinID(sort {$a<=>$b} keys %segCount){
	if($segCount{$efbinID}<$binThresh){next;}
	print EF $segID2seq[$efbinID];
	foreach my $geneID(sort keys %geneInfor){
		if($geneCount{$geneID}<$gThresh){next;}
		if(exists $segRecord{$geneID}{$efbinID}){
			print EF ',',$segRecord{$geneID}{$efbinID};
		} else {
			print EF ',',0;
		}
	}
	print EF "\n";
}
close(EF);


#loop over to have the segmentation 
sub vnetSeg{
	my $file = shift;
	my $area = shift;
	my ($xlimLeft,$xlimRight,$ylimLeft,$ylimRight) = split(/\,/,$area);
	
	my %gIndex2coord;
	my %gPair2Weight;
	my %node2node;
	my %nodeList;
	
	#my $record = 0;
    #loading postion and feature files
	open(FF,"$file") || die;
	while(my $ln=<FF>){
		chomp $ln;
		if($.==1){
			$weightTitle = $ln;
			next;
		}
		my @dat = split(/\t/,$ln);
		my $regionDetermine = 0;
		my ($x1,$y1);
		
		#any spots within the region are considered, including its connected spots
		if($dat[2]>$xlimLeft and $dat[2]<$xlimRight and $dat[3]>$ylimLeft and $dat[3]<$ylimRight){
			$x1=$dat[2];
			$y1=$dat[3];
			$regionDetermine++;
		}
		if($dat[4]>$xlimLeft and $dat[4]<$xlimRight and $dat[5]>$ylimLeft and $dat[5]<$ylimRight){
			$x1=$dat[4];
			$y1=$dat[5];
			$regionDetermine++;
		}
		if($regionDetermine==0){
			next;
		}
		
		#weight cutoff
		#my $base = (log2($dyd[13]*$dyd[8]/$dyd[9]+1)+log2($dyd[14]*$dyd[7]/$dyd[10]+1)+$dyd[12]);
		my $base = log2($dat[13]*$dat[8]/$dat[9]+1)+log2($dat[14]*$dat[7]/$dat[10]+1)+log2($dat[12]+1);
		#my $base = $dyd[13]*$dyd[8]/$dyd[9]+$dyd[14]*$dyd[7]/$dyd[10]+$dyd[12];
		#weight cutoff 
		if($base<=0){
			next;
		}
		my $weight = $dat[6]/$base;
		if($weight>3){
				next;
		}

		#$record++;
		$gIndex2coord{$dat[0]} = $dat[2].'_'.$dat[3];
		$gIndex2coord{$dat[1]} = $dat[4].'_'.$dat[5];
		$gPair2Weight{$dat[0]}{$dat[1]}{'U'} = $dat[6];
		$gPair2Weight{$dat[0]}{$dat[1]}{'P'} = $dat[7];
		$coord2UMI{$dat[2]}{$dat[3]} = $dat[8];
		$coord2UMI{$dat[4]}{$dat[5]} = $dat[9];
		#four possible combination
		my $a1 = int(floor($x1)/$gridR);
		my $b1 = int(floor($y1)/$gridC);

		$featureUMI{$a1.'_'.$b1}{$dat[0]} = $dat[7];
		$featureUMI{$a1.'_'.$b1}{$dat[1]} = $dat[8];
			
		$node2node{$dat[0]}{$dat[1]} = $ln."\t".$weight;
		$node2node{$dat[1]}{$dat[0]} = $ln."\t".$weight;

		$nodeList{$dat[0]} = $dat[7];
		$nodeList{$dat[1]} = $dat[8];
	}
	close(FF);
	#print $record,"\n";
	my $mergeCount = 0;
	my $edgeCount = 0;
	my $Label = 1;
	foreach my $key(sort keys %featureUMI){
		open(TT,'>regionU3Select.out') || die;
		print TT $weightTitle,"\t",'weight',"\n";
		my @vexID = sort {$featureUMI{$key}{$b} <=> $featureUMI{$key}{$a}} keys %{$featureUMI{$key}};
		my $half = int(($#vexID+1)/2);
		
		#choose half vetex based UMI order
		my %edgeSelect;
		for(my $i=0;$i<$half;$i++){
			my $kid = $vexID[$i];
			#level1
			foreach my $skid(keys %{$node2node{$kid}}){
				$edgeSelect{$kid}{$skid} = $node2node{$kid}{$skid};
				#level2
				foreach my $tkid(keys %{$node2node{$skid}}){
					$edgeSelect{$skid}{$tkid} = $node2node{$skid}{$tkid};
					#level3
					foreach my $fkid(keys %{$node2node{$tkid}}){
						$edgeSelect{$tkid}{$fkid} = $node2node{$tkid}{$fkid};
						#level4
						foreach my $vkid(keys %{$node2node{$fkid}}){
							$edgeSelect{$fkid}{$vkid} = $node2node{$fkid}{$vkid};
						}
					}
				}
			}
		}
		
		my %uniqueIdList;
		my %uniquePrint;
		foreach my $sky(sort keys %edgeSelect){
			foreach my $tky(sort keys %{$edgeSelect{$sky}}){
			#reassign label;
				my @dyd = split(/\t/,$edgeSelect{$sky}{$tky});
				#my $base = (log2($dyd[13]*$dyd[8]/$dyd[9]+1)+log2($dyd[14]*$dyd[7]/$dyd[10]+1)+$dyd[12]);
				#my $base = log2($dyd[13]*$dyd[8]/$dyd[9]+1)+log2($dyd[14]*$dyd[7]/$dyd[10]+1)+log2($dyd[12]+1);
				#my $base = $dyd[13]*$dyd[8]/$dyd[9]+$dyd[14]*$dyd[7]/$dyd[10]+$dyd[12];
				#weight cutoff 
				#if($base<=0){
				#	next;
				#}
				#my $weight = $dyd[6]/$base;
				#if($weight>2){
				#	next;
				#}
				
				if(not exists $uniqueIdList{$dyd[0]}){
					$uniqueIdList{$dyd[0]} = $Label;
					$Label++;
				}
				if(not exists $uniqueIdList{$dyd[1]}){
					$uniqueIdList{$dyd[1]} = $Label;
					$Label++;
				}
				if(not exists $uniquePrint{$uniqueIdList{$dyd[0]}}{$uniqueIdList{$dyd[1]}}){
					$uniquePrint{$uniqueIdList{$dyd[0]}}{$uniqueIdList{$dyd[1]}} = 1;
					print TT $uniqueIdList{$dyd[0]},"\t",$uniqueIdList{$dyd[1]},"\t",join("\t",@dyd[2..$#dyd]),"\n";
					$edgeCount++;
				}
			}
		}
		close(TT);
		print 'Edge:',$edgeCount,"\t";
		if($edgeCount>3000){
			system("Rscript --vanilla vNetseq3.R $key regionU3Select.out $search");
			maskMerge('cellAsign3.csv',$mergeCount);
			my $datestring = localtime();
			print $datestring,'> seg:     ',$mergeCount,"\t",$key,"\n";
		} else {
			print "skipped\n";
		}
		$mergeCount++;

		$edgeCount=0;
		$Label=1;
	}
}


#merge the mask from subregion into a whole matrix
sub maskMerge{
	my $file = shift;
	my $counter = shift;
	open(FF,"$file") || die;
	my %tempSeg;
	my %tempCount;
	my %tempList;
	while(my $ln=<FF>){
		if($.==1){next;}
		chomp $ln;
		$ln=~s/\"//g;
		my @dat = split(/\,/,$ln);
		if($counter==0){
			$cellSeg{12000-$dat[3]}{$dat[2]}=$dat[4].'|'.$dat[5];
			$cellBcount{$dat[4].'|'.$dat[5]}++;
			if(not exists $cell2blist{$dat[4].'|'.$dat[5]}){
				$cell2blist{$dat[4].'|'.$dat[5]} = (12000-$dat[2]).'|'.$dat[1];
			} else {
				$cell2blist{$dat[4].'|'.$dat[5]} = $cell2blist{$dat[4].'|'.$dat[5]}."\t".(12000-$dat[2]).'|'.$dat[1];
			}
			
		} else {
			$tempSeg{12000-$dat[3]}{$dat[2]}=$dat[4].'|'.$dat[5];
			$tempCount{$dat[4].'|'.$dat[5]}++;
			if(not exists $tempList{$dat[4].'|'.$dat[5]}){
				$tempList{$dat[4].'|'.$dat[5]} = (12000-$dat[3]).'|'.$dat[2];
			} else {
				$tempList{$dat[4].'|'.$dat[5]} = $tempList{$dat[4].'|'.$dat[5]}."\t".(12000-$dat[3]).'|'.$dat[2];
			}
		}
	}
	close(FF);
	
	#start merging two segmentation
	#rules: remove small, keep big, flush with new set
	if($counter==0){
		foreach my $key(sort keys %cellBcount){
			if($cellBcount{$key}<$lowercut){
				delete($cellBcount{$key});
				my @coordSet = split(/\t/,$cell2blist{$key});
				foreach my $coord(@coordSet){
					my($y,$x) = split(/\|/,$coord);
					delete($cellSeg{$y}{$x});
				}
				delete($cell2blist{$key});
			}
		}
	} else {
		foreach my $tkey(sort keys %tempCount){
			if($tempCount{$tkey}>=$lowercut){
				$cellBcount{$tkey} = $tempCount{$tkey};
				my @coordSet = split(/\t/,$tempList{$tkey});
				foreach my $coord(@coordSet){
					my($y,$x) = split(/\|/,$coord);
					$cellSeg{$y}{$x}=$tempSeg{$y}{$x};
				}
				$cell2blist{$tkey} = $tempList{$tkey};
			}
		}
	}
}


sub outputSeg{
	
	my %cellID;
	my $counter = 1;
	open(IDA, ">idAssignBarcode2.txt") || die;
	foreach my $key (sort keys %coord2UMI){
		foreach my $sky (sort keys %{$coord2UMI{$key}}){
			my ($v1,$v2) = ($key,$sky);
			if(exists $cellSeg{12000-$v2}{$v1}){
				if(exists $cellID{$cellSeg{12000-$v2}{$v1}}){
					print IDA $key,"\t",$sky,"\t",$cellSeg{12000-$v2}{$v1},"\t",$coord2UMI{$key}{$sky},"\n";
				} else {
					$cellID{$cellSeg{12000-$v2}{$v1}} = $counter;
					print IDA $key,"\t",$sky,"\t",$cellSeg{12000-$v2}{$v1},"\t",$coord2UMI{$key}{$sky},"\n";
					$counter++;
				}
			} else {
				print IDA $key,"\t",$sky,"\t",0,"\t",0,"\n";
			}
		}
	}
	close(IDA);
}


sub cellIDsort{
	my %cellid;
	my $idcounter = 1;
	foreach my $key(sort keys %cellSeg){
		foreach my $sky(sort keys %{$cellSeg{$key}}){
			if(not exists $cellid{$cellSeg{$key}{$sky}}){
				$cellid{$cellSeg{$key}{$sky}} = $idcounter;
				$cellSeg{$key}{$sky} = $idcounter;
				$idcounter++;
			} else {
				$cellSeg{$key}{$sky} = $cellid{$cellSeg{$key}{$sky}};
			}
		}
	}
}


sub log2 {
   my $n = shift;
   return log($n)/log(2);
}
