#!/usr/bin/perl -w
use POSIX;
use Getopt::Long;
use strict;

my $usage = <<USAGE;

perl edgeWeigh.pl <parameters>
    -distance_file|d 	<str>	  the raw distance file
    -coord_file|c 	  	<str>	  the raw coord file
    -final_file|f  		<str>	  the raw final file
    -output_file|o 		<str>	  the output file name

USAGE

my $distance_file = "";
my $coord_file = "";
my $final_file = "";
my $output_file = "";
GetOptions ("distance_file|d=s" => \$distance_file,
            "coord_file|c=s" => \$coord_file,
            "final_file|f=s" => \$final_file,
            "output_file|o=s" => \$output_file
            );

if ($distance_file eq "" or $coord_file eq "" or $final_file eq "" or $output_file eq "") {
    print $usage,"\n";
    exit;
}

my %pair;
my %weight;
my $dcutoff = 7;
my %index2coord;
my %index2umi;
my %index2mt;
my %index2intron;
my %index2totUMI;
my %nodeFeature;
my %index2MYH7;
my %index2DCN;
my %index2APOD;
my %index2FABP4;
my %index2MALAT1;

loadPair($distance_file);
loadUMI($coord_file);
loadFinal($final_file);
projectUMI();
edgeList($output_file);



sub edgeList{
	my $output = shift;
	open(OUT,'>', $output) or die;
	foreach my $key(sort keys %pair){
		foreach my $sky(sort keys %{$pair{$key}}){
			print OUT $key,"\t",$sky,"\t",$index2coord{$key},"\t",$index2coord{$sky},"\t";
			
			#shared umi, feature, 
			my $sharedFeature = 0;
			my $sharedCounts = 0;
			foreach my $umiG(keys %{$nodeFeature{$index2coord{$key}}}){
				if(exists $nodeFeature{$index2coord{$sky}}{$umiG}){
					$sharedFeature++;
					if($nodeFeature{$index2coord{$sky}}{$umiG} < $nodeFeature{$index2coord{$key}}{$umiG}){
						$sharedCounts +=$nodeFeature{$index2coord{$sky}}{$umiG};
					} else {
						$sharedCounts +=$nodeFeature{$index2coord{$key}}{$umiG};
					}
				}
			}
			
			
			
			
			my $dist = $pair{$key}{$sky};
			my $vdist = $pair{$key}{$sky}/($index2umi{$key}*$index2umi{$sky}/$index2totUMI{$key}+$index2umi{$sky}*$index2umi{$key}/$index2totUMI{$sky});
			if(($index2MYH7{$key}+$index2MYH7{$sky})>0){
				$dist = $dist * 0.5;
				$vdist = $vdist * 0.5;
			}
			if(($index2DCN{$key}+$index2DCN{$sky})>0){
				$dist = $dist * 0.5;
				$vdist = $vdist * 0.5;
			}
			if(($index2APOD{$key}+$index2APOD{$sky})>0){
				$dist = $dist * 0.5;
				$vdist = $vdist * 0.5;
			}
			if(($index2FABP4{$key}+$index2FABP4{$sky})>0){
				$dist = $dist * 0.5;
				$vdist = $vdist * 0.5;
			}
			if(($index2MALAT1{$key}+$index2MALAT1{$sky})>0){
				$dist = $dist * 0.5;
				$vdist = $vdist * 0.5;
			}
			# my $tuDist = $pair{$key}{$sky}/($index2umi{$key}+$index2umi{$sky});
			# my $mtDist = 0;
			my $umiMinus = abs($index2umi{$key}*$index2umi{$sky}/$index2totUMI{$key}-$index2umi{$sky}*$index2umi{$key}/$index2totUMI{$sky});
			# if(($index2mt{$key}+$index2mt{$sky})>0){
			# 	$mtDist = $pair{$key}{$sky}/($index2mt{$key}+$index2mt{$sky});
			# }
			# my $intronDist = 0;
			# if(($index2intron{$key}+$index2intron{$sky})>0){
			# 	$intronDist = $pair{$key}{$sky}/($index2intron{$key}+$index2intron{$sky});
			# }
			
			print OUT $vdist,"\t",$dist,"\t",$index2umi{$key},"\t",$index2umi{$sky};
			#print "\t",0.09*$dist+0.0018*$vdist-0.3*$tuDist+0.0175*$mtDist-0.00045*$intronDist-0.42,"\n";
			#print "\t",$mtDist*0.035-$intronDist*0.00031+$sharedFeature*0.053-$sharedCounts*0.01+$umiMinus*0.025-0.0956;
			#print "\t",0.6145 - 0.02475 * $dist - 0.0004703 * $vdist + 0.07878 * $tuDist - 0.005235 * $mtDist + 0.00009689 * $intronDist;
			#print "\t",0.9141-$dist*0.01382-$vdist*0.0002693+$tuDist*0.04951-$mtDist*0.002391+$intronDist*0.00005748-$sharedFeature*0.003785+$sharedCounts*0.001343-$umiMinus*0.002403;
			my $cvdist = $pair{$key}{$sky}/($index2umi{$key}*$index2umi{$sky}/$index2totUMI{$key}+$index2umi{$sky}*$index2umi{$key}/$index2totUMI{$sky}+$sharedFeature);
			print OUT "\t",$cvdist;
			print OUT "\t",$sharedFeature,"\t",$sharedCounts;
			print OUT "\t",$umiMinus,"\n";
			
		}
	}
	close(OUT);
}

sub loadPair{
	my $file = shift;
	open(FF,"$file") || die;
	while(my $line=<FF>){
		chomp $line;
		if($.==1){next;}
		my @dat = split(/\,/,$line);
		my $parent = $dat[1];
		my @child = @dat[2..20];
		my @dist = @dat[22..40];
		my $a;
		my $b;
		
		for(my $i = 0;$i<=$#dist;$i++){
			if($dist[$i]<=$dcutoff){
				if($parent>$child[$i]){
					$a = $child[$i];
					$b = $parent;
				} else {
					$a = $parent;
					$b = $child[$i];
				}
				$pair{$a}{$b} = $dist[$i];
			} else {
				last;
			}
		}
	}
	close(FF);
}



#take the order as unique index to assign pairs
sub loadUMI{
	my $file = shift;
	open(DD,"$file") || die;
	while(my $ln=<DD>){
		chomp $ln;
		if($.==1){next;}
		my @dat = split(/\,/,$ln);
		$index2coord{$.-1} = $dat[1]."\t".$dat[2];
		$index2umi{$.-1} = $dat[3];
		$index2mt{$.-1} = $dat[4];
		$index2intron{$.-1} = $dat[5];
		$index2MYH7{$.-1} = $dat[4];
		$index2DCN{$.-1} = $dat[5];
		$index2APOD{$.-1} = $dat[6];
		$index2FABP4{$.-1} = $dat[7];
		$index2MALAT1{$.-1} = $dat[8];
	}
	close(DD);
}


sub projectUMI{
	foreach my $key(sort keys %pair){
		foreach my $sky(sort keys %{$pair{$key}}){
			$index2totUMI{$key} += $index2umi{$sky};
			$index2totUMI{$sky} += $index2umi{$key};
		}
	}
}


sub loadFinal{
	my $file = shift;
	#my $record = 0;
    #loading postion and feature files
	open(FF,"$file") || die;
	while(my $ln=<FF>){
		chomp $ln;
		my @dat = split(/\t/,$ln);
		my ($x1,$y1) = ($dat[2],$dat[3]);
		$nodeFeature{$x1."\t".$y1}{$dat[12]}++;
	}
	close(FF);
}
