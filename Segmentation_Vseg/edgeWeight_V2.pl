#!/usr/bin/perl -w
use POSIX;
use strict;
use List::MoreUtils qw(uniq);

my %pair;
my %edge;
my %weight;
my $dcutoff = 7;
my %index2coord;
my %index2umi;
my %index2mt;
my %index2pro;
my %index2intron;
my %index2totUMI;
my %index2totPro;
my %index2marker;
my %nodeFeature;
my %coord2index;
my %index2totMrk;


loadPair("OB36_nnAll.txt");
loadUMI("OB36suniq.coord");
loadFinal('OB36.final');
projectUMI();
edgeList();



sub edgeList{
	print 'index1',"\t",'index2',"\t","coordX_1","\t","coordY_1","\t","coordX_2","\t","coordY_2","\t","dist","\t","UMI1","\t","UMI2","\t","k1_UMI1","\t","k1_UMI2","\t","sharedGene","\t","sharedUMI","\t";
	print 'pUMI1',"\t",'pUMI2',"\t",'k1_pUMI1',"\t",'k1_pUMI2',"\t",'K1marker1',"\t",'K1marker2',"\n";
	foreach my $key(sort keys %pair){
		foreach my $sky(sort keys %{$pair{$key}}){
			if(not exists $index2totUMI{$key} or not exists $index2totUMI{$sky} or not exists $index2umi{$key} or not exists $index2umi{$sky}){
				next;
			}
			
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
			
			my $weight = $pair{$key}{$sky}/($index2umi{$key}*$index2umi{$sky}/$index2totUMI{$key} + $index2umi{$key}*$index2umi{$sky}/$index2totUMI{$sky});
			if($weight<=4){
				print $key,"\t",$sky,"\t",$index2coord{$key},"\t",$index2coord{$sky},"\t";
				print $pair{$key}{$sky},"\t",$index2umi{$key},"\t",$index2umi{$sky},"\t",$index2totUMI{$key},"\t",$index2totUMI{$sky};
				print "\t",$sharedFeature,"\t",$sharedCounts,"\t",$index2pro{$key},"\t",$index2pro{$sky},"\t",$index2totPro{$key},"\t",$index2totPro{$sky},"\t",$index2totMrk{$key},"\t",$index2totMrk{$sky},"\n";
			}
		}
	}
}

sub loadPair{
	my $file = shift;
	open(FF,"$file") || die;
	my %distK;
	while(my $line=<FF>){
		chomp $line;
		if($.==1){next;}
		$line=~s/\"//g;
		my @dat = split(/\,/,$line);
		my $parent = $dat[1];
		my @child = @dat[2..20];
		my @dist = @dat[22..40];
		my $a;
		my $b;
		# the 6th neighbor
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
				$edge{$b}{$a} = $dist[$i];
				$edge{$a}{$b} = $dist[$i];
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
		$coord2index{$dat[1]."\t".$dat[2]} = ($.-1);
		$index2coord{$.-1} = $dat[1]."\t".$dat[2];
		
		$index2umi{$coord2index{$dat[1]."\t".$dat[2]}} = $dat[3];
		$index2mt{$coord2index{$dat[1]."\t".$dat[2]}} = $dat[4];
		$index2intron{$coord2index{$dat[1]."\t".$dat[2]}} = $dat[5];
		$index2marker{$coord2index{$dat[1]."\t".$dat[2]}} = $dat[6];
		$index2pro{$coord2index{$dat[1]."\t".$dat[2]}} = $dat[7];

	}
	close(DD);
}


sub projectUMI{
	foreach my $key(sort keys %edge){
		
		my @level1 = keys %{$edge{$key}};
	
		#k1 level UMI and protein
		foreach my $inkey(@level1){
			if(exists $index2umi{$inkey}){
				$index2totUMI{$key} +=$index2umi{$inkey};
			} else {
				$index2totUMI{$key} +=0;
			}
		
			if(exists $index2pro{$inkey}){
				$index2totPro{$key} +=$index2pro{$inkey};
			} else {
				$index2totPro{$key} +=0;
			}
			
			if(exists $index2marker{$inkey}){
				$index2totMrk{$key} +=$index2marker{$inkey};
			} else {
				$index2totMrk{$key} +=0;
			}
		}
		#initialize
		@level1 = ();

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
