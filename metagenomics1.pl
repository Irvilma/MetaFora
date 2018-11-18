#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

BEGIN {
	use Getopt::Long;
	die "You have an *old* version of Getopt::Long, please ",
		"update it asap!!\n"
		unless Getopt::Long->VERSION >= 2.4;
	use File::Basename;
	die "You have an *old* version of File::Basename, please ",
		"update it asap!!\n"
		unless File::Basename->VERSION >= 2.84;
}

my $Usage = "Usage: perl metagenomic1.pl -t <data type, PE = Paired-End or SE = Single-End>\n"
."-q <Trimmomatic AVGQUAL value>\n"
."-tr <Trimmomatic TRAILING value>.\n";

my ($file,$file1,$file2,$type,$qual,$trailing);
my(@file1,@file2);

my $path=$ENV{'PWD'}; 
# print "$path\n";

GetOptions(
	't=s' => \$type,          # string - "PE" or "SE"
	'q=s' => \$qual,          # numeric
	'tr=s' => \$trailing,     # numeric
	) or die "$Usage\n";

my @fullname = `find ~/Reads/ -name *.fastq.gz`; # users have to create Reads directory with the files fastq.gz
my $dirfile = dirname(@fullname);
#print "$dirfile\n";
opendir(DIR, $dirfile)|| die "can't opendir $dirfile: $!";
my @files = readdir(DIR);

if ($type eq "SE"){
	mkdir "fastqc_SE";
	mkdir "preprocessed_reads";
	mkdir "fastqc_SE_post_edit";
	my @fastqc1 = `find ~/ -name fastqc`;
	my $fastqc = $fastqc1[0];
	chomp $fastqc;
	my @trimmomatic = `find ~/ -name trimmomatic-0.36.jar`;
	my $trimm = $trimmomatic[0];
	chomp $trimm;
	foreach $file (@files){
		if ($file=~/R1/){
			$file1=$file;
			push @file1, $file1;
			if(scalar(@fastqc1)>=0){
				system ("$fastqc $dirfile/$file1 -outdir fastqc_SE/");
			}
			else{
				print STDERR "FastQC is not in the path\n";
			exit;
			}
			if (scalar (@trimmomatic)<1){
				die "SEQCOUNT ERROR :: system args failed: $? : Is Trimmomatic-0.36 installed and exported to \$PATH ?";
			}
			system ("java -jar $trimm SE $dirfile/$file1  preprocessed_reads/$file1 AVGQUAL:$qual TRAILING:$trailing");
			system ("$fastqc preprocessed_reads/$file1 -outdir fastqc_SE_post_edit/");
		}
	}
}

if ($type eq "PE"){
	mkdir "fastqc_PE";
	mkdir "paired";
	mkdir "unpaired";
	mkdir "fastqc_PE_post_edit";
	mkdir "preprocessed_reads";
	my @fastqc1 = `find ~/ -name fastqc`;
	my $fastqc = $fastqc1[0];
	chomp $fastqc;
	my @trimmomatic = `find ~/ -name trimmomatic-0.36.jar`;
	my $trimm = $trimmomatic[0];
	chomp $trimm;
	my @merge = `find ~/ -name pear`;
	my $pear = $merge[0];
	chomp $pear;
	foreach $file (@files) {
		if ($file=~/R1/){
			my $file1=$file;
			push @file1,$file1;
			@file1= sort { lc($a) cmp lc($b) } @file1; #order the array by alphabetic order ( orden alfabético, ignorando mayúsculas y minúsculas)
			system ("$fastqc $dirfile/$file1 -outdir fastqc_PE/");
		}	
		if ($file=~/R2/){
			my $file2=$file;
			push @file2,$file2;
			@file2= sort { lc($a) cmp lc($b) } @file2;
			system ("$fastqc $dirfile/$file2 -outdir fastqc_PE/");
		}		
	 } 
	
	my $posicion=0;
	foreach my $R1 (@file1){
		system("java -jar $trimm PE $dirfile/$R1 $dirfile/$file2[$posicion] paired/trimm_$R1 unpaired/$R1 paired/trimm_$file2[$posicion] unpaired/$file2[$posicion] AVGQUAL:$qual"); # editing the reads with low quality with Trimmomatic
		system ("$fastqc paired/trimm_$R1 paired/trimm_$file2[$posicion] -outdir fastqc_PE_post_edit/");
		system ("$pear -f paired/trimm_$R1 -r paired/trimm_$file2[$posicion] -o preprocessed_reads/$R1"); #merge the pair-end reads
	$posicion++;
	}
	
	system ( "rm ./preprocessed_reads/*.discarded.fastq ./preprocessed_reads/*.unassembled.forward.fastq ./preprocessed_reads/*.unassembled.reverse.fastq" ); #remove unnecessary files
}
closedir(DIR);	
