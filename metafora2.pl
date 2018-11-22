#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Statistics::R;
use metagenomics_mod;

BEGIN {
	use Getopt::Long;
	die "You have an *old* version of Getopt::Long, please ",
		"update it asap!!\n"
		unless Getopt::Long->VERSION >= 2.4;
	use File::Basename;
	die "You have an *old* version of File::Basename, please ",
		"update it asap!!\n"
		unless File::Basename->VERSION >= 2.84;
	use Statistics::R;
	die "You have an *old* version of Statistics::R, please ",
		"update it asap!!\n"
		unless Statistics::R->VERSION >= 0.34;
}

my ($db,$jobs,$level,$mapping,$n,$s);

my $Usage = "Usage: perl metagenomic2.pl -m <Filename of mapping file> -db <database: gg, RDP or SILVA> -jobs <number of jobs for working in parallel> -level <taxonomic level>\n"
."-count <minimum number of count for retaining a OUT> -samples <minimum sample in which the OUT should be present>\n";

GetOptions(
	'db=s' => \$db,			 		# string - "gg", "RDP", "SILVA"
	'jobs=s'=>\$jobs,				# numeric - number of jobs for workin in parallel
	'm=s' => \$mapping,				# string - "name of mapping file without extension"
	'level=s' =>\$level,			# numeric - from "L1" to "L7"
	'count=s'=>\$n,					 # numeric - minimum number of count for retaining a OUT						
	'samples=s'=>\$s,				# numeric - minimum sample in which the OUT should be present
	) or die "$Usage\n";

my $path=$ENV{'PWD'};


## Validating mapping_file:
$mapping = join ("", "$mapping",".txt");
my $a = validate($mapping)."\n";
print "Mapping file validated\n";

my @reads = `find preprocessed_reads/`; 
my $dir_reads = dirname(@reads);
$s = spliting($dir_reads);
print "Generating a single file\n";

## Downloading database and editing parameters file:
if ($db eq "gg"){
	open (PAR,">params.txt");
		print PAR "pick_otus:enable_rev_strand_match	True\n";
		print PAR "alpha_diversity:metrics	observed_otus,observed_species,shannon,simpson,chao1\n";
		print PAR "plot_taxa_summary:chart_type	bar\n";
	close PAR;
	
	# Processing the data.
	# Step 1. Picking open references OTUs:
	# $pick = "params.txt";
	my $p = picking($jobs)."\n";
	print "Picking open-reference otu's completed\n"; 
	# Step 2. Removing chimeras (usearch61):
	my $chimeras = join ("","$db","_chimera_checking"); 
	my $c = chimeras($chimeras)."\n";
	print "Removing chimeras completed\n";
}

if ($db eq "RDP"){
	my $RDP = '';
	system ("wget https://www.mothur.org/w/images/d/dc/Trainset16_022016.rdp.tgz");
	system ("tar -xvzf Trainset16_022016.rdp.tgz");
	my @path = `find -name trainset16_022016.rdp.fasta`;
	$RDP = "@path";
	chomp $RDP;
	my @pathtax = `find -name trainset16_022016.rdp.tax`;
	my $RDPtax = "@pathtax";
	chomp $RDPtax;
	
	open (PAR,">params.txt");
		print PAR "pick_otus:enable_rev_strand_match	True\n";
		print PAR "pick_open_reference_otus:pick_otus_reference_seqs_fp $RDP\n";
		print PAR "assign_taxonomy:id_to_taxonomy_fp $RDPtax\n";
		print PAR "assign_taxonomy:reference_seqs_fp $RDP\n";
		print PAR "identify_chimeric_seqs:reference_seqs_fp	$RDP\n";
		print PAR "alpha_diversity:metrics	observed_otus,observed_species,shannon,simpson,chao1\n";
		print PAR "plot_taxa_summary:chart_type	area,bar,pie\n";
	close PAR;
	
	# Processing the data.
	# Step 1. Picking open references OTUs:
	my $REF = $RDP;
	my $p = picking_ref($REF,$jobs)."\n";
	print "Picking open-reference otu's completed\n";
	
	# Step 2. Removing chimeras (usearch61):
	my $chimeras = join ("","$db","_chimera_checking"); 
	my $c = chimeras($chimeras)."\n";
	print "Removing chimeras completed\n";
}

if ($db eq "SILVA"){
	my $SILVA = '';
	system ("wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_128_release.tgz");
	system ("tar -xvzf Silva_128_release.tgz");
	my @path2=`find -name 97_otus_16S.fasta `;
	$SILVA = "@path2";
	chomp $SILVA;
	my @path3 = `find ./SILVA_128_QIIME_release/taxonomy/16S_only/97 -name consensus_taxonomy_7_levels.txt `;
	my $SILVAtax = "@path3";
	chomp $SILVAtax;
		
	open (PAR,">params.txt");
		print PAR "pick_otus:enable_rev_strand_match	True\n";
		print PAR "pick_open_reference_otus:pick_otus_reference_seqs_fp $SILVA\n";
		print PAR "assign_taxonomy:id_to_taxonomy_fp $SILVAtax\n";
		print PAR "assign_taxonomy:reference_seqs_fp $SILVA\n";
		print PAR "identify_chimeric_seqs:reference_seqs_fp	$SILVA\n";
		print PAR "alpha_diversity:metrics	observed_otus,observed_species,shannon,simpson,chao1\n";
		print PAR "plot_taxa_summary:chart_type	area,bar,pie\n";
	close PAR; 

	# Processing the data.
	# Step 1. Picking open references OTUs:
	my $REF = $SILVA;
	my $pref = picking_ref($REF,$jobs)."\n";
	print "Picking open-reference otu's completed\n";
	
	# Step 2. Removing chimeras (usearch61):
	my $chimeras = join ("","$db","_chimera_checking"); 
	my $c = chimeras($chimeras)."\n";
	print "Removing chimeras completed\n";
}

## 2b. Filtering chimeras from files
my $chimeras = join ("","$db","_chimera_checking");	
my $f = filter_chimeras($chimeras)."\n";
print "Filtering chimeras from files completed\n";

# 2c. Making new phylogenetic tree without chimeras
my $newtree = "$path/4_pick_otus/rep_set_aligned_chimerafree_pfiltered.fasta";
my $ph = phylogeny($newtree)."\n";
print "Making new phylogenetic tree without chimeras\n";

# 2d. Making new OTU table without chimeras 
my $otu_table = join("","otu_table_", "$db", "_nochimera.biom"); #my $nm = join("","otu_table_", "$db", "_nochimera.biom");
my $o = otu_table($otu_table,$chimeras)."\n";
print "Making new OTU table without chimeras\n";

# Step 3. Summarizing Sample OTU Count's
my $stats= join("","otu_table_","$db","_nochimera_stats.txt");
my $sum = summarize($otu_table,$stats)."\n";

# Step 4. Taxonomic Composition Analysis
my $Ltax = join("","otu_table_","$db","_nochimera_","$level",".biom");
my $t = taxonomy($otu_table,$mapping,$Ltax,$level)."\n";
print "Taxa summary plots\n";

##########################
### Diversity analyses ###
##########################

## Step5. Alpha diversity analysis 
my $al = alpha($otu_table,$mapping)."\n";
print "Alpha diversity analysis calculated\n";

## Step6. Beta diversity analysis
my $bc= join("", "bray_curtis_otu_table_","$db","_nochimera.txt");
my $be = beta($otu_table,$mapping, $bc);
print "Beta diversity analysis";

## Clustering samples from matrix distances.
## Calling R script
my ($dist);
my $path_to_bdiv = "7_bdiv_plots/";
opendir(DIR2, $path_to_bdiv)|| die "can't opendir $path_to_bdiv: $!";
my @dm = readdir(DIR2);
foreach my $distance (@dm){
	if ($distance eq '.' or $distance eq '..') {
			next;
		}
	if ($distance=~/(dm.txt)|(bray_curtis)/){
		my $dist = $distance;
		# print "$dist\n";
		
	my $R = Statistics::R->new();
	$R->startR;
	# Here-doc with multiple R commands:
my $cmds=<<EOF; #unquoted for expanding $variable from script! (use \$ for non expanding it, Example:hr\$height)
library(dplyr)
library(purrr)

data<-read.table("$path_to_bdiv/$dist", row.names =1)
hr<-hclust(as.dist(data), method = "complete")
plot(hr)
abline(h= max(hr\$height/1.5),col = "blue") # a partir de casi la mitad del peso máximo del dendograma (ver figura) agrupa los terminales
## draw rectangules in the dendogram
clust <- rect.hclust(hr, h=max(hr\$height/1.5), border=0)
NumMemb = sapply(clust, length)
clust <- rect.hclust(hr, h=max(hr\$height/1.5), which=which(NumMemb>1))
## export info of cluster 
clust <- clust %>%
	setNames(paste0('', 1:length(clust)))
clust_df <- clust %>% 
	map2(.y = names(clust), 
		~ names(.x) %>% 
		as.data.frame() %>% 
		setNames("Sample") %>% 
		mutate(Group = .y)) %>% 
	bind_rows
		  
## add cluster_df to mapping file
mapping<-read.table("mapping_file.txt",sep="\t")
names(mapping)<-c("Sample","BarcodeSequence",	"LinkerPrimerSequence",	"Description")
clust_sort<-clust_df[order(match(clust_df[,1],mapping[,1])),]
mapping_cluster<-mapping %>% left_join(clust_sort, by = "Sample")
names(mapping_cluster)<-c("#SampleID","BarcodeSequence",	"LinkerPrimerSequence",	"Description", "Group")
mapping_cluster<-mapping_cluster [ , c(1,2,3,5,4)] #order columns by index
mapping_cluster[is.na(mapping_cluster)] <- " "
	
write.table(mapping_cluster, "mapp_cluster_$dist", quote = FALSE, row.names=FALSE, sep = "\t")

print('ok')
EOF
	my $out2 = $R->run($cmds);
	# print "$out2\n";
	$R->stopR() ;
	}		
}

## Validating new mapping_file post clusters from each distance matrix:
my $mapping_unweig = "mapp_cluster_unweighted_unifrac_dm.txt";
my $mapping_weig = "mapp_cluster_weighted_unifrac_dm.txt";
my $mapping_bc = join("","mapp_cluster_bray_curtis_otu_table_","$db","_nochimera.txt");
my $v1 = validate($mapping_unweig);
my $v2 =validate ($mapping_weig);
my $v3 =validate ($mapping_bc);
print "Validating new mapping_file post clusters from each distance matrix\n";

my $mak = making_plots($mapping_unweig,$mapping_weig,$mapping_bc,$bc);
print "Creating 2D PCoA Beta Diversity Plots from unweighted, weighted and Bray-Curtis distance matrix\n";
print "Making Distance Boxplots among clusters from unweighted, weighted and Bray-Curtis distance matrix\n";
	
## Step7. Abundance Significance 
my $ab = abundance($otu_table,$n,$s,$chimeras,$mapping,$mapping_weig);
print "Making abundance significance: Kruskal-wallis test among samples and among the clusters (weigthed unifrac distance\n";

my $plot = plots_abundance($mapping,$mapping_weig);
print "Making abundance heatmap and cytoscape networks from filtered otu table\n";

## Make a new phylogenetic tree after filtering otus table:
my $new_tree= `find -path $path -name otu_table_no_singletons.biom`;
my $nw = new_phylo($new_tree);
print "Making new filetered phylogenetic tree without chimeras and OTUs which are found in a determined percentage of samples \n";

## Step8. Correlation analysis: all & groups
my $ob = observation($otu_table,$mapping_weig);

## Converting .biom files to .txt format.
## declarar variable para la tabla de OTUS sin filtrar
my $rawOTU= "4_pick_otus/otu_table_mc2_w_tax_no_pynast_failures.biom";
## declarar variable para la tabla de OTUs sin chimeras
my $otu = "4_pick_otus/$otu_table";
## declarar variable para la tabla OTUS filtrada por parametros
my $filter_otu = "otu_table_no_singletons.biom";
my $cv = convert ($rawOTU);
my $cv2 = convert($otu);
my $cv3 = convert($filter_otu);
