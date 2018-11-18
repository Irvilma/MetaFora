package metagenomics_mod;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Statistics::R;
use Exporter;
our @ISA=qw(Exporter);
our @EXPORT=qw(spliting validate picking picking_ref chimeras filter_chimeras phylogeny otu_table summarize taxonomy alpha beta making_plots abundance plots_abundance new_phylo observation convert);

sub spliting{
	my ($fna) = shift;
	my @split=`which multiple_split_libraries_fastq.py`;
	if (scalar(@split)<1){
		die "SEQCOUNT ERROR :: system args failed: $? : Is multiple_split_libraries_fastq.py installed and exported to \$PATH ?";
	}
	system ("multiple_split_libraries_fastq.py -i $fna -o 2_split_libraries --demultiplexing_method sampleid_by_file ");
}
sub validate{
	my ($map) = shift; 
	my $script = `which validate_mapping_file.py`;
	system ("validate_mapping_file.py -m $map -o 3_validate_map/");
}

sub picking{
	my ($params) = shift;
	my @pick=`which pick_open_reference_otus.py `;
	if(scalar(@pick)<1){
		die "SEQCOUNT ERROR :: system args failed: $? : Is pick_open_reference_otus.py installed and exported to \$PATH ?";
	}
	system ("pick_open_reference_otus.py -i 2_split_libraries/seqs.fna -o 4_pick_otus -p $params -a -O 8");
}
sub picking_ref{
	my ($params,$ref) = @_;
	my @pick=`which pick_open_reference_otus.py `;
	if(scalar(@pick)<1){
		die "SEQCOUNT ERROR :: system args failed: $? : Is pick_open_reference_otus.py installed and exported to \$PATH ?";
	}
	system ("pick_open_reference_otus.py -i 2_split_libraries/seqs.fna -o 4_pick_otus -p $params -r $ref -m uclust -a -O 8");
}
sub chimeras{
	my ($outChi) = shift;
	my @ch=`which identify_chimeric_seqs.py `;		
	if(scalar(@ch)<1){
		die "SEQCOUNT ERROR :: system args failed: $? : Is identify_chimeric_seqs.py installed and exported to \$PATH ?";
	}
	system("identify_chimeric_seqs.py -i 2_split_libraries/seqs.fna -m usearch61 -o $outChi/ --suppress_usearch61_ref");
}	
sub filter_chimeras{
	my ($newChi) = shift;
	my @filterfasta=`which filter_fasta.py `;		
	if(scalar(@filterfasta)<1){
		die "SEQCOUNT ERROR :: system args failed: $? : Is filter_fasta.py installed and exported to \$PATH ?";
	}
	system ("filter_fasta.py -f 4_pick_otus/pynast_aligned_seqs/rep_set_aligned_pfiltered.fasta -o 4_pick_otus/pynast_aligned_seqs/rep_set_aligned_chimerafree.fasta -s $newChi/chimeras.txt --negate");
	my @filteralign=`which filter_alignment.py `;	
	if(scalar(@filteralign)<1){
		die "SEQCOUNT ERROR :: system args failed: $? : Is filter_alignment.py installed and exported to \$PATH ?";
	}
	system ("filter_alignment.py -i 4_pick_otus/pynast_aligned_seqs/rep_set_aligned_chimerafree.fasta -o 4_pick_otus/");
}
sub phylogeny{
	my ($phy) = shift;
	my @phylo=`which make_phylogeny.py `;
	if(scalar(@phylo)<1){
		die "SEQCOUNT ERROR :: system args failed: $? : Is make_phylogeny.py installed and exported to \$PATH ?";
	}
	system ("make_phylogeny.py -i $phy -o 4_pick_otus/rep_set_chimerafree.tre");	
}
sub otu_table{
	my ($nm,$outChi) = @_;
	my @table= `which make_otu_table.py`;
	if (scalar(@table)<1){
		die "SEQCOUNT ERROR:: system args failed: $? : Is make_out_table.py isntalled and exported to \$PATH?";
	}
	system ("make_otu_table.py -i 4_pick_otus/final_otu_map_mc2.txt -o 4_pick_otus/$nm -t 4_pick_otus/uclust_assigned_taxonomy/rep_set_tax_assignments.txt -e $outChi/chimeras.txt");			
}
sub summarize{
	my ($table,$results) = @_;
	my @summarize=`which biom `;
	if (scalar(@summarize)<1){
		die "SEQCOUNT ERROR::system args failed: $?: IS biom summarize-table exported to \$PATH?";
	}
	system ("biom summarize-table -i 4_pick_otus/$table -o 4_pick_otus/$results");
}
sub taxonomy{
	my ($tab,$map,$tax,$lev) = @_;
	my @taxa=`which summarize_taxa_through_plots.py`;
	if (scalar(@taxa)<1){
		die "SEQCOUNT ERROR::system args failed: $? : Is summarize_taxa_through_plots.py exported to \$PATH?";
	}
	system ("summarize_taxa_through_plots.py -i 4_pick_otus/$tab -o 5_taxa_summary_plots -m $map -p params.txt");
	## Run summarize taxa, which will generate relative abundance text files for each level of phylogeny.
	my @sum_taxa=`which summarize_taxa.py`;
	if (scalar(@sum_taxa)<1){
		die "SEQCOUNT ERROR::system args failed: $?: Is summarize_taxa.py exported to \$PATH?";
	}
	system ("summarize_taxa.py -i 4_pick_otus/$tab -o 5b_summarize_taxa");
	system ("make_otu_heatmap.py -i 5b_summarize_taxa/$tax -o 5_taxa_summary_plots/otu_table_$lev\_heatmap.pdf -m $map");
}
sub alpha{
	my ($tab,$map) = @_;
	my @rarefaction=`which alpha_rarefaction.py `;
	if (scalar(@rarefaction)<1){
		die "SEQCOUNT ERROR::system args failed: $?: Is alpha_rarefaction.py exported to \$PATH?";
	}
	system ("alpha_rarefaction.py -i 4_pick_otus/$tab -o 6_alpha_output_folder -m $map -t 4_pick_otus/rep_set_chimerafree.tre -p params.txt");
}
sub beta{
	my($tab,$map,$bray) = @_;
	my @beta=`which beta_diversity_through_plots.py`;
	if (scalar(@beta)<1){
		die "SEQCOUNT ERROR::system args failed: $?: Is beta_diversity_through_plots.py exported to \$PATH?";
	}
	system ("beta_diversity_through_plots.py -i 4_pick_otus/$tab -o 7_bdiv_plots/ -m $map -t 4_pick_otus/rep_set_chimerafree.tre -p params.txt ");
	## Beta diversity non phylogenetic (Bray-curtis)
	my @beta_nonphy = `which beta_diversity.py`;
	if (scalar(@beta_nonphy)<1){
		die "SEQCOUNT ERROR::system args failed: $?: Is beta_diversity.py exported to \$PATH?";
	}
	system ("beta_diversity.py -i 4_pick_otus/$tab -m bray_curtis -o 7_bdiv_plots/");
	## make_emperor.py does not accept a distance matrix, so it is neccesary to run principal_coordinates.py on the distance matrix that you generated from beta_diversity.py
	my @coordinates = `which principal_coordinates.py`;
	if (scalar(@coordinates)<1){
		die "SEQCOUNT ERROR::system args failed: $?: Is principal_coordinates.py exported to \$PATH?";
	}
	system("principal_coordinates.py -i 7_bdiv_plots/$bray -o 7_bdiv_plots/bdiv_braycurtis_coords.txt");
}
sub making_plots{
	my ($map_un,$map_w,$map_bc,$bray) = @_;
	# Making 2D PCoA plots from distance matrices:
	my @plot2d=`which make_2d_plots.py`;
	if (scalar(@plot2d)<1){
		die "SEQCOUNT ERROR::system args failed: $?: Is make_2d_plots.py exported to \$PATH?";
	}
	## Unweighted
	system ("make_2d_plots.py -i 7_bdiv_plots/unweighted_unifrac_pc.txt -o 8_bdiv_2d_plot/ -m $map_un");
	## weighted distance matrices:
	system ("make_2d_plots.py -i 7_bdiv_plots/weighted_unifrac_pc.txt -o 8_bdiv_2d_plot/ -m $map_w");
	## Bray-Curtis distance matrices:
	system ("make_2d_plots.py -i 7_bdiv_plots/bdiv_braycurtis_coords.txt -o 8_bdiv_2d_plot/ -m $map_bc");
	## Making Distance Boxplots 
	my @boxplots=`which make_distance_boxplots.py`;
	if (scalar(@boxplots)<1){
		die "SEQCOUNT ERROR::system args failed: $?: Is make_distance_boxplots.py exported to \$PATH?";
	}
	## Unweighted
	system ("make_distance_boxplots.py -d 7_bdiv_plots/unweighted_unifrac_dm.txt -o 7_bdiv_plots/unweighted_distance_boxplot -m $map_un -f Group -n 999 --save_raw_data");
	## weighted distance matrices:
	system ("make_distance_boxplots.py -d 7_bdiv_plots/weighted_unifrac_dm.txt -o 7_bdiv_plots/weighted_distance_boxplot -m $map_w -f Group -n 999 --save_raw_data");
	## Bray-Curtis distance matrices:
	system ("make_distance_boxplots.py -d 7_bdiv_plots/$bray -o 7_bdiv_plots/bc_distance_boxplot -m $map_bc -f Group -n 999 --save_raw_data");
}
sub abundance{
	my ($tab,$count,$samples,$outChi,$map,$map_w) = @_;
	my @filter_otus = `which filter_otus_from_otu_table.py`;
	if (scalar(@filter_otus)<1){
		die "SEQCOUNT ERROR::system args failed: $? : Is filter_otus_from_otu_table.py exported to \$PATH?";
	}
	system ("filter_otus_from_otu_table.py -i 4_pick_otus/$tab -o otu_table_no_singletons.biom -n $count -s $samples -e $outChi/chimeras.txt");
	my @significance=`which group_significance.py`;
	if (scalar(@significance)<1){
		die "SEQCOUNT ERROR::system args failed: $? : Is group_significance.py exported to \$PATH?";
	}
	# comparing among samples (univariante analysis):
	system ("group_significance.py -i otu_table_no_singletons.biom -o kruskal_wallis_samples_test.txt -m $map -c Description -s kruskal_wallis");
	# comparing among clusters (univariante analysis):
	system ("group_significance.py -i otu_table_no_singletons.biom -o kruskal_wallis_cluster_weighted_test.txt -m $map_w -c Group -s kruskal_wallis");
}
sub new_phylo{
	my $tree = shift;
	system ("filter_fasta.py -f 4_pick_otus/pynast_aligned_seqs/rep_set_aligned.fasta -o 4_pick_otus/pynast_aligned_seqs/rep_set_aligned_chimerafree_filterotus.fasta -b otu_table_no_singletons.biom");
	my @filter_tree = `which filter_tree.py`;
	if (scalar(@filter_tree)<1){
		die "SEQCOUNT ERROR::system args failed: $? : Is filter_tree.py exported to \$PATH?";
	} 
	system ("filter_tree.py -i 4_pick_otus/rep_set_chimerafree.tre -o 4_pick_otus/rep_set_chimerafree_filterotus.tre -f 4_pick_otus/pynast_aligned_seqs/rep_set_aligned_chimerafree_filterotus.fasta");
}
sub plots_abundance{
	my ($map,$map_w)= @_;
	## Abundance Heatmap
	my @heatmap=`which make_otu_heatmap.py`;
	if (scalar(@heatmap)<1){
		die "SEQCOUNT ERROR::system args failed: $? : Is make_otu_heatmap.py exported to \$PATH?";
	}
	system ("make_otu_heatmap.py -i otu_table_no_singletons.biom -o otu_table_no_singletons_heatmap.pdf -m $map -t 4_pick_otus/rep_set_chimerafree_filterotus.tre");	

	## Making Cytoscape Networks: This script generates the otu network files to be passed into cytoscape and statistics for those networks. It uses the OTU fileand the user metadata mapping file.
	my @network=`which make_otu_network.py`;
	if (scalar(@network)<1){
		die "SEQCOUNT ERROR::system args falied: $?: Is make_otu_network.py exported to \$PATH?";
	}
	system ("make_otu_network.py -i otu_table_no_singletons.biom -m $map_w -o 9_otu_network");
}
sub observation{
	my($tab,$map_w) = @_;
	my @correlation = `which observation_metadata_correlation.py`;
	if (scalar (@correlation)<1){
		die "SEQCOUNT ERROR::system args falied: $?: Is observation_metadata_correlation.py exported to \$PATH?";
	}
	#correlation among all samples:
	system ("observation_metadata_correlation.py -i 4_pick_otus/$tab -m $map_w -c Description -s spearman -o spearman_otu_gradient_samples.txt");
	#correlation among clusters:
	system ("observation_metadata_correlation.py -i 4_pick_otus/$tab -m $map_w -c Group -s spearman -o spearman_otu_gradient_clusters.txt");
}
sub convert{
	my ($tab) = shift;
	my @convertion = `which biom`;
	if (scalar (@convertion)<1){
		die "SEQCOUNT ERROR::system args falied: $?: Is biom exported to \$PATH?";
	}
	## Convert .biom files to .txt 
	system("biom convert -i $tab -o $tab.txt --to-tsv");

}

1;