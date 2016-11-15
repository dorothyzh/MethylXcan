##!/usr/bin/perl 
#use strict;
my $start_run = time();

$ex_probe_list = shift;   ###probes
$me_probe_annotation=shift;    ###M450_annotation.csv
$ex_dataset=shift;      ###MuTHER_Fat_normalized_31032010_uncondensed_Ids.txt
$me_dataset=shift;      ###MuTHER_Fat_450K_norm_AE_030913.txt
$methylation_annotation=shift;   ###Methylation_annotation.txt
$gene_annotation=shift;   ###Gene_annotation.txt
$path="run";

die "Usage: $0 ex_probe_list me_anno ex_dataset me_dataset methylation_annotation gene_annotation\n" unless $me_probe_annotation;


&run("mkdir run");
&run("head -n 1 $ex_dataset > ./data/ex.header.txt");
&run("head -n 1 $me_dataset > ./data/me.header.txt");
open(GL, $ex_probe_list);
while(<GL>)
{
    chomp;
    ($ex_probe, $gene, $NM, $loc) = split /\t/;

    &run("mkdir $path/$gene");
    &run("cat data/ex.header.txt > $path/$gene/ex.$ex_probe.txt");
    &run("grep $ex_probe $ex_dataset >> $path/$gene/ex.$ex_probe.txt");
    &run("head -1 $gene_annotation > $path/$gene/gene_anno.txt");
    &run("awk \'\$5==\"$gene\"\' $gene_annotation >> $path/$gene/gene_anno.txt");	
# methylation data
    if (not -e "$path/$gene/me.txt") {
        &run("head -1 $methylation_annotation > $path/$gene/cpg_anno.txt");
	&run("cat data/me.header.txt > $path/$gene/me.txt");
        &run("grep \",$gene,\" $me_probe_annotation | cut -f 1 -d \',\' > $path/$gene/me.probes.txt");
        open(PB, "$path/$gene/me.probes.txt");
        while(<PB>)
        {
            chomp;
            ($me_probe) = split;
            &run("grep $me_probe $me_dataset >> $path/$gene/me.txt ");            
            &run("grep $me_probe $methylation_annotation >> $path/$gene/cpg_anno.txt ");
            
        }
        close PB;
    }

    # do regression
    &run("module load R; Rscript --vanilla script/run_lm.r --args $path/$gene/ex.$ex_probe.txt $path/$gene/me.txt");
    &run("rm -rf $path/$gene");

}


sub run
{
    $cmd = shift;
    print "$cmd\n";
    system $cmd;
}


&run("cat run/* > final.txt");
&run("head -n 1 final.txt > ./MethylXcan.txt");
&run("cat final.txt | sed '/CpG/d' >> ./MethylXcan.txt");
&run("rm final.txt");




my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";
