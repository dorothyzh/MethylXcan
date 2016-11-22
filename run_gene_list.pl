#!/usr/bin/perl 
#use strict;
my $start_run = time();

$ex_probe_list = shift;   ###probes
$ex_dataset=shift;      ###MuTHER_Fat_normalized_31032010_uncondensed_Ids.txt
$me_dataset=shift;      ###MuTHER_Fat_450K_norm_AE_030913.txt
$methylation_annotation=shift;   ###Methylation_annotation.txt
$gene_annotation=shift;   ###Gene_annotation.txt
$path="run";

die "Usage: $0 ex_probe_list ex_dataset me_dataset methylation_annotation gene_annotation\n" unless $ex_probe_list $ex_dataset $me_dataset $methylation_annotation $gene_annotation;


system("mkdir run");
system("head -n 1 $ex_dataset > ./data/ex.header.txt");
system("head -n 1 $me_dataset > ./data/me.header.txt");
open(GL, $ex_probe_list);
while(<GL>)
{
    chomp;
    ($ex_probe, $gene, $NM, $loc) = split /\t/;

    system("mkdir $path/$gene");
    system("cat data/ex.header.txt > $path/$gene/ex.$ex_probe.txt");
    system("grep $ex_probe $ex_dataset >> $path/$gene/ex.$ex_probe.txt");
    system("head -1 $gene_annotation > $path/$gene/gene_anno.txt");
    system("awk \'\$5==\"$gene\"\' $gene_annotation >> $path/$gene/gene_anno.txt");	
# methylation data
    if (not -e "$path/$gene/me.txt") {
        system("head -1 $methylation_annotation > $path/$gene/cpg_anno.txt");
        system("cat data/me.header.txt > $path/$gene/me.txt");
        system("grep \"$gene\" $methylation_annotation | cut -f 1 > $path/$gene/me.probes.txt");
        open(PB, "$path/$gene/me.probes.txt");
        while(<PB>)
        {
            chomp;
            ($me_probe) = split;
            system("grep $me_probe $me_dataset >> $path/$gene/me.txt ");            
            system("grep $me_probe $methylation_annotation >> $path/$gene/cpg_anno.txt ");
            
        }
        close PB;
    }

    # do regression
    system("module load R; Rscript --vanilla script/run_lm.r --args $path/$gene/ex.$ex_probe.txt $path/$gene/me.txt");
    system("rm -rf $path/$gene");

}


sub run
{
    $cmd = shift;
    print "$cmd\n";
    system $cmd;
}


system("cat run/* > final.txt");
system("head -n 1 final.txt > ./MethylXcan.txt");
system("cat final.txt | sed '/CpG/d' >> ./MethylXcan.txt");
system("rm final.txt");




my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";

