

#$id = 'CRC1';

$id = $ARGV[0];
$pheno = $ARGV[1];
$phenod = $ARGV[2];
$diffd = $ARGV[3];
$leavep = $ARGV[4];
$outdir = $ARGV[5];
$dir = "$phenod/$pheno/$id/sp";
$diff = "$diffd/$pheno/$id/$id.abundance.wilcox_testing.tsv";



@group = ('Health',$pheno);
$n = 30;

open OUT,">$outdir/$id.keystone.species.list";
print "$id\n";


open LS,"Name_short";
while(<LS>){
    chomp;
    ($tax,$newtax) = split /\t/;
    $tax =~ s/[\r\n]//g;     
    $newtax =~ s/[\r\n]//g;   
    $newtax{"s__$tax"} = $newtax;
}


foreach $g(@group){
    $in = "$dir/cluster_$g/keystone_node.tsv";
    # $eigen_c = `cat $in |awk -F '\\t' '\$2==0 && \$4 == "True" {print \$6}'`;
    open IN,$in;
    while(<IN>){
        chomp;
        s/\r//;
        
        ($spe,$l,$pr,$eigen,$c,$p) = split /\t/;
        #print "$p\t$eigen_c\t$l\n";
        if($eigen eq 'True' && $l == 0){
            $eigen_c = $p
        }
    }


    chomp $eigen_c;
    $eigen_c =~ s/\r//;

    open IN,$in;
    while(<IN>){
        chomp;
        s/\r//;
        
        ($spe,$l,$pr,$eigen,$c,$p) = split /\t/;
        #print "$p\t$eigen_c\t$l\n";
        if($p eq $eigen_c && $l == 0){
            if($eigen eq 'True'){
                $in{$spe}{$g} = 2;
                print "$spe\t$g\n";
            }else{
                $in{$spe}{$g} = 1;
            }
        }
        if($eigen eq 'True' && $l == 0){print OUT "$newtax{$spe}\n";}
        $pr{$spe}{$g} = sprintf("%.4f",$pr);
        if($spe =~ /s__/){
            $total_pr{$spe} += $pr;
        }
    }
}
foreach $g(@group){
    foreach $spe(keys %in){
        if(!exists $in{$spe}{$g}){$in{$spe}{$g} = "0"}
    }
}

$count = 0;
foreach $spe(sort {$total_pr{$b}<=>$total_pr{$a}} keys %total_pr){
    
    $count ++;
    if($count >$n*2){next}
    $rank_in{$spe} = 1;
    #print "$spe\n";
        
}



open OU,">$outdir/$id.group_PR.tsv";
open OT,">$outdir/$id.diff.tsv";

open LS,$diff;
<LS>;
print OU "Taxa\tPR_control\tPR_case\teigen_control\teigen_case\n";
print OT "Taxa\tFeature\tRes\n";
while(<LS>){
    chomp;
    @l = split /\t/;
    $tax = $l[0];
   
    unless(exists $rank_in{$tax}|| exists $in{$tax}){next}
  
    $fdr = $l[-1];
    ($case,$control) = ($l[2],$l[3]);
    if($fdr < 0.05){
        if($case > $control){$sig = 'case'}else{$sig = 'control'}
    }else{
        $sig = 'Other';
    }
    # print OT "$newtax{$tax}\tAbundance\t$sig\n";
    if(exists $in{$tax}){
        $eigen = "$in{$tax}{$group[0]}\t$in{$tax}{$group[1]}";
        #print "$tax\t$pr{$tax}{$group[0]}\t$pr{$tax}{$group[1]}\t$fdr\t$control\t$case\n";
    }else{
        $eigen = "0\t0";
    }
    print OU "$newtax{$tax}\t$pr{$tax}{$group[0]}\t$pr{$tax}{$group[1]}\t$eigen\n";
}

open LS,"phylum.color";
while(<LS>){
    chomp;
    ($phy,$color) = split /\t/;
    $in_phy{$phy} = 1;
    #print "$phy\n";

}
open LS,"../data/cMD.select_2008.species_phylum.tsv";
while(<LS>){
    chomp;
    ($tax,$phy) = split /\t/;
    unless(exists $rank_in{$tax}|| exists $in{$tax}){next}

    #print "$phy\n";
    if(exists $in_phy{$phy}){
        $out_phy = $phy;
        $out_phy =~ s/p__//;
        print OT "$newtax{$tax}\tPhylum\t$out_phy\n";
    }else{
        print OT "$newtax{$tax}\tPhylum\tOther\n";
    }
}


open LS,"$leavep";
<LS>;
while(<LS>){
    chomp;
    s/-/_/g;
    s/\tNA/\tOther/;
    ($tax,$clu,$super) = split /\t/;
    unless(exists $rank_in{$tax}|| exists $in{$tax}){next}
    print OT "$newtax{$tax}\tCluster\t$clu\n";
    
}


