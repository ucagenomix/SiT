use strict;
use warnings;
use POSIX;

my $sample = "hip142079";
my $nb = 2500;

#my $sample = "hip142078";
#my $nb = 2560;

#my $sample = "mob164315";
#my $nb = 918;

open(IN, $sample.".snp100838.RN3.QV6_snpmatrix.txt") || die $_;
my $header=<IN>;
$header=~ s/\t/,/g;

my %ai=();
while(my $ligne=<IN>){
	chomp($ligne);
	$ligne=~ s/\t/,/g;
	my @tab = split(",",$ligne);
	
	if($tab[1] =~ /G$/){ $ai{$tab[0].",".$tab[1]} = $ligne; }
	elsif($tab[1] =~ /A$/){ $ai{$tab[0].",".$tab[1]} = $ligne; }
	#else{
	#  print $tab[0].",".$tab[1]."\n";
	#}
}
close(IN);

print values(%ai)."\n";

my %already_in=();
my $tot=0;
open(OUT, "> ".$sample.".snp100838.RN3.QV6_snpmatrix.csv") || die $_;
print OUT $header;
foreach my $id (keys %ai){
  
  if($id =~ /G$/){
    $tot++;
    my $key_A = $id;
    $key_A =~ s/G$/A/;
    
    if($ai{$key_A}){
      print OUT $ai{$id}."\n";
      print OUT $ai{$key_A}."\n";
    }
    else{
      print OUT $ai{$id}."\n";
      print OUT $key_A.",na";
      for(my $i=0;$i<$nb;$i++){ print OUT ",0"; }
      print OUT "\n";
    }
    
    $already_in{$id} = 1;
  }
  elsif($id =~ /A$/){
    
    my $key_G = $id;
    $key_G =~ s/A$/G/;
    
    if(!$ai{$key_G}){
      $tot++;
      $tot++;
      print OUT $ai{$id}."\n";
      print OUT $key_G.",na";
      for(my $i=0;$i<$nb;$i++){ print OUT ",0"; }
      print OUT "\n";
      
    }
    
  }
    
}

close(OUT);

print $tot."\n";
