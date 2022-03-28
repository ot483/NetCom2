#!/usr/local/bin/perl -w


#----------- files ----------

my $file ="NTC_G210-ReadName_to_KeggName_Ultimate.txt";

open IN ,$file; 
open (OUT,">protein_EC_codes_total_NTC_G210.txt");



while (<IN>){
	chomp $_;
	if ($_=~/^(\S+).+\[EC:(.+)\]/){
		my @ecs=split(/\s/,$2);
		foreach my $ec(@ecs){
			print OUT "$1\t$ec\n";
			}
	
	}
	

}
close IN;


