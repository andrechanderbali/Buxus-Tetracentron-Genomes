#!/usr/bin/perl -w

use strict;

print "START - "; system "date";

my $count = $ARGV[0];
my $infile = $ARGV[1];
my $outfile = $ARGV[2];

# read the blast reports from ath and get the best hit
my %hit = ();
print "Reading blast reports\n";
open IN, $infile;
while (<IN>) {
	chomp;
	my $line = $_;
	my @field = split(/\t/,$line);
	my ($query, $hit, $evalue) = ($field[0], $field[1], $field[10]);
	#print "$query\t$hit\t$evalue\n";
	$hit{$query}{$evalue}{$hit} = $line;
}
close IN;

# scroll through and get the best hit(s) for each query and print to summary file 
print "Creating summary report\n"; system "date";
open OUT, ">$outfile";
foreach my $id (sort keys %hit){
	my $num_hits = 0;
	foreach my $evalue (sort {$a<=>$b} keys %{$hit{$id}}) {
		foreach my $hit (sort keys %{$hit{$id}{$evalue}}) {
			#print OUT "$id\t$hit\t$evalue\n" if $num_hits < $count;	# print out only query, hit, and evalue
			print OUT "$hit{$id}{$evalue}{$hit}\n" if $num_hits < $count;	# print out entire line
			$num_hits++;
		}
	}
}
close OUT;

print "STOP - "; system "date";
exit;
