#!/usr/bin/env perl
use strict;
use warnings;
use v5.10;

my $i=0;
my %store;
my @caller;
my @chr;
my %chrstore;
for my $i (0..$#ARGV/2){
	push @caller , $ARGV[$i+1+$#ARGV/2];
    open I, '<'.$ARGV[$i] or die $!;
    while(<I>){
    	chomp;
    	next if /^(#|\s+$)/;
    	my @l = split /\s+/;
    	unless (exists $chrstore{$l[0]}){
    		push @chr , $l[0];
    		$chrstore{$l[0]} = 1;
    	}
    	@{${$store{$l[0]}{$l[1]}}[$i]} = @l;
    }
    close I;
}

print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".join("\t",@caller)."\n";

for my $chr (@chr){
	for my $pos (sort {$a <=> $b} keys %{$store{$chr}}){
		my @gt;
		my @info;
		my @scores;
		my @refalt;
		my %tv;
		my %t;
		my $score='.';
		my @sortt;
		my $ref;
		my $alt;
        my @set;
		for my $i (0..$#caller){
            my $e = ${$store{$chr}{$pos}}[$i];
			unless ($e){
				$gt[$i] = './.';
			} else {
                push @set , $caller[$i];
                my @l = @{$e};
				for(split /;/ , $l[-3]){
					next if $_=~/indel/i || $_=~/type=/i;
					push @info , $caller[$i].':'.$_;
				}
				push @scores, $caller[$i].':SCORE='.$l[5];
				push @refalt, $caller[$i].':REF='.$l[3].';'.$caller[$i].':ALT='.$l[4];
				my @t = split /:/ , $l[-2];
				my @v = split /:/ , $l[-1];
				for (0..$#t){
					$tv{$i}{$t[$_]} = $v[$_];
					$t{$t[$_]} = 1;
				}
				if ($i==0){
					$score = $l[5] if $l[5]=~/\d/;
				}
				if ($#sortt==-1){
					@sortt = @t;
					$ref = $l[3];
					$alt = $l[4];
				}
			}
		}
		delete $t{$_} for @sortt;
		push @sortt , sort {$a cmp $b} keys %t;
		for my $i (0..$#caller){
			unless ($gt[$i]){
				my @v;
				push @v , exists $tv{$i}{$_} ? $tv{$i}{$_} : '.' for @sortt;
				$gt[$i] = join(':',@v);
			}
		}
		print $chr."\t".$pos."\t.\t".$ref."\t".$alt."\t".$score."\t".join(';',(@scores,@refalt,@info)).';set='.($#set==$#caller ? 'Intersection' : join('-',@set))."\t".join(':',@sortt)."\t".join("\t",@gt)."\n";
	}
}
