#! /usr/bin/env perl
use strict;
use warnings;
use v5.10;
use List::MoreUtils qw(indexes);
use List::Util qw(sum);

my @mem;
my $chr='';
my $pos=0;
while(<>){
    my $l = $_;
    chomp $l;
    my @l = split /\s+/,$l;
    if ($l[0]=~/^#/){
        say $l;
        next;
    }
    next if $l[-1] =~ /^\.\//;
    my ($dp) = grep { $_ =~ /^DP$/i } split /:/,$l[-2];
    unless ($dp){ 
        my @t = split /:/,$l[-2];
        my @v = split /:/,$l[-1];
        if ($l[-3]=~/[\s;]TC=([^\s;]+)/){ #platypus fix
            $dp = $1;
            splice @t, 1, 0, 'DP';
            splice @v, 1, 0, $dp;
        } elsif($l[-3]=~/[\s;]DP=([^\s;]+)/){ #haplotypecaller fix
            $dp = $1;
            splice @t, 1, 0, 'DP';
            splice @v, 1, 0, $dp;
        } elsif($l[-3]=~/[\s;]TLOD=([^\s;]+)/){ #mutect2 fix
            my $gq = $1;
            my ($i) = indexes { $_ =~ /^AD$/i } @t;
            $dp = sum split /,/,$v[$i];
            splice @t, $i, 0, 'DP';
            splice @v, $i, 0, $dp;
            splice @t, $i, 0, 'GQ';
            splice @v, $i, 0, $gq;
        } else {
            say "foo";
            next;
        }
        $l[-1] = join ':' , @v;
        $l[-2] = join ':' , @t;
        $l = join "\t" , @l;
    }
    if ($#mem == -1){
    	push @mem,$l;
	} else {
        my @m = split /\s+/,$mem[-1];
        unless ($l[0] eq $m[0] && $l[1] == $m[1]){
            &select;
            @mem=();
        }
        push @mem,$l;
    }
}
&select;

sub select (){
    my %gq;
    for my $m (@mem){
        my @m = split /\s+/ , $m;
        next unless $m[-1]=~/^(0|1)\/(0|1)/;
        my ($i) = indexes { $_ =~ /^GQ$/i } split /:/,$m[-2];
        my ($j) = indexes { $_ =~ /^DP$/i } split /:/,$m[-2];
        my @v = split /:/,$m[-1];
        $gq{$v[$i]}{$v[$j]} = $m;
    }
    for my $i (reverse sort {$a <=> $b} keys %gq){
        for my $j (reverse sort {$a <=> $b} keys %{$gq{$i}}){
            say $gq{$i}{$j};
            last;
        }
        last;
    }
}
