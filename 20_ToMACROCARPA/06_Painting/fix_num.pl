#!/usr/bin/perl
use strict;
use warnings;

my $header = 1;

while (<>) {
    chomp;
    my @cols = split /\t/;

    # Print header without processing
    if ($header) {
        print join("\t", @cols), "\n";
        $header = 0;
        next;
    }

    # Skip row if column 1 is not an integer divisible by 10
    if (!defined $cols[1] || $cols[1] !~ /^\d+$/ || $cols[1] % 2 != 0) {
        next;
    }

    next unless defined $cols[2] && defined $cols[3];
    next unless $cols[2] =~ /^[01]$/ && $cols[3] =~ /^[01]$/;

    # Process only if column 3 is defined
    if (defined $cols[2]) {
        if ($cols[2] eq '1') {
            for my $i (4..529) {
                if (defined $cols[$i]) {
                    if ($cols[$i] eq '1') {
                        $cols[$i] = '0';
                    } elsif ($cols[$i] eq '0') {
                        $cols[$i] = '1';
                    } else {
                        $cols[$i] = 'NA';
                    }
                }
            }
        } elsif ($cols[2] eq '0') {
            for my $i (4..529) {
                if (defined $cols[$i]) {
                    if ($cols[$i] ne '1' && $cols[$i] ne '0') {
                        $cols[$i] = 'NA';
                    }
                }
            }
        }
        # else: column 3 not 0 or 1, row remains unchanged
    }

    print join("\t", @cols), "\n";
}