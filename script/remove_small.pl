####download perl script from http://www.bioinformatics-made-simple.com/2012/07/how-to-filter-sequence-by-their-length.html######
#!/usr/bin/perl 
use strict;
use warnings;
 
my $minlen = shift or die "Error: `minlen` parameter not provided\n";
{
    local $/=">";
    while(<>) {
        chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;
        print ">$_" if($seqlen >= $minlen);
    }
    local $/="\n";
}
