#!/usr/bin/perl -w
# Author : Charles VAN GOETHEM

use strict;
use warnings;

use HTTP::Tiny;

use JSON;

use Cwd;

use Getopt::Long;
use Pod::Usage;

use Data::Dumper;

my $makeSourcesFiles_VER = '0.0.1b';

################################################################################
################################################################################

my $pwd = getcwd();

# config server
my $http = HTTP::Tiny->new();
my $server = 'https://rest.ensembl.org';

# Mandatory arguments
my ($bed_file,$comment);

# Optional arguments
my ($output,$example);

# Supports arguments
my ($opt_help, $opt_man, $opt_versions);

GetOptions(
    'e|example'	    => \$example,
	'b|bed=s'		=> \$bed_file,
    'c|comment=s'	=> \$comment,
    'o|output=s'	=> \$output,
    'help|?'		=> \$opt_help,
    'version!'		=> \$opt_versions,
    'man!'			=> \$opt_man
) or pod2usage(-verbose => 1) && exit;

pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 2) if defined $opt_man;

# check args

# check if launch example :
if($example) {
    $bed_file = "test/test.bed";
	$comment = "test_V1";
	$output="test/";
}

# Mandatory
if (not defined $bed_file) {
	pod2usage(-verbose => 1, -message => "Error : parameters \"-b|--bed\" is undefined.");
	exit;
} elsif ($bed_file !~ /(.+\.bed)$/o) {
	pod2usage(-verbose => 1, -message => "Error : parameters \"-b|--bed\" not correctly format (expected *.bed)");
	exit;
} else {
    $bed_file = $1;
}

# Optional output
if ($output ne "") {
	pod2usage(-verbose => 1, -message => "Error : '$output' does not exist !\n") unless (-e $output);
	pod2usage(-verbose => 1, -message => "Error : '$output' is not a repertory !\n") unless (-d $output);
} else {
	$output = $pwd;
	print STDERR "Warning : All the outputs will be create here : '$output'";
}

################################################################################
################################################################################

open(OUT , ">".$output."/".$comment.".csv") or die "Could not open file '".$output."/".$comment."' $!";

open(BED, $bed_file) or die "Cannot open file $bed_file $!";

my $bed_array;
my $j = 0;
while(<BED>)
{
	print STDERR ++$j."\n";
    chomp;
	# CHR	EX	GENE	ORIENTATION	PCR_NAME	PCR_SIZE	PCR_START	PCR_STOP	SEQ_SIZE	SEQ_START	SEQ_STOP	AMORCE1	INSERT	AMORCE2
    my ($chr, $start, $end, $name, $size5, $size3) = split;
	my ($gene, $exon);
	if ($name =~ m/([A-Za-z0-9]+)Ex([0-9]+)[a-z]*/g) {
		$gene = $1;
		$exon =$2;
	} else {
		$gene = $name;
		$exon = 0;
	}

    # get sequence of amorce 1 (chromosomal orientation)
    my $amorce1 = "/sequence/region/human/$chr:$start..".($start + $size5 - 1).":1?coord_system_version=GRCh37";
    my $response1 = $http->get($server.$amorce1, {
        headers => {
            'Content-type' => 'text/plain'
        }
    });

    # get sequence of amorce 2 (inverse of chromosomal orientation)
    my $amorce2 = "/sequence/region/human/$chr:".($end - $size3 + 1)."..".$end.":1?coord_system_version=GRCh37";
    my $response2 = $http->get($server.$amorce2, {
        headers => {
            'Content-type' => 'text/plain'
        }
    });

    # get sequence of insert (chromosomal orientation)
    my $insert = "/sequence/region/human/$chr:".($start + $size5)."..".($end - $size3).":1?coord_system_version=GRCh37";
    my $response3 = $http->get($server.$insert, {
        headers => {
            'Content-type' => 'text/plain'
        }
    });

    # get the strand of gene
    my $strand = "/overlap/region/human/$chr:".$start."..".$end.":1?feature=gene;coord_system_version=GRCh37";
    # print STDERR $strand;
    my $response4 = $http->get($server.$strand, {
        headers => {
            'Content-type' => 'application/json'
        }
    });

	my $hash_for_gene  = decode_json $response4->{content};
    # check if gene is known else force strand +
	my $found_gene = 0;
	foreach my $ensembl_gene (@{$hash_for_gene}) {
		if($ensembl_gene->{external_name} eq $gene) {
			if($ensembl_gene->{strand} > 0) {
				$strand = '+';
			} else {
				$strand = '-';
			}
			$found_gene=1;
			last;
		}
	}
	if(!$found_gene) {
		$strand= '+';
	}

    # Create data for this amplicon
	my $amplicon = {
		'chr'   => $chr,
		'exon'     => int($exon),
		'gene'     => $gene,
		'orientation' => $strand,
		'PCR_name'  => $name,
		'PCR_size' => ($end - $start + 1),
		'PCR_start' => $start,
		'PCR_stop'   => $end,
		'SEQ_size'   => (($end - $size3) - ($start + $size5) + 1),
		'SEQ_start'   => ($start + $size5),
		'SEQ_stop'   => ($end - $size3),
		'amorce1'  => $response1->{content},
		'insert'  => $response3->{content},
		'amorce2'  => $response2->{content},
		'size5' => $size5,
		'size3' => $size3
	};
    push(@{$bed_array}, $amplicon);


    # force script sleeping each 3 iteration because only 15 request per seconds
    # are authorized by ensembl rest api
    # https://github.com/Ensembl/ensembl-rest/wiki/Rate-Limits
    if(($j%3) == 0) {
        sleep(1);
    }
}

# Split amplicon by group in case of overlapping
# https://stackoverflow.com/questions/10395383/sorting-an-array-of-hash-by-multiple-keys-perl
@{$bed_array} = sort {
	$a->{chr} cmp $b->{chr} or
    $a->{PCR_start} <=> $b->{PCR_start}
} @{$bed_array};

my $i = "a";

my %order = (
	$i => []
);

my $valid =1;
$j=0;


# print results
print OUT "COMMENT;AMPLICON;CHR;EX;GENE;ORIENTATION;PCR_NAME;PCR_SIZE;PCR_START;PCR_STOP;SEQ_SIZE;SEQ_START;SEQ_STOP;AMORCE1;INSERT;AMORCE2\n";

foreach my $amp (@{$bed_array}) {
	$j++;
	foreach my $lvl (sort keys(%order)) {
		$valid = 1;
		foreach my $amplicon (@{$order{$lvl}}) {
			if ( %{$amp}{"PCR_start"} >= %{$amplicon}{"PCR_start"} && %{$amp}{"PCR_start"} <= %{$amplicon}{"PCR_stop"} && %{$amp}{"chr"} eq %{$amplicon}{"chr"}) {
				$valid = 0;
			}
		}
		if($valid) {
			my $s =  @{$order{$lvl}};
			@{$order{$lvl}}[$s] = $amp;
			print OUT $comment."_".$j.";";
			print OUT "$lvl;";
			print OUT $amp->{chr}.";";
			print OUT $amp->{exon}.";";
			print OUT $amp->{gene}.";";
			print OUT $amp->{orientation}.";";
			print OUT $amp->{PCR_name}.";";
			print OUT $amp->{PCR_size}.";";
			print OUT $amp->{PCR_start}.";";
			print OUT $amp->{PCR_stop}.";";
			print OUT $amp->{SEQ_size}.";";
			print OUT $amp->{SEQ_start}.";";
			print OUT $amp->{SEQ_stop}.";";
			print OUT $amp->{amorce1}.";";
			print OUT $amp->{insert}.";";
			print OUT $amp->{amorce2}."\n";
			last;
		}
	}
	if(!$valid) {
		$i++;
		$order{$i}[0] = $amp;
		print OUT $comment."_".$j.";";
		print OUT $i.";";
		print OUT $amp->{chr}.";";
		print OUT $amp->{exon}.";";
		print OUT $amp->{gene}.";";
		print OUT $amp->{orientation}.";";
		print OUT $amp->{PCR_name}.";";
		print OUT $amp->{PCR_size}.";";
		print OUT $amp->{PCR_start}.";";
		print OUT $amp->{PCR_stop}.";";
		print OUT $amp->{SEQ_size}.";";
		print OUT $amp->{SEQ_start}.";";
		print OUT $amp->{SEQ_stop}.";";
		print OUT $amp->{amorce1}.";";
		print OUT $amp->{insert}.";";
		print OUT $amp->{amorce2}."\n";
	}
}

# print Dumper \%order;
close(OUT);

# print results for overlapping
foreach my $key (sort keys(%order)) {
	open(OUT, ">".$output."/".$comment."-igv_".$key.".bed") or die "Could not open file '".$output."/".$comment."' $!";
	foreach my $amp (@{$order{$key}}) {
		print OUT $amp->{chr}."\t";
		print OUT $amp->{PCR_start}."\t";
		print OUT $amp->{PCR_stop}."\t";
		print OUT $amp->{PCR_name}."\t";
		print OUT $amp->{exon}."\t";
		print OUT $amp->{orientation}."\t";
		print OUT $amp->{gene}."\n";
	}
	close(OUT);
}

##########################################################################################
##########################################################################################

BEGIN {
	delete @ENV{'IFS', 'CDPATH', 'ENV', 'BASH_ENV', 'PATH'};
}

END{
  if(defined $opt_versions){
    print
      "\nModules, Perl, OS, Program info:\n",
	  "    strict                 $strict::VERSION\n",
	  "    warnings               $warnings::VERSION\n",
	  "    Getopt::Long           $Getopt::Long::VERSION\n",
	  "    Pod::Usage             $Pod::Usage::VERSION\n",
	  "    HTTP::Tiny             $HTTP::Tiny::VERSION\n",
	  "    JSON                   $JSON::VERSION\n",
	  "    Cwd                    $Cwd::VERSION\n",
      "    Perl                   $]\n",
      "    OS                     $^O\n",
	  "    makeSourcesFiles.pl    $makeSourcesFiles_VER\n";
  }
}

##########################################################################################
##########################################################################################

__END__

=pod

=encoding UTF-8

=head1 NAME

makeSourcesFiles.pl

=head1 VERSION

version 0.0.1a

=head1 SYNOPSIS

makeSourcesFiles.pl -b input.bed -c comment [-o output.vcf]

Get a bed of amplicons and transform it on sources files for SARDINe

=head1 DESCRIPTION

This script transform a bed of sources files for SARDINe

=head1 ARGUMENTS

=head2 General

    -h,--help       Print this help

    -m,--man        Open man page

=head2 Mandatory arguments

    -b,--bed        input.bed     Specify the bed on input

    -c,--comment    comment       String for the comment column (define version)

=head2 Optional arguments

    -o,--output     /path/output/repertory    You can specify an output repertory (default: ./)

    -e,--example    launch example (equivalent to: 'perl addAF.pl -v test/test.bed -c test/')

=head1 AUTHORS

=over 4

=item -

Charles VAN GOETHEM

=back

=head1 TODO

=over 4

=item -

Generate bed files !

=item -

Add option to change genome assembly

=item -

Improved speed (algorithm really not optimized)

=item -

Clean and comment code

=back

=head1 LICENCE

MIT License

Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=cut
