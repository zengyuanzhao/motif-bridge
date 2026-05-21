#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw($GunzipError);

my $input       = '';
my $db          = 'NA';
my $motif_name  = '';
my $extract     = '';
my $bg          = 0.25;
my $t_offset    = 4;
my $output_fmt  = 'homer';
my $alphabet    = 'ACGT';

GetOptions(
    'i=s' => \$input,
    'j=s' => \$db,
    'k=s' => \$motif_name,
    'e=s' => \$extract,
    'b=f' => \$bg,
    't=f' => \$t_offset,
    'f=s' => \$output_fmt,
    'format=s' => \$output_fmt,
    'alphabet=s' => \$alphabet,
    'h'   => sub { usage() },
) or usage();

usage() unless $input;
die "Error: -b must be in (0, 1].\n" unless $bg > 0 && $bg <= 1;
die "Error: unknown alphabet: $alphabet\n" unless $alphabet =~ /^(ACGT|ACGU|PROTEIN)$/;

my %ALPHABETS = (
    'ACGT' => 'ACGT',
    'ACGU' => 'ACGU',
    'PROTEIN' => 'ACDEFGHIKLMNPQRSTVWY'
);
my $expected_cols = length($ALPHABETS{$alphabet} || $alphabet);

die "Error: unknown format: $output_fmt\n" unless $output_fmt eq 'homer' || $output_fmt eq 'json';

my $fh;
if ($input eq '-') {
    $fh = \*STDIN;
} elsif ($input =~ /\.gz$/) {
    $fh = IO::Uncompress::Gunzip->new($input)
        or die "Cannot open $input: $GunzipError";
} else {
    open $fh, '<', $input or die "Cannot open $input: $!";
}

my $in_motif  = 0;
my $in_matrix = 0;
my $motif_id  = '';
my $description = '';
my @matrix;
my @motifs;

while (<$fh>) {
    chomp;

    if (/^MOTIF\s+(\S+)(?:\s+(.*))?/) {
        if ($in_motif && @matrix) {
            if ($output_fmt eq 'json') {
                push @motifs, {
                    id => $motif_id,
                    description => $description,
                    matrix => [map { [@$_] } @matrix],
                    alphabet => $alphabet,
                };
            } else {
                my $score = calculate_score(\@matrix, $bg, $t_offset);
                print_motif($motif_id, $description, $score, \@matrix);
            }
        }

        $in_motif  = 1;
        $in_matrix = 0;
        $motif_id  = $1;
        my $original_name = defined $2 ? $2 : $1;
        $original_name =~ s/\s+/ /g;
        $original_name =~ s/^\s+|\s+$//g;

        $description = $motif_name ? "$motif_name/$db" : "$original_name/$db";

        if ($extract && $motif_id ne $extract && $original_name ne $extract) {
            $in_motif = 0;
            @matrix   = ();
            next;
        }

        @matrix = ();
        next;
    }

    next unless $in_motif;

    if (/^letter-probability matrix:/) {
        $in_matrix = 1;
        if (/alength=\s*(\d+)/) {
            $expected_cols = $1;
        }
        next;
    }

    if (/^URL/) {
        next;
    }
    if (/^\/\//) {
        if (@matrix) {
            if ($output_fmt eq 'json') {
                push @motifs, {
                    id => $motif_id,
                    description => $description,
                    matrix => [map { [@$_] } @matrix],
                    alphabet => $alphabet,
                };
            } else {
                my $score = calculate_score(\@matrix, $bg, $t_offset);
                print_motif($motif_id, $description, $score, \@matrix);
            }
        }
        $in_motif  = 0;
        $in_matrix = 0;
        @matrix    = ();
        next;
    }

    if ($in_matrix && /^\s*[\d.]/) {
        s/^\s+//;
        my @row = split /\s+/;
        if (scalar(@row) == $expected_cols) {
            push @matrix, \@row;
        } else {
            warn "Warning: skipping malformed matrix row (expected $expected_cols cols, got "
                 . scalar(@row) . "): $_\n";
        }
    }
}

if ($in_motif && @matrix) {
    if ($output_fmt eq 'json') {
        push @motifs, {
            id => $motif_id,
            description => $description,
            matrix => [map { [@$_] } @matrix],
                    alphabet => $alphabet,
        };
    } else {
        my $score = calculate_score(\@matrix, $bg, $t_offset);
        print_motif($motif_id, $description, $score, \@matrix);
    }
}

if ($input ne '-') {
    if (ref $fh && $fh->isa('IO::Uncompress::Gunzip')) {
        $fh->close() or warn "Error closing gz file: $GunzipError";
    } else {
        close $fh or warn "Error closing file: $!";
    }
}

if ($output_fmt eq 'json') {
    print_json(\@motifs);
}

# ---------------------------------------------------------------------------

sub calculate_score {
    my ($matrix_ref, $bg, $t_offset) = @_;
    my $score = 0;
    foreach my $row (@$matrix_ref) {
        my $max_p = 0;
        foreach my $p (@$row) {
            $max_p = $p if $p > $max_p;
        }
        $score += log2($max_p / $bg) if $max_p > 0;
    }
    $score -= $t_offset;
    $score = 0 if $score < 0;
    return sprintf('%.6f', $score);
}

sub log2 {
    my ($x) = @_;
    return log($x) / log(2);
}

sub print_motif {
    my ($id, $desc, $score, $matrix_ref) = @_;
    print ">$id\t$desc\t$score\t0\t0\t0\n";
    foreach my $row (@$matrix_ref) {
        print join("\t", map { sprintf('%.6f', $_) } @$row) . "\n";
    }
}

sub escape_json {
    my ($s) = @_;
    $s =~ s/\\/\\\\/g;
    $s =~ s/"/\\"/g;
    $s =~ s/\n/\\n/g;
    $s =~ s/\r/\\r/g;
    $s =~ s/\t/\\t/g;
    return $s;
}

sub print_json {
    my ($motifs_ref) = @_;
    print "{\n";
    print "  \"version\": \"1.0\",\n";
    print "  \"source\": \"meme\",\n";
    print "  \"motifs\": [\n";
    for my $mi (0 .. $#{$motifs_ref}) {
        my $m = $motifs_ref->[$mi];
        print "    {\n";
        print "      \"id\": \"" . escape_json($m->{id}) . "\",\n";
        print "      \"description\": \"" . escape_json($m->{description}) . "\",\n";
        print "      \"alphabet\": \"" . escape_json($m->{alphabet}) . "\",\n" if $m->{alphabet};
        print "      \"matrix\": [\n";
        for my $ri (0 .. $#{$m->{matrix}}) {
            my $row = $m->{matrix}->[$ri];
            my $vals = join(", ", map { sprintf("%.6f", $_) } @$row);
            if ($ri < $#{$m->{matrix}}) {
                print "        [$vals],\n";
            } else {
                print "        [$vals]\n";
            }
        }
        print "      ]\n";
        if ($mi < $#{$motifs_ref}) {
            print "    },\n";
        } else {
            print "    }\n";
        }
    }
    print "  ]\n";
    print "}\n";
}

sub usage {
    print <<EOF;
Usage: $0 -i <input_file> [OPTIONS]

Convert MEME format to HOMER motif format.

Options:
    -i <file>    Input MEME format file (or '-' for stdin, supports .gz)
    -j <string>  Database name (default: NA)
    -k <string>  Motif name to use (default: name from MEME file)
    -e <string>  Extract only specified motif by id or name
    -b <float>   Background probability (default: 0.25, uniform)
    -t <float>   Threshold offset in log2 bits (default: 4)
    -f, --format <fmt>  Output format: homer (default) or json
    --alphabet <str> Alphabet: ACGT (DNA, default), ACGU (RNA), or PROTEIN
    -h           Show this help

Examples:
    $0 -i raw/motifs.meme -j JASPAR2026 > results/motifs.homer
    $0 -i raw/motifs.meme.gz -e MA0021.1
    $0 -i motifs.meme -b 0.25 -t 6
    $0 -i motifs.meme -f json > motifs.json
    cat motifs.meme | $0 -i -

EOF
    exit 0;
}
