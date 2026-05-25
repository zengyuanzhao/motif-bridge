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
my $do_rc       = 0;
my $trim_edges  = 0;
my $min_ic      = 0;

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
    'rc'    => \$do_rc,
    'trim-edges=f' => \$trim_edges,
    'min-ic=f' => \$min_ic,
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
my %config = (
    db => $db,
    motif_name => $motif_name,
    extract => $extract,
    bg => $bg,
    t_offset => $t_offset,
    output_fmt => $output_fmt,
    alphabet => $alphabet,
    do_rc => $do_rc,
    trim_edges => $trim_edges,
    min_ic => $min_ic,
    expected_cols => length($ALPHABETS{$alphabet} || $alphabet),
);

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
            process_motif($motif_id, $description, \@matrix, \%config, \@motifs);
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
            $config{expected_cols} = $1;
        }
        next;
    }

    if (/^URL/) {
        next;
    }
    if (/^\/\//) {
        if (@matrix) {
            process_motif($motif_id, $description, \@matrix, \%config, \@motifs);
        }
        $in_motif  = 0;
        $in_matrix = 0;
        @matrix    = ();
        next;
    }

    if ($in_matrix && /^\s*[\d\.-]/) {
        s/^\s+//;
        my @row = split /\s+/;
        if (scalar(@row) == $config{expected_cols}) {
            if (grep { $_ < 0 } @row) {
                warn "Warning: negative value in matrix row (expected probabilities): $_\n";
            }
            push @matrix, \@row;
        } else {
            warn "Warning: skipping malformed matrix row (expected $config{expected_cols} cols, got "
                 . scalar(@row) . "): $_\n";
        }
    }
}

if ($in_motif && @matrix) {
    process_motif($motif_id, $description, \@matrix, \%config, \@motifs);
}

if ($input ne '-') {
    if (ref $fh && $fh->isa('IO::Uncompress::Gunzip')) {
        $fh->close() or warn "Error closing gz file: $GunzipError";
    } else {
        close $fh or warn "Error closing file: $!";
    }
}

if ($config{output_fmt} eq 'json') {
    print_json(\@motifs);
}

# ---------------------------------------------------------------------------

sub process_motif {
    my ($id, $desc, $matrix_ref, $config, $motifs_ref) = @_;

    my @mat = @$matrix_ref;
    my $alphabet = $config->{alphabet};

    if ($config->{do_rc}) {
        if ($alphabet eq 'ACGT' || $alphabet eq 'ACGU') {
            @mat = reverse_complement(\@mat, \$id);
        } else {
            warn "Warning: reverse complement not supported for alphabet: $alphabet\n";
        }
    }

    if ($config->{trim_edges} > 0) {
        @mat = trim_edges(\@mat, $config->{trim_edges}, $alphabet);
        if (!@mat) {
            warn "Warning: motif '$id' trimmed to empty matrix (IC threshold=$config->{trim_edges})\n";
            return;
        }
    }

    if ($config->{min_ic} > 0 && total_ic(\@mat, $alphabet) < $config->{min_ic}) {
        return;
    }

    if ($config->{output_fmt} eq 'json') {
        push @$motifs_ref, {
            id => $id,
            description => $desc,
            matrix => [map { [@$_] } @mat],
            alphabet => $alphabet,
        };
    } else {
        my $score = calculate_score(\@mat, $config->{bg}, $config->{t_offset});
        print_motif($id, $desc, $score, \@mat);
    }
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

sub calculate_ic {
    my ($matrix_ref, $alphabet) = @_;
    my $max_ic = ($alphabet eq 'PROTEIN') ? log2(20) : 2.0;
    my @ic_list;
    foreach my $row (@$matrix_ref) {
        my $h = 0;
        foreach my $p (@$row) {
            $h -= $p * log2($p) if $p > 0;
        }
        my $ic = $max_ic - $h;
        $ic = 0 if $ic < 0;
        push @ic_list, $ic;
    }
    return @ic_list;
}

sub total_ic {
    my ($matrix_ref, $alphabet) = @_;
    my @ic = calculate_ic($matrix_ref, $alphabet);
    my $sum = 0;
    $sum += $_ for @ic;
    return $sum;
}

sub reverse_complement {
    my ($matrix_ref, $id_ref) = @_;
    my @rev = reverse @$matrix_ref;
    my @rc;
    for my $row (@rev) {
        push @rc, [$row->[3], $row->[2], $row->[1], $row->[0]];
    }
    $$id_ref .= '_RC';
    return @rc;
}

sub trim_edges {
    my ($matrix_ref, $threshold, $alphabet) = @_;
    my @ic = calculate_ic($matrix_ref, $alphabet);
    my $start = 0;
    while ($start < scalar(@ic) && $ic[$start] < $threshold) {
        $start++;
    }
    my $end = scalar(@ic);
    while ($end > $start && $ic[$end - 1] < $threshold) {
        $end--;
    }
    if ($start < $end) {
        return @$matrix_ref[$start .. $end - 1];
    } else {
        return ();
    }
}

sub print_motif {
    my ($id, $desc, $score, $matrix_ref) = @_;
    print ">$id\t$desc\t$score\t0\t0\t0\n";
    foreach my $row (@$matrix_ref) {
        print join("\t", map { sprintf('%.6f', $_) } @$row) . "\n";
    }
}

sub print_json {
    my ($motifs_ref) = @_;
    require JSON::PP;
    my $encoder = JSON::PP->new->allow_nonref;
    print "{\n";
    print "  \"version\": \"1.0\",\n";
    print "  \"source\": \"meme\",\n";
    print "  \"motifs\": [\n";
    for my $mi (0 .. $#{$motifs_ref}) {
        my $m = $motifs_ref->[$mi];
        print "    {\n";
        print "      \"id\": " . $encoder->encode($m->{id}) . ",\n";
        print "      \"description\": " . $encoder->encode($m->{description}) . ",\n";
        print "      \"alphabet\": " . $encoder->encode($m->{alphabet}) . ",\n" if $m->{alphabet};
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
    --rc                Output the reverse complement of the motif (DNA/RNA only)
    --trim-edges <float> Trim edges with information content below threshold
    --min-ic <float>    Filter out motifs with total information content below threshold
    -h           Show this help

Examples:
    $0 -i raw/motifs.meme -j JASPAR2026 > results/motifs.homer
    $0 -i raw/motifs.meme.gz -e MA0021.1
    $0 -i motifs.meme -b 0.25 -t 6
    $0 -i motifs.meme -f json > motifs.json
    $0 -i motifs.meme --rc
    cat motifs.meme | $0 -i -

EOF
    exit 0;
}
