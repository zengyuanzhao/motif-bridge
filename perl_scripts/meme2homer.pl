#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw($GunzipError);

our $VERSION = '0.3.1';

my $input       = '';
my $db          = 'NA';
my $motif_name  = '';
my $extract     = '';
my $bg          = '0.25';
my $t_offset    = 4;
my $output_fmt  = 'homer';
my $alphabet    = '';
my $do_rc       = 0;
my $trim_edges  = 0;
my $min_ic      = 0;
my $renormalize = 0;
my $keep_threshold = 0;
my $show_version = 0;
my $strict_mode = 0;

GetOptions(
    'i=s' => \$input,
    'j=s' => \$db,
    'k=s' => \$motif_name,
    'e=s' => \$extract,
    'b=s' => \$bg,
    't=f' => \$t_offset,
    'f=s' => \$output_fmt,
    'format=s' => \$output_fmt,
    'alphabet=s' => \$alphabet,
    'rc'    => \$do_rc,
    'trim-edges=f' => \$trim_edges,
    'min-ic=f' => \$min_ic,
    'renormalize' => \$renormalize,
    'keep-threshold' => \$keep_threshold,
    'version' => \$show_version,
    'strict' => \$strict_mode,
    'h'   => sub { usage() },
) or usage();

if ($show_version) {
    print "meme2homer $VERSION\n";
    exit 0;
}

usage() unless $input;
my @background_values = parse_background($bg);
my $alphabet_override = ($alphabet ne '');
if ($alphabet_override) {
    die "Error: unknown alphabet: $alphabet\n" unless $alphabet =~ /^(ACGT|ACGU|PROTEIN)$/;
}

my %ALPHABETS = (
    'ACGT' => 'ACGT',
    'ACGU' => 'ACGU',
    'PROTEIN' => 'ACDEFGHIKLMNPQRSTVWY'
);
my $effective_alphabet = $alphabet_override ? $alphabet : 'ACGT';
my %config = (
    db => $db,
    motif_name => $motif_name,
    extract => $extract,
    bg => \@background_values,
    t_offset => $t_offset,
    output_fmt => $output_fmt,
    alphabet => $effective_alphabet,
    do_rc => $do_rc,
    trim_edges => $trim_edges,
    min_ic => $min_ic,
    renormalize => $renormalize,
    keep_threshold => $keep_threshold,
    strict => $strict_mode,
    expected_cols => length($ALPHABETS{$effective_alphabet} || $effective_alphabet),
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
binmode($fh, ':encoding(UTF-8)');
binmode(STDOUT, ':encoding(UTF-8)');
binmode(STDERR, ':encoding(UTF-8)');

my $in_motif  = 0;
my $in_matrix = 0;
my $motif_id  = '';
my $description = '';
my $motif_nsites;
my $motif_evalue;
my @matrix;
my @motifs;

while (<$fh>) {
    chomp;

    if (/^ALPHABET=\s*(\S+)/) {
        if (!$alphabet_override) {
            my $detected = $1;
            $config{alphabet} = $detected;
            $config{expected_cols} = length($ALPHABETS{$detected} || $detected);
        }
        next;
    }

    if (/^MOTIF\s+(\S+)(?:\s+(.*))?/) {
        if ($in_motif && @matrix) {
            process_meme_motif($motif_id, $description, \@matrix, \%config, \@motifs, $motif_nsites, $motif_evalue);
        }

        $in_motif  = 1;
        $in_matrix = 0;
        $motif_id  = $1;
        my $original_name = defined $2 ? $2 : $1;
        $original_name =~ s/\s+/ /g;
        $original_name =~ s/^\s+|\s+$//g;

        $description = $motif_name ? "$motif_name/$db" : "$original_name/$db";
        $motif_nsites = undef;
        $motif_evalue = undef;

        if ($extract && $motif_id ne $extract && $original_name ne $extract) {
            $in_motif = 0;
            @matrix   = ();
            next;
        }

        @matrix = ();
        next;
    }

    next unless $in_motif;

    if (/^MOTIF\S/) {
        $in_matrix = 0;
        next;
    }

    if (/^letter-probability matrix:/) {
        $in_matrix = 1;
        if (/nsites=\s*(\d+)/) {
            $motif_nsites = $1 + 0;
        }
        if (/E=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)/) {
            $motif_evalue = $1 + 0;
        }
        if (/alength=\s*(\d+)/) {
            my $alength = $1;
            if ($alength != $config{expected_cols}) {
                my $message = "alength=$alength conflicts with alphabet $config{alphabet} "
                   . "(expected $config{expected_cols} cols); using alphabet-derived width";
                die "Error: $message\n" if $config{strict};
                warn "Warning: $message\n";
            } else {
                $config{expected_cols} = $alength;
            }
        }
        next;
    }

    if (/^URL/) {
        next;
    }
    if (/^\/\//) {
        if (@matrix) {
            process_meme_motif($motif_id, $description, \@matrix, \%config, \@motifs, $motif_nsites, $motif_evalue);
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
            if ($config{strict}) {
                validate_probability_row(\@row, "MEME matrix row for motif '$motif_id'");
            } elsif (grep { $_ < 0 } @row) {
                warn "Warning: negative value in matrix row (expected probabilities): $_\n";
            }
            push @matrix, \@row;
        } else {
            my $message = "skipping malformed matrix row (expected $config{expected_cols} cols, got "
                 . scalar(@row) . "): $_";
            die "Error: $message\n" if $config{strict};
            warn "Warning: $message\n";
        }
    }
}

if ($in_motif && @matrix) {
    process_meme_motif($motif_id, $description, \@matrix, \%config, \@motifs, $motif_nsites, $motif_evalue);
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

sub process_meme_motif {
    my ($id, $desc, $matrix_ref, $config, $motifs_ref, $nsites, $evalue) = @_;

    my @mat = @$matrix_ref;
    my $alphabet = $config->{alphabet};

    if ($config->{do_rc}) {
        if ($alphabet eq 'ACGT' || $alphabet eq 'ACGU') {
            @mat = reverse_complement(\@mat, \$id);
        } else {
            warn "Warning: skipping motif '$id': reverse complement not supported for alphabet: $alphabet\n";
            return;
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
            (defined $nsites ? (nsites => $nsites) : ()),
            (defined $evalue ? (evalue => $evalue) : ()),
        };
    } else {
        my $score = calculate_score(\@mat, $config->{bg}, $config->{t_offset}, $config->{renormalize});
        if ($score == 0) {
            warn "Warning: HOMER threshold for motif '$id' was clipped to 0 "
               . "(t_offset=$config->{t_offset}); scanning may be very permissive\n";
        }
        print_motif($id, $desc, $score, \@mat, $config->{renormalize});
    }
}

# ---------------------------------------------------------------------------

sub calculate_score {
    my ($matrix_ref, $bg_ref, $t_offset, $renormalize) = @_;
    my $score = 0;
    foreach my $row (@$matrix_ref) {
        my @values = $renormalize ? renormalized_row($row) : @$row;
        my @bg = background_for_width($bg_ref, scalar(@values));
        my $best_idx = 0;
        for my $i (0 .. $#values) {
            $best_idx = $i if $values[$i] > $values[$best_idx];
        }
        my $max_p = $values[$best_idx];
        $score += log2($max_p / $bg[$best_idx]) if $max_p > 0;
    }
    $score -= $t_offset;
    $score = 0 if $score < 0;
    return $score;
}

sub log2 {
    my ($x) = @_;
    return log($x) / log(2);
}

sub parse_background {
    my ($value) = @_;
    my @parts = split /,/, $value, -1;
    die "Error: -b must contain at least one value.\n" unless @parts;
    my @values;
    foreach my $part (@parts) {
        die "Error: invalid -b value: $value\n" unless $part =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
        my $v = $part + 0;
        die "Error: -b values must be in (0, 1].\n" unless $v > 0 && $v <= 1;
        push @values, $v;
    }
    if (scalar(@values) > 1) {
        my $sum = 0;
        $sum += $_ for @values;
        die "Error: -b vector must sum to 1.0.\n" if abs($sum - 1.0) > 1e-3;
    }
    return @values;
}

sub background_for_width {
    my ($bg_ref, $width) = @_;
    if (scalar(@$bg_ref) == 1) {
        return (($bg_ref->[0]) x $width);
    }
    die "Error: background length " . scalar(@$bg_ref) . " does not match row width $width\n"
        unless scalar(@$bg_ref) == $width;
    return @$bg_ref;
}

sub renormalized_row {
    my ($row_ref) = @_;
    my $sum = 0;
    $sum += $_ for @$row_ref;
    return @$row_ref unless $sum > 0;
    return map { $_ / $sum } @$row_ref;
}

sub validate_probability_row {
    my ($row_ref, $context) = @_;
    foreach my $v (@$row_ref) {
        die "Error: $context values must be in [0, 1]\n" if $v < 0 || $v > 1;
    }
    my $sum = 0;
    $sum += $_ for @$row_ref;
    die "Error: $context must sum to 1.0, got " . sprintf("%.6f", $sum) . "\n"
        if abs($sum - 1.0) > 1e-3;
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
    my ($id, $desc, $score, $matrix_ref, $renormalize) = @_;
    print ">$id\t$desc\t" . sprintf('%.6f', $score) . "\t0\t0\t0\n";
    foreach my $row (@$matrix_ref) {
        my @values = $renormalize ? renormalized_row($row) : @$row;
        print join("\t", map { sprintf('%.6f', $_) } @values) . "\n";
    }
}

sub print_json {
    my ($motifs_ref) = @_;
    require JSON::PP;
    my $encoder = JSON::PP->new->allow_nonref->ascii(0);
    print "{\n";
    print "  \"version\": \"1.0\",\n";
    print "  \"format\": \"motif-bridge-json\",\n";
    print "  \"motifs\": [\n";
    for my $mi (0 .. $#{$motifs_ref}) {
        my $m = $motifs_ref->[$mi];
        print "    {\n";
        print "      \"id\": " . $encoder->encode($m->{id}) . ",\n";
        print "      \"description\": " . $encoder->encode($m->{description}) . ",\n";
        print "      \"alphabet\": " . $encoder->encode($m->{alphabet}) . ",\n" if $m->{alphabet};
        print "      \"nsites\": $m->{nsites},\n" if defined $m->{nsites};
        print "      \"evalue\": " . sprintf("%.6f", $m->{evalue}) . ",\n" if defined $m->{evalue};
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
    -b <float[,float...]> Background probability scalar or vector (default: 0.25)
    -t <float>   Threshold offset in log2 bits (default: 4)
    -f, --format <fmt>  Output format: homer (default) or json
    --alphabet <str> Alphabet override: ACGT (DNA), ACGU (RNA), or PROTEIN
    --version        Show version and exit
    --rc                Output the reverse complement of the motif (DNA/RNA only)
    --trim-edges <float> Trim edges with information content below threshold
    --min-ic <float>    Filter out motifs with total information content below threshold
    --renormalize        Renormalize each row before writing HOMER output
    --keep-threshold     Keep an existing motif threshold when present; plain MEME input has none
    --strict             Fail on malformed or invalid probability matrix rows
    -h           Show this help

Examples:
    $0 -i raw/motifs.meme -j JASPAR2026 > results/motifs.homer
    $0 -i raw/motifs.meme.gz -e MA0021.1
    $0 -i motifs.meme -b 0.25 -t 6
    $0 -i motifs.meme -b 0.29,0.21,0.21,0.29
    $0 -i motifs.meme -f json > motifs.json
    $0 -i motifs.meme --rc
    cat motifs.meme | $0 -i -

EOF
    exit 0;
}
