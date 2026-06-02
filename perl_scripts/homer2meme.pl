#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw($GunzipError);

our $VERSION = '0.2.0';

my $input       = '';
my $extract     = '';
my $pseudocount = 0.01;
my $background  = '0.25';
my $input_fmt   = 'homer';
my $matrix_type = 'auto';
my $do_rc       = 0;
my $trim_edges  = 0;
my $min_ic      = 0;
my $nsites;
my $evalue;
my $renormalize = 0;
my $alphabet    = 'ACGT';
my $show_version = 0;

GetOptions(
    'i=s' => \$input,
    'e=s' => \$extract,
    'a=f' => \$pseudocount,
    'b=s' => \$background,
    'f=s' => \$input_fmt,
    'format=s' => \$input_fmt,
    'input-format=s' => \$matrix_type,
    'alphabet=s' => \$alphabet,
    'rc'    => \$do_rc,
    'trim-edges=f' => \$trim_edges,
    'min-ic=f' => \$min_ic,
    'nsites=i' => \$nsites,
    'evalue=f' => \$evalue,
    'renormalize' => \$renormalize,
    'version' => \$show_version,
    'h'   => sub { usage() },
) or usage();

if ($show_version) {
    print "homer2meme $VERSION\n";
    exit 0;
}

usage() unless $input;
die "Error: -a must be > 0.\n" unless $pseudocount > 0;
my @background_values = parse_background($background);
die "Error: --nsites must be > 0.\n" if defined $nsites && $nsites <= 0;
die "Error: --evalue must be >= 0.\n" if defined $evalue && $evalue < 0;
die "Error: unknown format: $input_fmt\n" unless $input_fmt eq 'homer' || $input_fmt eq 'json';
die "Error: unknown input-format: $matrix_type\n" unless $matrix_type eq 'auto' || $matrix_type eq 'logodds' || $matrix_type eq 'probability';
die "Error: unknown alphabet: $alphabet\n" unless $alphabet =~ /^(ACGT|ACGU|PROTEIN)$/;

my %ALPHABETS = (
    'ACGT' => 'ACGT',
    'ACGU' => 'ACGU',
    'PROTEIN' => 'ACDEFGHIKLMNPQRSTVWY'
);
my %config = (
    extract => $extract,
    pseudocount => $pseudocount,
    background => \@background_values,
    input_fmt => $input_fmt,
    matrix_type => $matrix_type,
    do_rc => $do_rc,
    trim_edges => $trim_edges,
    min_ic => $min_ic,
    nsites => $nsites,
    evalue => $evalue,
    renormalize => $renormalize,
    alphabet => $alphabet,
    expected_cols => length($ALPHABETS{$alphabet} || $alphabet),
);

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

if ($input_fmt eq 'json') {
    parse_and_convert_json($fh, \%config);
} else {
    parse_and_convert_homer($fh, \%config);
}

if ($input ne '-') {
    if (ref $fh && $fh->isa('IO::Uncompress::Gunzip')) {
        $fh->close() or warn "Error closing gz file: $GunzipError";
    } else {
        close $fh or warn "Error closing file: $!";
    }
}

# ---------------------------------------------------------------------------

sub is_logodds {
    my ($row_ref, $matrix_type) = @_;
    return 1 if $matrix_type eq 'logodds';
    return 0 if $matrix_type eq 'probability';
    my $sum = 0;
    $sum += $_ for @$row_ref;
    my $is_probability = ($sum >= 0.98 && $sum <= 1.02);
    my $all_nonnegative = 1;
    for my $v (@$row_ref) {
        if ($v < 0) {
            $all_nonnegative = 0;
            last;
        }
    }
    if (!$is_probability && $all_nonnegative && $sum >= 0.90 && $sum <= 1.10) {
        warn sprintf(
            "Warning: row sum %.4f is close to the auto-detection boundary; use --input-format logodds or --input-format probability for reproducibility\n",
            $sum
        );
    }
    return !$is_probability;
}

sub logodds_to_prob {
    my ($row_ref, $pc, $bg_ref) = @_;
    my @bg = background_for_width($bg_ref, scalar(@$row_ref));
    my @raw;
    for my $i (0 .. $#{$row_ref}) {
        push @raw, 2 ** $row_ref->[$i] * $bg[$i];
    }
    my $total = $pc * scalar(@raw);
    $total += $_ for @raw;
    return map { ($_ + $pc) / $total } @raw;
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

sub calculate_ic {
    my ($matrix_ref, $alphabet) = @_;
    my $max_ic = ($alphabet && $alphabet eq 'PROTEIN') ? log(20) / log(2) : 2.0;
    my @ic_list;
    foreach my $row (@$matrix_ref) {
        my $h = 0;
        foreach my $p (@$row) {
            $h -= $p * (log($p) / log(2)) if $p > 0;
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

sub process_homer_motif {
    my ($id, $desc, $matrix_ref, $header_ref, $header_alphabet_ref, $config) = @_;

    my @mat = @$matrix_ref;

    if ($config->{do_rc}) {
        if ($config->{alphabet} eq 'ACGT' || $config->{alphabet} eq 'ACGU') {
            @mat = reverse_complement(\@mat, \$id);
        } else {
            warn "Warning: skipping motif '$id': reverse complement not supported for alphabet: $config->{alphabet}\n";
            return;
        }
    }

    if ($config->{trim_edges} > 0) {
        @mat = trim_edges(\@mat, $config->{trim_edges}, $config->{alphabet});
        if (!@mat) {
            warn "Warning: motif '$id' trimmed to empty matrix (IC threshold=$config->{trim_edges})\n";
            return;
        }
    }

    if ($config->{min_ic} > 0 && total_ic(\@mat, $config->{alphabet}) < $config->{min_ic}) {
        return;
    }

    if ($$header_ref) {
        if ($$header_alphabet_ref ne $config->{alphabet}) {
            warn "Warning: skipping motif '$id' with alphabet $config->{alphabet} (header uses $$header_alphabet_ref)\n";
            return;
        }
    } else {
        print_meme_header($config->{alphabet});
        $$header_ref = 1;
        $$header_alphabet_ref = $config->{alphabet};
    }
    print_meme_motif($id, $desc, \@mat, $config);
}

sub background_line {
    my ($alphabet) = @_;
    my $letters = $ALPHABETS{$alphabet} || $alphabet;
    return '' unless length($letters);
    my $count = length($letters);
    my $freq = 1 / $count;
    my $fmt = ($count == 4 || $count == 20) ? "%.2f" : "%.6f";
    my @parts;
    foreach my $ch (split //, $letters) {
        push @parts, sprintf("%s $fmt", $ch, $freq);
    }
    return join(" ", @parts);
}

sub print_meme_header {
    my ($alphabet) = @_;
    print "MEME version 4\n";
    print "\n";
    print "ALPHABET= $alphabet\n";
    print "\n";
    if ($alphabet eq 'ACGT' || $alphabet eq 'ACGU') {
        print "strands: + -\n";
        print "\n";
    }
    print "Background letter frequencies\n";
    print background_line($alphabet) . "\n";
    print "\n";
}

sub print_meme_motif {
    my ($id, $desc, $matrix_ref, $config) = @_;
    my $expected_cols = $config->{expected_cols};
    my $nsites = defined $config->{nsites} ? $config->{nsites} : 20;
    my $evalue = defined $config->{evalue} ? sprintf("%.6f", $config->{evalue}) : "0";
    my $width = scalar @$matrix_ref;
    print "MOTIF $id $desc\n";
    print "\n";
    print "letter-probability matrix: alength= $expected_cols w= $width nsites= $nsites E= $evalue\n";
    foreach my $row (@$matrix_ref) {
        my @values = $config->{renormalize} ? renormalized_row($row) : @$row;
        print "  " . join("  ", map { sprintf("%.6f", $_) } @values) . "\n";
    }
    print "\n";
}

sub parse_and_convert_homer {
    my ($fh, $config) = @_;
    my $header_printed = 0;
    my $header_alphabet = '';
    my $in_motif = 0;
    my $motif_id = '';
    my $description = '';
    my @matrix;
    my $extract = $config->{extract};
    my $pseudocount = $config->{pseudocount};
    my $background = $config->{background};
    my $matrix_type = $config->{matrix_type};

    while (<$fh>) {
        chomp;
        next unless length($_);

        if (/^>(.*)/) {
            my $rest = $1;

            if ($in_motif && @matrix) {
                process_homer_motif($motif_id, $description, \@matrix, \$header_printed, \$header_alphabet, $config);
            }
            @matrix = ();

            my @parts = split /\t/, $rest;
            my $mid  = defined $parts[0] ? $parts[0] : 'motif';
            my $desc = defined $parts[1] ? $parts[1] : $mid;

            if ($extract && $mid ne $extract && $desc ne $extract) {
                $in_motif = 0;
                next;
            }

            $motif_id    = $mid;
            $description = $desc;
            $in_motif    = 1;
            next;
        }

        next unless $in_motif;

        my @tokens = split /\s+/;
        my $all_numeric = 1;
        for my $t (@tokens) {
            unless ($t =~ /^-?[\d.]+([eE][+-]?\d+)?$/) {
                $all_numeric = 0;
                last;
            }
        }
        next unless $all_numeric && @tokens;

        my @row = map { $_ + 0 } @tokens;
        if (scalar(@row) != $config->{expected_cols}) {
            warn "Warning: skipping malformed matrix row (expected $config->{expected_cols} cols, got "
                 . scalar(@row) . "): $_\n";
            next;
        }

        if (is_logodds(\@row, $matrix_type)) {
            @row = logodds_to_prob(\@row, $pseudocount, $background);
        }
        push @matrix, \@row;
    }

    if ($in_motif && @matrix) {
        process_homer_motif($motif_id, $description, \@matrix, \$header_printed, \$header_alphabet, $config);
    }
}

sub parse_and_convert_json {
    my ($fh, $config) = @_;
    my $content = '';
    while (<$fh>) {
        $content .= $_;
    }

    require JSON::PP;
    my $data = eval { JSON::PP::decode_json($content) };
    if ($@) {
        die "Error: Invalid JSON: $@\n";
    }

    my $header_printed = 0;
    my $header_alphabet = '';
    my @motifs;
    my $extract = $config->{extract};
    my $pseudocount = $config->{pseudocount};
    my $background = $config->{background};
    my $matrix_type = $config->{matrix_type};

    if (ref $data eq 'HASH' && ref $data->{motifs} eq 'ARRAY') {
        for my $m (@{$data->{motifs}}) {
            my $id = $m->{id} || 'motif';
            my $desc = $m->{description} || $id;
            my $motif_alphabet = $m->{alphabet} || 'ACGT';
            my $expected_cols = length($ALPHABETS{$motif_alphabet} || $motif_alphabet);

            if ($extract && $id ne $extract && $desc ne $extract) {
                next;
            }

            my @matrix;
            if (ref $m->{matrix} eq 'ARRAY') {
                for my $row (@{$m->{matrix}}) {
                    if (ref $row eq 'ARRAY' && scalar(@$row) == $expected_cols) {
                        my @vals = map { $_ + 0 } @$row;
                        if (is_logodds(\@vals, $matrix_type)) {
                            @vals = logodds_to_prob(\@vals, $pseudocount, $background);
                        }
                        push @matrix, \@vals;
                    } else {
                        warn "Warning: skipping malformed matrix row (expected $expected_cols cols)\n";
                    }
                }
            }

            if (@matrix) {
                push @motifs, {
                    id => $id,
                    desc => $desc,
                    matrix => \@matrix,
                    alphabet => $motif_alphabet,
                    expected_cols => $expected_cols,
                };
            }
        }
    }

    if (@motifs) {
        for my $m (@motifs) {
            my %local_config = %$config;
            $local_config{alphabet} = $m->{alphabet};
            $local_config{expected_cols} = $m->{expected_cols};
            process_homer_motif(
                $m->{id},
                $m->{desc},
                $m->{matrix},
                \$header_printed,
                \$header_alphabet,
                \%local_config
            );
        }
    }
}

sub usage {
    print <<EOF;
Usage: $0 -i <input_file> [OPTIONS]

Convert HOMER motif format to MEME format.

Options:
    -i <file>   Input HOMER motif file (or '-' for stdin, supports .gz)
    -e <string> Extract only specified motif by id or description
    -a <float>  Pseudocount for log-odds to probability conversion (default: 0.01)
    -b <float[,float...]>  Background probability scalar or vector (default: 0.25)
    -f, --format <fmt>  Input format: homer (default) or json
    --input-format <fmt>  Matrix type: auto (default), logodds, or probability
    --alphabet <str> Alphabet: ACGT (default), ACGU, or PROTEIN
    --version        Show version and exit
    --rc                Output the reverse complement of the motif (DNA/RNA only)
    --trim-edges <float> Trim edges with information content below threshold
    --min-ic <float>    Filter out motifs with total information content below threshold
    --nsites <int>       Override MEME nsites metadata in output
    --evalue <float>     Override MEME E metadata in output
    --renormalize        Renormalize each row before writing MEME output
    -h          Show this help

Examples:
    $0 -i results/motifs.homer > raw/motifs.meme
    $0 -i results/motifs.homer.gz > raw/motifs.meme
    $0 -i results/motifs.homer -e "CTCF/Jaspar"
    $0 -i motifs.json -f json > motifs.meme
    $0 -i motifs.homer --input-format logodds
    $0 -i motifs.homer --input-format logodds -b 0.29,0.21,0.21,0.29
    $0 -i motifs.homer --rc
    cat motifs.homer | $0 -i -

EOF
    exit 0;
}
