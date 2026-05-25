#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw($GunzipError);

my $input       = '';
my $extract     = '';
my $pseudocount = 0.01;
my $background  = 0.25;
my $input_fmt   = 'homer';
my $matrix_type = 'auto';
my $do_rc       = 0;
my $trim_edges  = 0;
my $min_ic      = 0;
my $alphabet    = 'ACGT';

GetOptions(
    'i=s' => \$input,
    'e=s' => \$extract,
    'a=f' => \$pseudocount,
    'b=f' => \$background,
    'f=s' => \$input_fmt,
    'format=s' => \$input_fmt,
    'input-format=s' => \$matrix_type,
    'alphabet=s' => \$alphabet,
    'rc'    => \$do_rc,
    'trim-edges=f' => \$trim_edges,
    'min-ic=f' => \$min_ic,
    'h'   => sub { usage() },
) or usage();

usage() unless $input;
die "Error: -a must be > 0.\n" unless $pseudocount > 0;
die "Error: -b must be in (0, 1].\n" unless $background > 0 && $background <= 1;
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
    background => $background,
    input_fmt => $input_fmt,
    matrix_type => $matrix_type,
    do_rc => $do_rc,
    trim_edges => $trim_edges,
    min_ic => $min_ic,
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
    return ($sum < 0.98 || $sum > 1.02);
}

sub logodds_to_prob {
    my ($row_ref, $pc, $bg) = @_;
    my @raw = map { 2 ** $_ * $bg } @$row_ref;
    my $total = $pc * scalar(@raw);
    $total += $_ for @raw;
    return map { ($_ + $pc) / $total } @raw;
}

sub calculate_ic {
    my ($matrix_ref) = @_;
    my $max_ic = 2.0;
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
    my ($matrix_ref) = @_;
    my @ic = calculate_ic($matrix_ref);
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
    my ($matrix_ref, $threshold) = @_;
    my @ic = calculate_ic($matrix_ref);
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

sub process_motif {
    my ($id, $desc, $matrix_ref, $header_ref, $header_alphabet_ref, $config) = @_;

    my @mat = @$matrix_ref;

    if ($config->{do_rc}) {
        @mat = reverse_complement(\@mat, \$id);
    }

    if ($config->{trim_edges} > 0) {
        @mat = trim_edges(\@mat, $config->{trim_edges});
        if (!@mat) {
            warn "Warning: motif '$id' trimmed to empty matrix (IC threshold=$config->{trim_edges})\n";
            return;
        }
    }

    if ($config->{min_ic} > 0 && total_ic(\@mat) < $config->{min_ic}) {
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
    print_meme_motif($id, $desc, \@mat, $config->{expected_cols});
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
    print "strands: + -\n";
    print "\n";
    print "Background letter frequencies\n";
    print background_line($alphabet) . "\n";
    print "\n";
}

sub print_meme_motif {
    my ($id, $desc, $matrix_ref, $expected_cols) = @_;
    my $width = scalar @$matrix_ref;
    print "MOTIF $id $desc\n";
    print "\n";
    print "letter-probability matrix: alength= $expected_cols w= $width nsites= 20 E= 0\n";
    foreach my $row (@$matrix_ref) {
        print "  " . join("  ", map { sprintf("%.6f", $_) } @$row) . "\n";
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
                process_motif($motif_id, $description, \@matrix, \$header_printed, \$header_alphabet, $config);
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
        process_motif($motif_id, $description, \@matrix, \$header_printed, \$header_alphabet, $config);
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
            process_motif(
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
    -b <float>  Background probability for log-odds conversion (default: 0.25)
    -f, --format <fmt>  Input format: homer (default) or json
    --input-format <fmt>  Matrix type: auto (default), logodds, or probability
    --alphabet <str> Alphabet: ACGT (default), ACGU, or PROTEIN
    --rc                Output the reverse complement of the motif (DNA/RNA only)
    --trim-edges <float> Trim edges with information content below threshold
    --min-ic <float>    Filter out motifs with total information content below threshold
    -h          Show this help

Examples:
    $0 -i results/motifs.homer > raw/motifs.meme
    $0 -i results/motifs.homer.gz > raw/motifs.meme
    $0 -i results/motifs.homer -e "CTCF/Jaspar"
    $0 -i motifs.json -f json > motifs.meme
    $0 -i motifs.homer --input-format logodds
    $0 -i motifs.homer --rc
    cat motifs.homer | $0 -i -

EOF
    exit 0;
}
