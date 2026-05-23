#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw($GunzipError);

my $input       = '';
my $extract     = '';
my $pseudocount = 0.01;
my $input_fmt   = 'homer';
my $matrix_type = 'auto';
my $do_rc       = 0;
my $trim_edges  = 0;
my $min_ic      = 0;

GetOptions(
    'i=s' => \$input,
    'e=s' => \$extract,
    'a=f' => \$pseudocount,
    'f=s' => \$input_fmt,
    'format=s' => \$input_fmt,
    'input-format=s' => \$matrix_type,
    'rc'    => \$do_rc,
    'trim-edges=f' => \$trim_edges,
    'min-ic=f' => \$min_ic,
    'h'   => sub { usage() },
) or usage();

usage() unless $input;
die "Error: -a must be > 0.\n" unless $pseudocount > 0;
die "Error: unknown format: $input_fmt\n" unless $input_fmt eq 'homer' || $input_fmt eq 'json';
die "Error: unknown input-format: $matrix_type\n" unless $matrix_type eq 'auto' || $matrix_type eq 'logodds' || $matrix_type eq 'probability';

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
    parse_and_convert_json($fh);
} else {
    parse_and_convert_homer($fh);
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
    my ($row_ref) = @_;
    return 1 if $matrix_type eq 'logodds';
    return 0 if $matrix_type eq 'probability';
    my $sum = 0;
    $sum += $_ for @$row_ref;
    return ($sum < 0.98 || $sum > 1.02);
}

sub logodds_to_prob {
    my ($row_ref, $pc) = @_;
    my $background = 0.25;
    my @raw = map { 2 ** $_ * $background } @$row_ref;
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
    my ($id, $desc, $matrix_ref, $header_ref) = @_;

    my @mat = @$matrix_ref;

    if ($do_rc) {
        @mat = reverse_complement(\@mat, \$id);
    }

    if ($trim_edges > 0) {
        @mat = trim_edges(\@mat, $trim_edges);
        if (!@mat) {
            warn "Warning: motif '$id' trimmed to empty matrix (IC threshold=$trim_edges)\n";
            return;
        }
    }

    if ($min_ic > 0 && total_ic(\@mat) < $min_ic) {
        return;
    }

    print_meme_header() unless $$header_ref;
    $$header_ref = 1;
    print_meme_motif($id, $desc, \@mat);
}

sub print_meme_header {
    print "MEME version 4\n";
    print "\n";
    print "ALPHABET= ACGT\n";
    print "\n";
    print "strands: + -\n";
    print "\n";
    print "Background letter frequencies\n";
    print "A 0.25 C 0.25 G 0.25 T 0.25\n";
    print "\n";
}

sub print_meme_motif {
    my ($id, $desc, $matrix_ref) = @_;
    my $width = scalar @$matrix_ref;
    print "MOTIF $id $desc\n";
    print "\n";
    print "letter-probability matrix: alength= 4 w= $width nsites= 20 E= 0\n";
    foreach my $row (@$matrix_ref) {
        print "  " . join("  ", map { sprintf("%.6f", $_) } @$row) . "\n";
    }
    print "\n";
}

sub parse_and_convert_homer {
    my ($fh) = @_;
    my $header_printed = 0;
    my $in_motif = 0;
    my $motif_id = '';
    my $description = '';
    my @matrix;

    while (<$fh>) {
        chomp;
        next unless length($_);

        if (/^>(.*)/) {
            my $rest = $1;

            if ($in_motif && @matrix) {
                process_motif($motif_id, $description, \@matrix, \$header_printed);
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
        if (scalar(@row) != 4) {
            warn "Warning: skipping malformed matrix row (expected 4 cols, got "
                 . scalar(@row) . "): $_\n";
            next;
        }

        if (is_logodds(\@row)) {
            @row = logodds_to_prob(\@row, $pseudocount);
        }
        push @matrix, \@row;
    }

    if ($in_motif && @matrix) {
        process_motif($motif_id, $description, \@matrix, \$header_printed);
    }
}

sub parse_and_convert_json {
    my ($fh) = @_;
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
    my @motifs;
    if (ref $data eq 'HASH' && ref $data->{motifs} eq 'ARRAY') {
        for my $m (@{$data->{motifs}}) {
            my $id = $m->{id} || 'motif';
            my $desc = $m->{description} || $id;

            if ($extract && $id ne $extract && $desc ne $extract) {
                next;
            }

            my @matrix;
            if (ref $m->{matrix} eq 'ARRAY') {
                for my $row (@{$m->{matrix}}) {
                    if (ref $row eq 'ARRAY' && scalar(@$row) == 4) {
                        my @vals = map { $_ + 0 } @$row;
                        if (is_logodds(\@vals)) {
                            @vals = logodds_to_prob(\@vals, $pseudocount);
                        }
                        push @matrix, \@vals;
                    } else {
                        warn "Warning: skipping malformed matrix row (expected 4 cols)\n";
                    }
                }
            }

            if (@matrix) {
                push @motifs, { id => $id, desc => $desc, matrix => \@matrix };
            }
        }
    }

    if (@motifs) {
        for my $m (@motifs) {
            process_motif($m->{id}, $m->{desc}, $m->{matrix}, \$header_printed);
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
    -f, --format <fmt>  Input format: homer (default) or json
    --input-format <fmt>  Matrix type: auto (default), logodds, or probability
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
