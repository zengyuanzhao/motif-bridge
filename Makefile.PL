use ExtUtils::MakeMaker;

WriteMakefile(
    NAME         => 'motif-bridge',
    VERSION_FROM => 'perl_scripts/meme2homer.pl',
    AUTHOR       => 'Zengyuan Zhao',
    ABSTRACT     => 'Cross-language toolkit for converting motifs between MEME and HOMER formats',
    EXE_FILES    => [
        'perl_scripts/meme2homer.pl',
        'perl_scripts/homer2meme.pl',
    ],
    PREREQ_PM    => {
        'IO::Uncompress::Gunzip' => 0,
        'Getopt::Long'           => 0,
    },
    LICENSE      => 'mit',
);
