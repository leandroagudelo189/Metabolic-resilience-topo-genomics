#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/";  # Ensure the script can locate local modules
use Hypergeometric;
use Class::Struct;

=head1 NAME

pge_script - Positional Gene Enrichment (PGE) Analysis

=head1 SYNOPSIS

    perl pge_script.pl -q query_IDs_file [options]

    Options:
        -h, --help              Display this help message
        -l, --list              List available reference datasets
        -g, --allow-gaps        Allow gaps in enriched regions
        -q, --query FILE        Query IDs file (required)
        -r, --reference DATASET Reference dataset (default: ensembl42)
        -a, --alpha FLOAT       P-value significance threshold (default: 0.05)
        -s, --ratio FLOAT       Redundancy filter ratio (default: 1.5)
        -m, --mapped-ids TYPE   Map to which IDs [ensembl|symbols] (default: ensembl)
        -c, --chromosome CHR    Restrict search to a particular chromosome
        -p, --path DIR          Path to datasets (default: /data2/pge/data)

    Example:
        perl pge_script.pl -q my_query_file.txt -r ensembl42

=head1 DESCRIPTION

This script performs a Positional Gene Enrichment (PGE) analysis to identify regions in the genome that are significantly enriched with genes of interest.

=head1 AUTHOR

Your Name

=head1 LICENSE AND COPYRIGHT

This program is released under the GNU General Public License version 3.

=cut

# Define Data Structures
struct Dataset => {
    name       => '$',
    loadMethod => '$',
    file       => '$',
};

struct ID => {
    id        => '$',
    mappedID  => '$',
    symbol    => '$',
    ensg      => '$',
    chr       => '$',
    start     => '$',
    end       => '$',
    rank      => '$',
};

struct Chromosome => {
    name       => '$',
    ids        => '@',
    size       => '$',
    query      => '@',
    sEnriched  => '@',  # Significantly Enriched
    dEnriched  => '@',
    location   => '%',  # Map of rank to gene
};

struct Hit => {
    name       => '$',
    pvalue     => '$',
    pvalueadj  => '$',
    chr        => '$',
    start      => '$',
    end        => '$',
    common     => '$',
    size       => '$',
    seqStart   => '$',
    seqEnd     => '$',
    pertinent  => '$',
};

# Initialize Variables and Default Parameters
my $dataDir        = '/data2/pge/data';
my %refDatasets    = ();
my %humanDatasets  = ();
my %mouseDatasets  = ();

# Initialize Reference Datasets
initialize_datasets();

# Default Parameters
my $dataset         = 'ensembl42';
my $alpha           = 0.05;
my $summarize       = 3;
my $min_ratio       = 1.5;
my $mappedIDs       = 'ensembl';
my $restrict_chr    = undef;
my $allow_gaps      = 0;

# Command-Line Options Parsing
my %opt;
GetOptions(
    'h|help'            => \$opt{h},
    'l|list'            => \$opt{l},
    'g|allow-gaps'      => \$opt{g},
    'q|query=s'         => \$opt{q},
    'r|reference=s'     => \$opt{r},
    'a|alpha=f'         => \$opt{a},
    's|ratio=f'         => \$opt{s},
    'm|mapped-ids=s'    => \$opt{m},
    'c|chromosome=s'    => \$opt{c},
    'p|path=s'          => \$opt{p},
) or pod2usage(2);

pod2usage(1) if $opt{h};

# List datasets if requested
list_datasets() if $opt{l};

# Ensure query file is provided
pod2usage("Error: Query file is required.\n") unless $opt{q};
my $queryFile = $opt{q};

# Set parameters from options
$dataset     = lc($opt{r}) if $opt{r};
$alpha       = $opt{a} if defined $opt{a};
$min_ratio   = $opt{s} if defined $opt{s} && $summarize == 3;
$mappedIDs   = $opt{m} if defined $opt{m};
$restrict_chr = lc($opt{c}) if defined $opt{c};
$allow_gaps  = 1 if $opt{g};
$dataDir     = $opt{p} if defined $opt{p};

# Validate dataset
unless (exists $refDatasets{$dataset}) {
    warn "Unknown reference dataset '$dataset'.\n";
    list_datasets();
}

# Initialize Hypergeometric Object
my $hg = Hypergeometric->new();

# Initialize Chromosomes
my %chrs;
initialize_chromosomes();

# Load Query Elements
my @queryElm;
my $querySize = load_query_elements($queryFile, \@queryElm);

# Load Reference IDs
my %ids = load_reference_ids();

# Compute Ranks of Genes on Chromosomes
my %ids_ranks;
compute_gene_ranks(\%ids_ranks);

# Retrieve Query Element Locations
my %queryElmHash;
my $query_symbols = '';
retrieve_query_locations(\%queryElmHash, \$query_symbols, \%ids_ranks);

# Perform Enrichment Analysis
my (@results, @pvalues);
perform_enrichment_analysis(\@results, \@pvalues);

# Adjust P-values for False Discovery Rate (FDR)
adjust_pvalues(\@results);

# Filter Significant Hits
my @sresults;
my $nbSignificant = filter_significant_hits(\@results, \@sresults);

# Filter Redundant Hits
filter_redundant_hits();

# Generate Output
generate_output(\@sresults, $query_symbols);

exit;

# Subroutine Definitions

sub initialize_datasets {
    # Human Datasets
    $refDatasets{'hgu95av2'}        = Dataset->new(
        name       => 'hgu95av2',
        loadMethod => \&loadAffy,
        file       => "$dataDir/HG_U95Av2.na22.annot.txt"
    );
    $humanDatasets{'hgu95av2'}      = 1;

    $refDatasets{'hgu133a'}         = Dataset->new(
        name       => 'hgu133a',
        loadMethod => \&loadAffy,
        file       => "$dataDir/HG-U133A.na22.annot.txt"
    );
    $humanDatasets{'hgu133a'}        = 1;

    $refDatasets{'hgu133a_subset'}  = Dataset->new(
        name       => 'hgu133a_subset',
        loadMethod => \&loadAffy,
        file       => "$dataDir/hgu133a_subset.txt"
    );
    $humanDatasets{'hgu133a_subset'} = 1;

    $refDatasets{'ensembl42'}       = Dataset->new(
        name       => 'ensembl42',
        loadMethod => \&load_ID_CHR_START_END_DESC,
        file       => "$dataDir/ensembl42.txt"
    );
    $humanDatasets{'ensembl42'}      = 1;

    $refDatasets{'entrez'}          = Dataset->new(
        name       => 'entrez',
        loadMethod => \&load_ID_CHR_START_END_DESC,
        file       => "$dataDir/entrez.txt"
    );
    $humanDatasets{'entrez'}         = 1;

    $refDatasets{'symbols'}         = Dataset->new(
        name       => 'symbols',
        loadMethod => \&load_ID_CHR_START_END_DESC,
        file       => "$dataDir/symbols.txt"
    );
    $humanDatasets{'symbols'}        = 1;

    $refDatasets{'refseq_dna'}      = Dataset->new(
        name       => 'refseq_dna',
        loadMethod => \&load_ID_CHR_START_END_DESC,
        file       => "$dataDir/refseq_dna.txt"
    );
    $humanDatasets{'refseq_dna'}     = 1;

    $refDatasets{'refseq_peptide'}  = Dataset->new(
        name       => 'refseq_peptide',
        loadMethod => \&load_ID_CHR_START_END_DESC,
        file       => "$dataDir/refseq_peptide.txt"
    );
    $humanDatasets{'refseq_peptide'} = 1;

    # Mouse Datasets
    $refDatasets{'mouse430_2'}      = Dataset->new(
        name       => 'mouse430_2',
        loadMethod => \&loadAffy,
        file       => "$dataDir/Mouse430_2.na24.annot.csv"
    );
    $mouseDatasets{'mouse430_2'}    = 1;
}

sub list_datasets {
    print "Available datasets:\n";
    foreach my $dataset_name (sort keys %refDatasets) {
        print "\t- $dataset_name\n";
    }
    exit;
}

sub initialize_chromosomes {
    if ($restrict_chr) {
        $chrs{$restrict_chr} = Chromosome->new(name => $restrict_chr, size => 0);
    }
    else {
        if (exists $humanDatasets{$dataset}) {
            foreach my $i (1..22, 'x', 'y') {
                $chrs{$i} = Chromosome->new(name => $i, size => 0);
            }
        }
        elsif (exists $mouseDatasets{$dataset}) {
            foreach my $i (1..19, 'x', 'y') {
                $chrs{$i} = Chromosome->new(name => $i, size => 0);
            }
        }
    }
}

sub load_query_elements {
    my ($file, $queryElm_ref) = @_;
    open my $fh, '<', $file or die "Cannot open query file '$file': $!";
    my $querySize = 0;
    while (my $line = <$fh>) {
        chomp $line;
        my @elements = split /\s+/, $line;
        foreach my $element (@elements) {
            next unless length $element;
            push @$queryElm_ref, lc($element);
            $querySize++;
        }
    }
    close $fh;
    return $querySize;
}

sub load_reference_ids {
    my $load_method = $refDatasets{$dataset}->loadMethod;
    my $file = $refDatasets{$dataset}->file;
    my %ids = $load_method->($file, $mappedIDs);
    return %ids;
}

sub loadAffy {
    my ($filename, $mappedIDs) = @_;
    my %ids;
    my %alreadyMappedIDs;
    open my $fh, '<', $filename or die "Cannot open file '$filename': $!";
    while (my $line = <$fh>) {
        next if $. == 1;
        chomp $line;
        my @fields = split /","/, $line;
        my ($refID) = $fields[0] =~ /.*"(\S+)/;
        $refID = lc($refID);
        my $loc = $fields[12];
        next if $loc eq '---';
        next if index($loc, "///") != -1;
        my ($chr, $start, $end) = $loc =~ /chr(\S+):(\d+)-(\d+)/;
        $chr = lc($chr);
        next unless exists $chrs{$chr};
        my $ensg = $fields[17];
        my $symbol = $fields[14] eq '---' ? $ensg : $fields[14];
        my $mappedID = $mappedIDs eq 'symbols' ? $symbol : $ensg;
        next if $mappedID eq '---' || index($mappedID, '///') != -1;
        my $id = ID->new(
            id       => $refID,
            mappedID => $mappedID,
            symbol   => $symbol,
            ensg     => $ensg,
            chr      => $chr,
            start    => $start,
            end      => $end
        );
        unless (exists $alreadyMappedIDs{$mappedID}) {
            push @{$chrs{$chr}->ids}, $id;
            $alreadyMappedIDs{$mappedID} = 1;
        }
        push @{$ids{$refID}}, $id;
    }
    close $fh;
    return %ids;
}

sub load_ID_CHR_START_END_DESC {
    my ($filename) = @_;
    my %ids;
    open my $fh, '<', $filename or die "Cannot open file '$filename': $!";
    while (my $line = <$fh>) {
        chomp $line;
        my ($refID, $chr, $start, $end, $desc) = split /\t/, $line;
        $refID = lc($refID);
        $chr = lc($chr);
        next unless exists $chrs{$chr};
        my $id = ID->new(
            id       => $refID,
            mappedID => $refID,
            symbol   => $refID,
            ensg     => $refID,
            chr      => $chr,
            start    => $start,
            end      => $end
        );
        push @{$chrs{$chr}->ids}, $id;
        push @{$ids{$refID}}, $id;
    }
    close $fh;
    return %ids;
}

sub compute_gene_ranks {
    my ($ids_ranks_ref) = @_;
    foreach my $chr_obj (values %chrs) {
        my $rank = 1;
        my %loc;
        foreach my $gene (sort { $a->start <=> $b->start } @{$chr_obj->ids}) {
            unless (exists $ids_ranks_ref->{$gene->mappedID}) {
                $gene->rank($rank);
                $loc{$rank} = $gene;
                $ids_ranks_ref->{$gene->mappedID} = $rank;
                $rank++;
            }
            else {
                $gene->rank($ids_ranks_ref->{$gene->mappedID});
            }
        }
        $chr_obj->size($rank - 1);
        $chr_obj->location(\%loc);
    }
}

sub retrieve_query_locations {
    my ($queryElmHash_ref, $query_symbols_ref, $ids_ranks_ref) = @_;
    foreach my $element (@queryElm) {
        if (!exists $ids{$element}) {
            next;
        }
        foreach my $id_obj (@{$ids{$element}}) {
            next if exists $queryElmHash_ref->{$id_obj->mappedID};
            push @{$chrs{$id_obj->chr}->query}, $id_obj;
            $queryElmHash_ref->{$id_obj->mappedID} = 1;
            $id_obj->rank($ids_ranks_ref->{$id_obj->mappedID}) unless defined $id_obj->rank;
            $$query_symbols_ref .= join(',', $element, $id_obj->mappedID, $id_obj->chr, $id_obj->start, $id_obj->end, $id_obj->ensg) . "\n";
        }
    }
}

sub perform_enrichment_analysis {
    my ($results_ref, $pvalues_ref) = @_;
    if ($allow_gaps) {
        foreach my $chr_obj (values %chrs) {
            my @query_ranks = map { $_->rank } sort { $a->rank <=> $b->rank } @{$chr_obj->query};
            my @chr_results = compare_pertinent_target_set(\@query_ranks, $querySize, scalar keys %ids, $chr_obj->name);
            push @$results_ref, @chr_results;
            push @$pvalues_ref, map { $_->pvalue } @chr_results;
        }
    }
    else {
        my $avg_gap = scalar(keys %ids) / $querySize;
        foreach my $chr_obj (values %chrs) {
            my @c_query = sort { $a->rank <=> $b->rank } @{$chr_obj->query};
            next unless @c_query > 1;
            my @queryInt;
            my $lastRank = 0;
            foreach my $gene (@c_query) {
                if (!$lastRank || $gene->rank - $lastRank > $avg_gap) {
                    if (@queryInt > 1) {
                        my @chr_results = compare_pertinent_target_set(\@queryInt, $querySize, scalar keys %ids, $chr_obj->name);
                        push @$results_ref, @chr_results;
                        push @$pvalues_ref, map { $_->pvalue } @chr_results;
                    }
                    @queryInt = ();
                }
                push @queryInt, $gene->rank;
                $lastRank = $gene->rank;
            }
            if (@queryInt > 1) {
                my @chr_results = compare_pertinent_target_set(\@queryInt, $querySize, scalar keys %ids, $chr_obj->name);
                push @$results_ref, @chr_results;
                push @$pvalues_ref, map { $_->pvalue } @chr_results;
            }
        }
    }
}

sub compare_pertinent_target_set {
    my ($qset_ref, $q_size, $popsize, $chr_name) = @_;
    my @loc = @$qset_ref;
    my @results;
    my $qschr = scalar @loc;
    my $nchr = scalar @{$chrs{$chr_name}->ids};
    for (my $i = 0; $i < $qschr - 1; $i++) {
        next if ($i > 0 && $loc[$i - 1] + 1 == $loc[$i]);
        for (my $j = $i + 1; $j < $qschr; $j++) {
            next if ($j < $qschr - 1 && $loc[$j] + 1 == $loc[$j + 1]);
            my $t = $loc[$j] - $loc[$i] + 1;
            my $c = $j - $i + 1;
            my $pvalue = $hg->compute($popsize, $q_size, $t, $c);
            my $hit = Hit->new(
                name       => "$chr_name,$loc[$i]:$loc[$j]",
                pvalue     => $pvalue,
                chr        => $chr_name,
                start      => $loc[$i],
                end        => $loc[$j],
                common     => $c,
                size       => $t,
                pertinent  => 1,
            );
            push @results, $hit;
        }
    }
    return @results;
}

sub adjust_pvalues {
    my ($results_ref) = @_;
    my @hits = sort { $a->pvalue <=> $b->pvalue } @$results_ref;
    my $n = scalar @hits;
    if ($n > 0) {
        $hits[$n - 1]->pvalueadj($hits[$n - 1]->pvalue);
        for (my $i = $n - 2; $i >= 0; $i--) {
            my $adj_pvalue = min($hits[$i + 1]->pvalueadj, ($n / ($i + 1)) * $hits[$i]->pvalue, 1);
            $hits[$i]->pvalueadj($adj_pvalue);
        }
    }
}

sub filter_significant_hits {
    my ($results_ref, $sresults_ref) = @_;
    my $nbSignificant = 0;
    foreach my $hit (sort { $a->pvalueadj <=> $b->pvalueadj } @$results_ref) {
        if ($hit->pvalueadj <= $alpha) {
            push @$sresults_ref, $hit;
            $nbSignificant++;
            push @{$chrs{$hit->chr}->sEnriched}, $hit;
        }
    }
    return $nbSignificant;
}

sub filter_redundant_hits {
    foreach my $chr_obj (values %chrs) {
        my @hits = sort { $a->start <=> $b->start } @{$chr_obj->sEnriched};
        my $nb = scalar @hits;
        for (my $i = 0; $i < $nb - 1; $i++) {
            my $h1 = $hits[$i];
            for (my $j = $i + 1; $j < $nb && $hits[$j]->start < $h1->end; $j++) {
                my $h2 = $hits[$j];
                if ($h2->end <= $h1->end) {
                    if ($h1->pvalue <= $h2->pvalue &&
                        $min_ratio * $h1->common / $h1->size > $h2->common / $h2->size) {
                        $h2->pertinent(0);
                    }
                    if ($h1->pvalue >= $h2->pvalue) {
                        $h1->pertinent(0);
                    }
                }
                elsif ($h1->start == $h2->start) {
                    if ($h1->end <= $h2->end) {
                        if ($h2->pvalue <= $h1->pvalue &&
                            $min_ratio * $h2->common / $h2->size > $h1->common / $h1->size) {
                            $h1->pertinent(0);
                        }
                        if ($h2->pvalue >= $h1->pvalue) {
                            $h2->pertinent(0);
                        }
                    }
                }
            }
        }
    }
}

sub generate_output {
    my ($sresults_ref, $query_symbols) = @_;
    my $bed = "track name=\"PGE\" description=\"Positional Gene Enrichment\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n";
    my $raw = '';
    my $adjustmentMethod = 'FDR';

    foreach my $chr_obj (values %chrs) {
        my @hits = sort { $a->pvalue <=> $b->pvalue } @{$chr_obj->sEnriched};
        foreach my $hit (@hits) {
            next unless $hit->pertinent;
            my $start = $chr_obj->location->{$hit->start}->start;
            my $end   = $chr_obj->location->{$hit->end}->end;
            $bed .= "chr" . uc($chr_obj->name) . "\t$start\t$end\t$hit->pvalue\n";
            $raw .= uc($chr_obj->name) . "\t$start\t$end\t$hit->pvalue\t$hit->pvalueadj\t$hit->common\t$hit->size\n";
        }
    }

    print "<mapping>$query_symbols</mapping><raw>$raw</raw><bed>$bed</bed><alpha>$alpha</alpha><adjust>$adjustmentMethod</adjust><reference>$dataset</reference>";
}

sub min {
    my @values = @_;
    my $min = shift @values;
    foreach (@values) {
        $min = $_ if $_ < $min;
    }
    return $min;
}
