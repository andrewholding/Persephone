#!/usr/bin/perl -w

use strict;
use DBI;
use Parallel::ForkManager;

use lib 'lib';
use Persephone::Data;
use Persephone::Config;
use Persephone::Import;
use Persephone::Proteins;
use Persephone::Scoring;

# Load settings from experiment definition file
my $experiment_definition_file = $ARGV[0];
my $settings_ref = import_experiment_definition_file($experiment_definition_file);
my %settings = %{$settings_ref};

# Save settings - currently saves in something that approximates to hekate's
# settings, in future is should just store the settings file.
my $experiment_id = import_settings(\%settings);

# Load data from Fasta file into database and store peptide sequences
import_fasta($experiment_id, \%settings);

# Digest sequences in protein table to generate peptide table
import_peptides($experiment_id, \%settings);

# Calculate linear peptide masses
calculate_peptide_masses($experiment_id, \%settings);

# Cross-link peptides and calculate masses 
calculate_crosslink_peptides( $experiment_id,  \%settings);

# Generate mono-links and calculate masses
generate_monolink_peptides($experiment_id, \%settings);

# Generate Modified Peptides
generate_modified_peptides($experiment_id, \%settings);

# Import data from mgf
import_mgf($experiment_id, \%settings);

# Find doublets in mass spec data
loaddoubletlist_db($experiment_id, \%settings);

# Match doublets to peptides and score them
matchpeaks($experiment_id, \%settings);

calculate_fdr($experiment_id, \%settings);