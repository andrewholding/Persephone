#!/usr/bin/perl -w

use strict;
use DBI;
use Parallel::ForkManager;

use lib 'lib';
use Persephone::Data;
use Persephone::Config;



# Load settings from experiment definition file
my $experiment_definition_file = $ARGV[0];
my $settings_ref = import_experiment_definition_file($experiment_definition_file);
my %settings = %{$settings_ref};

#Load data from Fasta file into database and store peptide sequences
my $experiment_id = import_fasta(\%settings);

#Calculate linear peptide masses
my $settings_dbh = connect_db(\%settings);
my %protein_residuemass = protein_residuemass($experiment_id, $settings_dbh);
calculate_peptide_masses($experiment_id, \%protein_residuemass, \%settings);

#Cross-link peptides and calculate masses 
calculate_crosslink_peptides( $experiment_id,  \%settings);

#Generate mono-links and calculate masses
generate_monolink_peptides($experiment_id, \%settings);

#Generate Modified Peptides
generate_modified_peptides($experiment_id, \%settings);