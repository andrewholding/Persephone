use strict;

package Persephone::Proteins;

use base 'Exporter';
use lib 'lib';
use Persephone::Config;
use Persephone::Data;

our @EXPORT = ( 'protein_residuemass',         'calculate_peptide_masses',
		'calculate_crosslink_peptides', 'generate_monolink_peptides',
		'generate_modified_peptides',   'digest_proteins', 'modifications');

sub modifications {

    my ($table, $settings_ref) = @_;

    my %settings = %{$settings_ref};
    my $mono_mass_diff = $settings{'monolink_mass'};
    my $xlinker_mass = $settings{'linker_mass'};
    my $reactive_site = $settings{'reactive_amino_acid'};
    my $dbh = connect_db(\%settings);
    
        if ($reactive_site =~ /[^,]/) {  $reactive_site = $reactive_site . ',' . $reactive_site};
    my @reactive_sites = split (',',$reactive_site);

    my %modifications = (

        MonoLink => {
                      Name     => 'mono link',
                      Location => $reactive_sites[0],
                      Delta    => $mono_mass_diff,
        },
        LoopLink => {
            Name     => 'loop link',
            Location => $reactive_sites[1],
            Delta    => -18.0105646
            ,    #Loop links are treated as a modified monolink (loop link on an xlink is too complicated, and wierd)
        },
        NoMod => {
                   Name     => ' ',
                   Location => '^.',
                   Delta    => 0,
        },
    );



    if (defined $table) {
        my $dynamic_mods = get_mods($table, 'dynamic', $dbh);
        while ((my $dynamic_mod = $dynamic_mods->fetchrow_hashref)) {
            $modifications{ $dynamic_mod->{'mod_id'} }{'Name'}     = $dynamic_mod->{'mod_name'};
            $modifications{ $dynamic_mod->{'mod_id'} }{'Location'} = $dynamic_mod->{'mod_residue'};
            $modifications{ $dynamic_mod->{'mod_id'} }{'Delta'}    = $dynamic_mod->{'mod_mass'};
        }
    }

    return %modifications;
}

sub generate_modified_peptides {

my ($results_table,$settings_ref) = @_;

  print "Calculating modified peptides ... \n";

my %settings = %{$settings_ref};
my %modifications       = modifications($results_table, \%settings);


foreach my $modification (sort(keys %modifications)) {


		    my $results_dbh = connect_db_results($results_table,0, \%settings);


		    if (   !($modifications{$modification}{Name} eq "loop link" )
			&& !($modifications{$modification}{Name} eq "mono link" ) 
			&& !($modifications{$modification}{Name} eq " ")
		      ) 
		    {
			for (my $n = 1; $n <= 3; $n++) {
			      my $modify = $results_dbh->prepare("
				    INSERT INTO peptides
				    SELECT 
					  null as rowid,
					  results_table as results_table,
					  sequence as sequence,
					  source as source,
					  linear_only as linear_only,
					  mass + ? as mass, 
					  ? as modifications,
					  monolink as monolink,
					  xlink as xlink,
					  ? as no_of_mods
					  FROM peptides
						    WHERE modifications = '' and  (LENGTH(sequence) - LENGTH(REPLACE(sequence, ?, ''))) >= (? + 0)
					  ");
			$modify->execute($modifications{$modification}{Delta}*$n,$modification, $n, $modifications{$modification}{Location}, $n+0)
			}
		    } elsif ($modifications{$modification}{Name} eq "loop link" ) {
			      my $monolinks = $results_dbh->prepare("
				    INSERT INTO peptides
				    SELECT 
					  null as rowid,
					  results_table as results_table,
					  sequence as sequence,
					  source as source,
					  linear_only as linear_only,
					  mass + ? as mass, 
					  ? as modifications,
					  monolink as monolink,
					  xlink as xlink,
					  ? as no_of_mods
					  FROM peptides
						    WHERE  sequence LIKE ? and xlink = 0 and monolink > 0 and modifications = ''
			      ");
			$monolinks->execute($modifications{$modification}{Delta},$modification, 1, "%".$modifications{$modification}{Location}."%");
		  }
		  $results_dbh->commit;
	       }    
}

sub generate_monolink_peptides {

  my ( $results_table, $settings_ref ) = @_;

  print "Calculating monolink ... \n";

  my %settings = %{$settings_ref};

  my $reactive_site = $settings{'reactive_amino_acid'};
  my @monolink_masses = split( ",", $settings{monolink_mass} );
  my $results_dbh = connect_db_results( $results_table, 1, \%settings );

  if ( length($reactive_site) < 2 ) {
    $reactive_site = $reactive_site . "," . $reactive_site;
  }

  my $monolinks = $results_dbh->prepare( "
	INSERT INTO peptides
	SELECT 
		 null as rowid,
		 results_table as results_table,
		 sequence as sequence,
		 source as source,
		 linear_only as linear_only,
		 mass + ? as mass, 
		 '' as modifications,
		 ? as monolink,
		 0 as xlink,
		 0 as no_of_mods
	FROM peptides
 	WHERE sequence LIKE ? and xlink = 0 and monolink = 0 and results_table = ?;
    " );

  my @reactive_sites = split( ',', $reactive_site );

  foreach my $monolink_mass (@monolink_masses) {
    $monolinks->execute( $monolink_mass, $monolink_mass,
      "%" . $reactive_sites[0] . "%",
      $results_table );
  }

  return;
}

sub calculate_crosslink_peptides {

  my ( $results_table, $settings_ref ) = @_;

  my %settings = %{$settings_ref};
  my $xlink;
  my $xlink_fragment_mass;
  my $xlink_fragment_sources;
  my $results_dbh = connect_db_results( $results_table, 1, \%settings );
  my $settings_dbh = connect_db( \%settings );

  my $reactive_site      = $settings{'reactive_amino_acid'};
  my $min_peptide_length = $settings{'enzyme_min_peptide_length'};
  my $xlinker_mass       = $settings{'linker_mass'};
  my $missed_clevages    = $settings{'max_missed_cleavages'};
  my $cut_residues     = $settings{'enzyme_cut_residues'};

  my $stop_duplicates = 'AND p1.rowid >= p2.rowid';

  if ( $reactive_site !~ m/,/ ) {
    $reactive_site = $reactive_site . "," . $reactive_site;
  }

  my @reactive_sites = split( ',', $reactive_site );

  if ( $reactive_sites[0] ne $reactive_sites[1] ) {
    $stop_duplicates = '';
  }

  my $peptidelist;

  my $index = $results_dbh->prepare(
    "CREATE INDEX peptide_index ON peptides (sequence(15));");
  $index->execute();

  my $sequencelist = $results_dbh->prepare(
    "select distinct source from peptides where xlink = 0");
  $sequencelist->execute;
  my $rows  = $sequencelist->rows;
  my $count = 0;

  my $threads = 0;
  $threads = 1;

  # 		 Threads = 1 to avoid a deadlock issue. Needs to be fixed
  my $pm = Parallel::ForkManager->new($threads);

  while ( my $source = $sequencelist->fetchrow_hashref ) {
    print "Crosslinking peptides  ... ",
      sprintf( "%.0f", ( $count / $rows ) * 100 ), "%\n";
    $count++;
    $pm->start and next;

    my $results_dbh_fork =
      connect_db_results( $results_table, 0, \%settings );

    $peptidelist = $results_dbh_fork->prepare( "
	      INSERT INTO peptides
	      SELECT
		    null as rowid,
		    p1.results_table as results_table,
		    concat (p1.sequence, '-', p2.sequence) as sequence,
		    concat (p1.source , '-' , p2.source) as source,
		    0 as linear_only,
		    p1.mass + p2.mass + ? as mass, 
		    '' as modifications,
		    0 as monolink,
		    1 as xlink,
		    0 as no_of_mods
	    
		FROM peptides p1 inner join peptides p2 on (p1.results_table = p2.results_table)
		WHERE p1.source = ? and p1.linear_only = '0' AND p2.linear_only = '0' AND p1.xlink ='0' and p2.xlink = '0' AND p1.sequence LIKE ? AND p2.sequence LIKE ?
		$stop_duplicates
	" );

    foreach my $reactive_site_chain_1 ( split //, $reactive_sites[0] ) {
      foreach my $reactive_site_chain_2 ( split //, $reactive_sites[1] )
      {

# 	warn "%".$reactive_site_chain_1."_%","%".$reactive_site_chain_2."_%";
	$peptidelist->execute(
	  $xlinker_mass, $source->{'source'},
	  "%" . $reactive_site_chain_1 . "_%",
	  "%" . $reactive_site_chain_2 . "_%"
	);
	$results_dbh_fork->commit;
      }
    }
    $peptidelist->finish;
    $results_dbh_fork->disconnect;

    $pm->finish;
  }
  $pm->wait_all_children;

  $sequencelist->finish;
  $results_dbh->disconnect;

}

sub protein_residuemass {

  my ( $table, $dbh ) = @_;

  my %protein_residuemass = (
    G => 57.02146,
    A => 71.03711,
    S => 87.03203,
    P => 97.05276,
    V => 99.06841,
    T => 101.04768,
    C => 103.00919,
    L => 113.08406,
    I => 113.08406,
    X => 113.08406,    # (L or I)
    N => 114.04293,
    O => 114.07931,
    B => 114.53494,    # (avg N+D)
    D => 115.02694,
    Q => 128.05858,
    K => 128.09496,
    Z => 128.55059,    #(avg Q+E)
    E => 129.04259,
    M => 131.04049,
    H => 137.05891,
    F => 147.06841,
    R => 156.10111,
    Y => 163.06333,
    W => 186.07931
  );

  if ( defined $table ) {
    my $fixed_mods = get_mods( $table, 'fixed', $dbh );
    while ( ( my $fixed_mod = $fixed_mods->fetchrow_hashref ) ) {
      $protein_residuemass{ $fixed_mod->{'mod_residue'} } =
	$protein_residuemass{ $fixed_mod->{'mod_residue'} } +
	$fixed_mod->{'mod_mass'};
    }

  }
  return %protein_residuemass;
}


sub digest_proteins    #Digest a string into an array of peptides
{
  my ( $missed_clevages, $protein_sequences, $cut_residues,
    $nocut_residues, $n_or_c )
    = @_;
  my @protein_fragments;
  $protein_sequences =~ s/[^A-Z]//g;
  my $protein = $protein_sequences;
  $protein =~ s/[^\w\d>]//g;
  if ( $nocut_residues eq '' ) { $nocut_residues = '9' }
  ; #Numbers don't appear in sequences so this just works, easier than having a second regex

#   print "No Cut:", $nocut_residues, "\n";
  my @digest;
  my @digest_not_for_crosslinking;
  if ( $n_or_c eq 'C' ) {

#     print "Protease type:", $n_or_c, "\n";
    @digest = $protein =~
m/(?:(?:[^$cut_residues]|[$cut_residues]$nocut_residues)*(?:[$cut_residues](?!$nocut_residues)|.(?=$))){1}/g;
    my @single_digest = @digest;
    for ( my $i = 2 ; $i < ( $missed_clevages * 2 + 1 ) + 2 ; $i++ ) {
      my @single_digest_trimmed = @single_digest
	; #need to include missed cleavages for each possible missed position
      my @parts = $protein =~
m/(?:(?:[^$cut_residues]|[$cut_residues]$nocut_residues)*(?:[$cut_residues](?!$nocut_residues)|.(?=$))){$i}/g;
      if ( $i < $missed_clevages + 2 ) {
	push( @digest, @parts );
      }
      else {
	push( @digest_not_for_crosslinking, @parts );
      }
      for ( my $j = 1 ; $j < $i ; $j++ ) {
	shift @single_digest_trimmed;
	@parts =
	  join( "", @single_digest_trimmed ) =~
m/(?:(?:[^$cut_residues]|[$cut_residues]$nocut_residues)*(?:[$cut_residues](?!$nocut_residues)|.(?=$))){$i}/g;
	if ( $i < $missed_clevages + 2 ) {
	  push( @digest, @parts );
	}
	else {
	  push( @digest_not_for_crosslinking, @parts );
	}
      }
    }

  }
  else {

    print "Protease type:", $n_or_c, "\n";
    @digest = $protein =~
m/(?:(?:[$cut_residues](?!$nocut_residues)|^.)(?:[^$cut_residues]|[$cut_residues]$nocut_residues)*){1}/g;
    my @single_digest = @digest;
    for ( my $i = 2 ; $i < ( $missed_clevages * 2 + 1 ) + 2 ; $i++ ) {
      my @single_digest_trimmed = @single_digest;
      my @parts                 = $protein =~
m/(?:(?:[$cut_residues](?!$nocut_residues)|^.)(?:[^$cut_residues]|[$cut_residues]$nocut_residues)*){1}/g;
      if ( $i < $missed_clevages + 2 ) {
	push( @digest, @parts );
      }
      else {
	push( @digest_not_for_crosslinking, @parts );
      }
      for ( my $j = 1 ; $j < $i ; $j++ ) {
	shift @single_digest_trimmed;
	@parts =
	  join( "", @single_digest_trimmed ) =~
m/(?:(?:[$cut_residues](?!$nocut_residues)|^.)(?:[^$cut_residues]|[$cut_residues]$nocut_residues)*){1}/g;
	if ( $i < $missed_clevages + 2 ) {
	  push( @digest, @parts );
	}
	else {
	  push( @digest_not_for_crosslinking, @parts );
	}
      }
    }
  }
  return \@digest, \@digest_not_for_crosslinking;
}



sub digest_proteins_masses    #Calculates the mass of a list of peptides
{
  my ( $protein_fragments_ref, $protein_residuemass_ref,
    $fragment_source_ref )
    = @_;
  my @protein_fragments   = @{$protein_fragments_ref};
  my %protein_residuemass = %{$protein_residuemass_ref};
  my %fragment_source     = %{$fragment_source_ref};
  my $peptide_mass        = 0;
  my $terminalmass        = 1.0078250 * 2 + 15.9949146 * 1;
  my %protein_fragments_masses;
  my @protein_fragments_masses;

  foreach my $peptide (@protein_fragments) {
    if ( $peptide =~ /[ARNDCEQGHILKMFPSTWYV]/ ) {
      my @residues = split //, $peptide;

      foreach my $residue (@residues)
      {    #split the peptide in indivual amino acids
	$peptide_mass =
	  $peptide_mass +
	  $protein_residuemass{$residue
	  };    #tally the masses of each amino acid one at a time
      }

      $protein_fragments_masses{$peptide} =
	$peptide_mass + $terminalmass;

#  	 warn "," ,$peptide," , " ,$peptide_mass +$terminalmass ," , " ,($peptide_mass +$terminalmass)/2 , ",",($peptide_mass +$terminalmass)/3 , "\n";
      $peptide_mass = 0;
    }

  }

  return %protein_fragments_masses;
}

sub calculate_peptide_masses {

  my ( $results_table, $settings_ref ) = @_;
  my %settings            = %{$settings_ref};
  my $peptide_mass        = 0;
  my $terminalmass        = 1.0078250 * 2 + 15.9949146 * 1;
  my $count               = 0;

  print "Calulating masses...  \n";

  my $results_dbh = connect_db_results( $results_table, 1, \%settings );
  my $settings_dbh = connect_db( \%settings );
  
  my %protein_residuemass = protein_residuemass($results_table, $settings_dbh);
  
  my $sequencelist = $results_dbh->prepare(
    "select distinct source from peptides where xlink = 0");
  $sequencelist->execute;
  my $rows = $sequencelist->rows;

  my $threads = 0;
  $threads = no_of_threads;
  my $pm = Parallel::ForkManager->new($threads);

  while ( my $source = $sequencelist->fetchrow_hashref ) {
    print "Calculating masses ... ",
      sprintf( "%.0f", ( $count / $rows ) * 100 ), "%\n";
    $count++;

    $pm->start and next;    # do the fork

    my $results_dbh_fork =
      connect_db_results( $results_table, 0, \%settings );
    my $peptidelist = $results_dbh_fork->prepare(
      "SELECT * FROM peptides WHERE results_table = ? AND source = ?");
    my $update_mass = $results_dbh_fork->prepare(
"UPDATE peptides SET mass = ? WHERE  results_table = ? AND sequence = ?;"
    );

    $peptidelist->execute( $results_table, $source->{'source'} );

    while ( my $peptides = $peptidelist->fetchrow_hashref ) {
      my $peptide = $peptides->{'sequence'};
      if ( $peptide =~ /[ARNDCEQGHILKMFPSTWYV]/ ) {
	my @residues = split //, $peptide;

	foreach my $residue (@residues)
	{    #split the peptide in indivual amino acids
	  $peptide_mass =
	    $peptide_mass +
	    $protein_residuemass{$residue
	    };    #tally the masses of each amino acid one at a time
	}

#          $protein_fragments_masses{$peptide} = $peptide_mass + $terminalmass;
	$update_mass->execute( $peptide_mass + $terminalmass,
	  $results_table, $peptide );
	$peptide_mass = 0;
      }

    }
    $peptidelist->finish;
    $update_mass->finish;
    $results_dbh_fork->commit;
    $results_dbh_fork->disconnect;
    $pm->finish;
  }
  $pm->wait_all_children;
  $sequencelist->finish;
  $results_dbh->disconnect;
}




1;