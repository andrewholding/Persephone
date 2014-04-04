use strict;

package Persephone::Data;

use base 'Exporter';
use lib 'lib';
use Persephone::Config;

our @EXPORT = (
  'save_settings',                'connect_db',
  'connect_db_results',		   'create_peptide_table', 
  
  'loaddoubletlist_db',	
  
  'add_peptide', 		  'get_mods',
);


sub loaddoubletlist_db    #Used to get mass-doublets from the data.
{

    my ($results_id, $settings_ref) = @_;

    my %settings = %{$settings_ref};
    
    my $dbh = connect_db_results($results_id, 0, \%settings);
    
    print "Run $results_id: Finding doublets...  \n";
 
    
    my $doublet_ppm_err  = $settings{'doublet_tollerance'};
    my $linkspacing = $settings{'number_labelled_atoms'};
    my $isotope = $settings{'isotope'};
    my $scan_width =  $settings{'max_scan_seperation'};
    my $match_charge =  $settings{'charge_match'};
    my $match_intensity =  $settings{'match_intensity'};
    my $ms1_intensity_ratio = $settings{'match_intensity_ratio'};
    
    
    my $constants_ref = load_constants;
    my %constants = %{$constants_ref};
    
    my $mass_of_deuterium = $constants{'mass_of_deuterium'};
    my $mass_of_hydrogen = $constants{'mass_of_hydrogen'};
    my $mass_of_carbon13 = $constants{'mass_of_carbon13'};
    my $mass_of_carbon12 = $constants{'mass_of_carbon12'};
    

    my $mass_seperation = 0;
    if ($isotope eq "deuterium") {
        $mass_seperation = $linkspacing * ($mass_of_deuterium - $mass_of_hydrogen);
    } elsif ($isotope eq "carbon-13") {
        $mass_seperation = $linkspacing * ($mass_of_carbon13 - $mass_of_carbon12);
    }

    my $average_peptide_mass  = 750;
    my $mass_seperation_upper = $mass_seperation + $average_peptide_mass * (0 + ($doublet_ppm_err / 1000000));
    my $mass_seperation_lower = $mass_seperation + $average_peptide_mass * (0 - ($doublet_ppm_err / 1000000));
    my $isopairs;
    my @peaklist;

    my  $masslist;

    $masslist = $dbh->do("alter table  msdata add index (monoisotopic_mw)");	  
	

    if ($isotope ne "none") {
        my $charge_match_string = "";
        if ($match_charge == "1") {
            $charge_match_string = "and d1.charge = d2.charge ";
        }
        my $intensity_match_string = "";
        if ($match_intensity == "1") {
            $intensity_match_string = "and (d1.abundance > d2.abundance * ? and d1.abundance < d2.abundance * ?)  ";
        }
        $masslist = $dbh->prepare(
            "create table doublets as select * from (SELECT d1.*,
				  d2.scan_num as d2_scan_num,
				  d2.mz as d2_mz,
				  d2.MSn_string as d2_MSn_string,
				  d2.charge as d2_charge,
				  d2.monoisotopic_mw as d2_monoisotopic_mw,
				  d2.title as d2_title,
				  d2.abundance as d2_abundance,
				  d2.precursor_scan as d2_precursor_scan
			  FROM msdata d1 inner join msdata d2 on (d2.monoisotopic_mw between d1.monoisotopic_mw + ? and d1.monoisotopic_mw + ? )
				  and d2.scan_num between d1.scan_num - ? 
				  and d1.scan_num + ? " . $charge_match_string . $intensity_match_string . "and d1.fraction = d2.fraction 
				  and d1.msorder = 2 and d2.msorder = 2
			  ORDER BY d1.scan_num ASC) as doublets"
        );
        
#          create table doublets as select sequence from (select * from peptides) as doubles limit 5;

        

#                warn "Exceuting Doublet Search\n";
        if ($match_intensity == "1") {
            if ($ms1_intensity_ratio == '0' or !defined $ms1_intensity_ratio) { $ms1_intensity_ratio = 1 }
            
                $masslist->execute($mass_seperation_lower, $mass_seperation_upper, $scan_width, $scan_width,
                                   $ms1_intensity_ratio, 1 / $ms1_intensity_ratio);
           
        } else {
             $masslist->execute($mass_seperation_lower, $mass_seperation_upper, $scan_width, $scan_width) ;
        }

        #       warn "Finished Doublet Search\n";
    } else {
        $masslist = $dbh->prepare(
            "create table doublets SELECT *
  			  FROM msdata 
  			  ORDER BY scan_num ASC "
        );
         $masslist->execute() ;
    }
    my $settings_dbh = connect_db(\%settings);
    my $doublets_found = $dbh->selectrow_array('select count(*) from doublets');
    my $settings_sql = $settings_dbh->prepare("UPDATE settings SET doublets_found = ? WHERE  name = ?;");
    $settings_sql->execute($doublets_found, $results_id);

     print "Doublets found: $doublets_found\n";
    
    $dbh->disconnect;

}

sub create_settings {

  my ($settings_dbh) = @_;

  my $row_id_type = '';

  $row_id_type = "INTEGER NOT NULL AUTO_INCREMENT, PRIMARY KEY (name)";

  $settings_dbh->do(
    "CREATE TABLE IF NOT EXISTS settings (
						      name " . $row_id_type . " ,
						      description TEXT,
						      cut_residues TEXT,
						      protein_sequences TEXT,
						      reactive_site TEXT,
						      mono_mass_diff TEXT,
						      xlinker_mass TEXT,
						      decoy TEXT,
						      ms2_da TEXT,
						      ms1_ppm FLOAT,
						      finished FLOAT,
						      isotoptic_shift TEXT,
						      threshold TEXT,
						      doublets_found TEXT,
						      charge_match NUMERIC,
						      intensity_match NUMERIC,
						      scored_ions TEXT,
						      amber TEXT,
						      time TEXT,
						      non_specific_digest NUMERIC,
						      no_enzyme_min TEXT,
						      no_enzyme_max TEXT,		
						      use_previous INTEGER
						) "
  );

}

sub save_settings {

  my ( $settings_dbh, $settings_ref ) = @_;

  my %settings = %{$settings_ref};

#     if (!defined $amber_codon) {$amber_codon = 0;};
#
#
#     if   (defined $match_charge && $match_charge == '1') { $match_charge = 'Yes' }
#     else                        { $match_charge = 'No' }
#     if   (defined $match_intensity && $match_intensity == '1') { $match_intensity = 'Yes' }
#     else                           { $match_intensity = 'No' }

  my $scored_ions = '';
  if ( defined $settings{'score_aions'} ) {
    $scored_ions = $scored_ions . 'A-ions, ';
  }
  if ( defined $settings{'score_bions'} ) {
    $scored_ions = $scored_ions . 'B-ions, ';
  }
  if ( defined $settings{'score_yions'} ) {
    $scored_ions = $scored_ions . 'Y-ions, ';
  }
  if ( defined $settings{'score_waterions'} ) {
    $scored_ions = $scored_ions . 'Water-loss ions, ';
  }
  if ( defined $settings{'score_ammoniaions'} ) {
    $scored_ions = $scored_ions . 'ammonia-loss ions, ';
  }

  create_settings($settings_dbh);

  my $settings_sql = $settings_dbh->prepare(
    "INSERT INTO settings 
						(
						      description,
						      cut_residues,
						      protein_sequences,
						      reactive_site,
						      mono_mass_diff,
						      xlinker_mass,
						      decoy,
						      ms2_da,
						      ms1_ppm,
						      finished,
						      isotoptic_shift,	     
						      threshold,
						      charge_match,
						      intensity_match,
						      scored_ions,
						      amber,
						      time,
						      non_specific_digest,
						      no_enzyme_min,
						      no_enzyme_max,		
						      use_previous
						 ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
  );

  $settings_sql->execute(
    $settings{'description'},
    $settings{'enzyme'},
    $settings{'protein_sequences'},
    $settings{'reactive_amino_acid'},
    $settings{'monolink_mass'},
    $settings{'linker_mass'},
    $settings{'decoy'},
    $settings{'accuracy_ms2'},
    $settings{'accuracy_ms1'},
    0,
    $settings{'number_labelled_atoms'},
    $settings{'threshold'},
    $settings{'charge_match'},
    $settings{'intensity_match'},
    $scored_ions,
    $settings{'amber_codon'},
    time,
    $settings{'non_specific_digest'},
    $settings{'no_enzyme_min'},
    $settings{'no_enzyme_max'},
    $settings{'use_previous'}
  );

  my $results_table;
  $results_table = $settings_dbh->{'mysql_insertid'};

  save_modifications(
    $settings{'modifications_fixed'},
    $settings{'modifications_dynamic'},
    $results_table, $settings_dbh
  );

  return $results_table;
}

sub connect_db {

  my ($settings_ref) = @_;
  my %settings = %{$settings_ref};
  my $settings_dbh;

  $settings_dbh =
    DBI->connect( "dbi:mysql:", $settings{'username'},
    $settings{'password'}, { RaiseError => 1, AutoCommit => 1 } );
  $settings_dbh->do("create database if not exists settings");
  $settings_dbh->disconnect;
  $settings_dbh = DBI->connect(
    "dbi:mysql:settings", $settings{'username'},
    $settings{'password'}, { RaiseError => 1, AutoCommit => 1 }
  );
  $settings_dbh->do(
    "CREATE TABLE IF NOT EXISTS status (
						      run_id INTEGER NOT NULL,
						      status TEXT,
						      percent INTEGER,  PRIMARY KEY (run_id)
						) "
  );

  return ($settings_dbh);
}

sub create_results {

  my ($results_dbh) = @_;

  $results_dbh->do(
    "CREATE TABLE IF NOT EXISTS results (
						      rowid INTEGER NOT NULL AUTO_INCREMENT, PRIMARY KEY (rowid),
						      name 		TEXT,
						      MSn_string	TEXT,
						      d2_MSn_string	TEXT,
						      mz		DOUBLE,
						      charge		NUMERIC,
						      fragment		TEXT,    
						      sequence1		TEXT,
						      sequence2		TEXT,
						      sequence1_name	TEXT,
						      sequence2_name	TEXT,
						      score 		REAL,
						      fraction		NUMERIC,
						      scan		NUMERIC,
						      d2_scan		NUMERIC,
						      modification	TEXT,
						      no_of_mods	TEXT,
						      best_x		INTEGER,
						      best_y		INTEGER,
						      unmodified_fragment	TEXT,
						      ppm			FLOAT,
						      top_10			TEXT,
						      d2_top_10			TEXT,
						      matched_abundance		TEXT,
						      d2_matched_abundance	TEXT,
						      total_abundance		TEXT,
						      d2_total_abundance	TEXT,
						      matched_common		TEXT,
     						      matched_crosslink		TEXT,
						      d2_matched_common		TEXT,
						      d2_matched_crosslink	TEXT,
						      monolink_mass		DOUBLE,
						      best_alpha 		REAL,
						      best_beta 		REAL,
						      min_chain_score		TEXT,
						      time			TEXT,
						      precursor_scan		TEXT,
						      FDR			TEXT) "
  );

}

sub create_peptide_table {

  my ($dbh) = @_;

  $dbh->do(
    "CREATE TABLE IF NOT EXISTS peptides ( 
					    rowid INTEGER NOT NULL AUTO_INCREMENT, PRIMARY KEY (rowid),					    results_table NUMERIC,
					    sequence	  TEXT,
					    source	  TEXT,
					    linear_only   NUMERIC,
					    mass 	  DOUBLE,
					    modifications TEXT,
					    monolink 	  DOUBLE,
					    xlink 	  NUMERIC,
					    no_of_mods 	  NUMERIC) "
  );

}

sub connect_db_results {
  my ( $name, $autocommit, $settings_ref ) = @_;
  my %settings = %{$settings_ref};
  if ( !defined $autocommit ) { $autocommit = 1 }
  my $results_dbh;

  $results_dbh =
    DBI->connect( "dbi:mysql:", $settings{'username'},
    $settings{'password'}, { RaiseError => 1, AutoCommit => 1 } );
  $results_dbh->do("create database if not exists results$name");
  $results_dbh->disconnect;
  $results_dbh =
    DBI->connect( "dbi:mysql:results$name", $settings{'username'},
    $settings{'password'},
    { RaiseError => 1, AutoCommit => $autocommit } );

  create_results($results_dbh);

  return ($results_dbh);
}

sub add_peptide {

  my (
    $dbh,           $table,       $sequence,
    $source,        $linear_only, $mass,
    $modifications, $monolink,    $xlink
  ) = @_;

  my $newline = $dbh->prepare(
"INSERT INTO peptides (results_table, sequence, source, linear_only, mass, modifications, monolink, xlink, no_of_mods) VALUES (?,?,?,?,?,?,?,?, 0)"
  );

  $newline->execute( $table, $sequence, $source, $linear_only, $mass,
    $modifications, $monolink, $xlink );
}

sub get_mods {

  my ( $table, $mod_type, $dbh ) = @_;
  my $sql = $dbh->prepare(
    "SELECT  * FROM modifications WHERE run_id = ? AND mod_type = ?");
  $sql->execute( $table, $mod_type );
  return $sql;

}

sub save_modifications {
  my ( $modifications_fixed, $modifications_dynamic, $results_table,
    $settings_dbh )
    = @_;

  $settings_dbh->do(
    "CREATE TABLE IF NOT EXISTS modifications (
						      run_id TEXT,
						      mod_id MEDIUMINT NOT NULL AUTO_INCREMENT, PRIMARY KEY (mod_id),
						      mod_name TEXT,
						      mod_mass DOUBLE,
						      mod_residue TEXT,
						      mod_type TEXT
						) "
  );

  my $settings_sql = $settings_dbh->prepare(
    "INSERT INTO modifications 
						    (
							  run_id,
							  mod_name,
							  mod_mass,
							  mod_residue,
							  mod_type
						    ) VALUES (?,?,?,?,?)"
  );
  my @fixed_mods = split( ',', $modifications_fixed );
  foreach my $mod (@fixed_mods) {
    my @values = split( ':', $mod );
    $settings_sql->execute( $results_table, $values[0], $values[1],
      $values[2], 'fixed' );
  }

  my @dynamic_mods = split( ',', $modifications_dynamic );
  foreach my $mod (@dynamic_mods) {
    my @values = split( ':', $mod );
    $settings_sql->execute( $results_table, $values[0], $values[1],
      $values[2], 'dynamic' );
  }

}



1;