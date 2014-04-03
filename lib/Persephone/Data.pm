use strict;

package Persephone::Data;

use base 'Exporter';
use lib 'lib';
use Persephone::Config;

our @EXPORT = (
  'import_fasta',                 'connect_db',
  'protein_residuemass',          'calculate_peptide_masses',
  'calculate_crosslink_peptides', 'generate_monolink_peptides',
  'generate_modified_peptides', 'import_mgf',
  'loaddoubletlist_db'
);

sub no_of_threads {
  return 16;
}


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

  print "Run $results_table: Calculating modified peptides ... \n";

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

  print "Run $results_table: Calculating monolink ... \n";

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
  my $missed_clevages    = $settings{'max_missed_cleavages'},
    my $cut_residues     = $settings{'enzyme_cut_residues'},

    my $stop_duplicates = 'AND p1.rowid >= p2.rowid';

  if ( length($reactive_site) < 2 ) {
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
    print "Run $results_table: Crosslinking peptides  ... ",
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

sub get_mods {

  my ( $table, $mod_type, $dbh ) = @_;
  my $sql = $dbh->prepare(
    "SELECT  * FROM modifications WHERE run_id = ? AND mod_type = ?");
  $sql->execute( $table, $mod_type );
  return $sql;

}



sub import_mgf    #Enters the uploaded MGF into a SQLite database
{

    my ($results_id, $settings_ref) = @_;

    my %settings = %{$settings_ref};
    my $filename = $settings{'directory'}."/".$settings{'input'};
    my $fraction = 1;
    my $dbh = connect_db_results($results_id, 0, \%settings);
    my %line;
    my $MSn_count = 0;
    my $dataset   = 0;
    my $MSn_string;

    $line{'fraction'} = $fraction;
      
    open FILE, "<$filename", or die "Error: Cannot open MGF."; 
  
    print "Importing $filename ... \n";
  
    my $masslist = $dbh->prepare("DROP TABLE IF EXISTS msdata;");
    $masslist->execute();
        $dbh->do(
	"CREATE TABLE msdata (
		scan_num 	NUMERIC,
		fraction 	NUMERIC, 
		title	 	TEXT, 
		charge 		NUMERIC,
		mz		DOUBLE,
		monoisotopic_mw DOUBLE,
		abundance 	FLOAT,
		MSn_string	TEXT,
		msorder 	NUMERIC,
		precursor_scan	NUMERIC) "
        );
    
    while (<FILE>) {
        if ($_ =~ "^BEGIN IONS") { $dataset = $dataset + 1; $line{'abundance'} = 0 }
        elsif ($_ =~ "^PEPMASS") {
             my $mystring = $_;	
            if ($mystring =~ m/=(.*?) /)      { $line{'mz'}        = $1;}
	    elsif ($mystring =~ m/=(.*?)[\r\n]/) { $line{'mz'}     = $1;}
            if ($mystring =~ m/ (.*?)[\r\n]/) { $line{'abundance'} = $1 ;}
        } elsif ($_ =~ "^SCANS") {
            my $mystring = $_;
            if ($mystring =~ m/=(.*?)[\r\n]/) { $line{'scan_num'} = $1 }
        } elsif ($_ =~ "^CHARGE") {
            my $mystring = $_;
            if ($mystring =~ m/=(.*?)\+/) { $line{'charge'} = $1 }
        } elsif ($_ =~ "^TITLE") {
            my $mystring = $_;
            if ($mystring =~ m/=(.*?)[\r\n]/) { $line{'title'} = $1 }
        }

        elsif ($_ =~ "^.[0-9]") {
            my $MSn_row = $_;
            $MSn_count  = $MSn_count + 1;
            $MSn_string = $MSn_string . $MSn_row;
            my @MSn_split = split(/ /, $MSn_row);
            my ($ms2_mz, $ms2_abundance) = @MSn_split;

            #       $scan_data->execute($line{'scan_num'},$ms2_mz, $ms2_abundance);
        }

        elsif ($_ =~ "^END IONS") {
            $line{'monoisoptic_mw'} = $line{'mz'} * $line{'charge'} - ($line{'charge'} * 1.00728);
            my $newline = $dbh->prepare(
"INSERT INTO msdata (scan_num, fraction, title, charge, mz, abundance, monoisotopic_mw, MSn_string, msorder) VALUES (? , ?, ?, ?, ?, ?, ?,?, 2)"
            );
                $newline->execute($line{'scan_num'}, $line{'fraction'}, $line{'title'}, $line{'charge'}, $line{'mz'},
                                  $line{'abundance'}, $line{'monoisoptic_mw'}, $MSn_string);

            # 	 warn "Scan imported \n";

            $line{'scan_num'} = $line{'monoisoptic_mw'} = $line{'abundance'} = $MSn_string = '';
            $MSn_count = 0;
        }
    }

    $dbh->commit;
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

  print "No Cut:", $nocut_residues, "\n";
  my @digest;
  my @digest_not_for_crosslinking;
  if ( $n_or_c eq 'C' ) {

    print "Protease type:", $n_or_c, "\n";
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

  my ( $results_table, $protein_residuemass_ref, $settings_ref ) = @_;
  my %protein_residuemass = %{$protein_residuemass_ref};
  my %settings            = %{$settings_ref};
  my $peptide_mass        = 0;
  my $terminalmass        = 1.0078250 * 2 + 15.9949146 * 1;
  my $count               = 0;

  print "Run $results_table: Calulating masses...  \n";

  my $results_dbh = connect_db_results( $results_table, 1, \%settings );
  my $settings_dbh = connect_db( \%settings );
  my $sequencelist = $results_dbh->prepare(
    "select distinct source from peptides where xlink = 0");
  $sequencelist->execute;
  my $rows = $sequencelist->rows;

  my $threads = 0;
  $threads = no_of_threads;
  my $pm = Parallel::ForkManager->new($threads);

  while ( my $source = $sequencelist->fetchrow_hashref ) {
    print "Run $results_table: Calculating masses ... ",
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

sub import_fasta {
  my ($settings_ref) = @_;
  my %settings = %{$settings_ref};

  my $filename = $settings{directory} . "/" . $settings{fasta};

  open( my $fh, $filename ) or die "Error: Could not open Fasta file\n";
  $settings{'protein_sequences'} = do { local ($/); <$fh> };

  $settings{'protein_sequences'} =~ s/\r//g;
  my @sequence_names = $settings{'protein_sequences'} =~ m/^>.*$/mg;
  $settings{'protein_sequences'} =~ s/^>.*$/>/mg;
  $settings{'protein_sequences'} =~ s/\n//g;
  $settings{'protein_sequences'} =~ s/^.//;
  $settings{'protein_sequences'} =~ s/ //g;
  my @sequences = split '>', $settings{'protein_sequences'};

  my ($settings_dbh) = connect_db( \%settings );
  my $results_table = save_settings( $settings_dbh, \%settings );
  my ($results_dbh) =
    connect_db_results( $results_table, 0, \%settings );

  print "Results ID:", $results_table, "\n";

  create_peptide_table($results_dbh);

  my $count = 0;

  foreach my $sequence (@sequences) {

    my $sequence_fragments_ref;
    my $sequence_fragments_linear_only_ref;
    my @sequence_fragments;
    my @sequence_fragments_linear_only;

    if ( !defined $settings{'non_specific_digest'}
      || $settings{'non_specific_digest'} == 0 )
    {
      ( $sequence_fragments_ref, $sequence_fragments_linear_only_ref ) =
	digest_proteins(
	$settings{'max_missed_cleavages'},
	$sequence,
	$settings{'enzyme_cut_residues'},
	$settings{'enzyme_no_cut_residues'},
	$settings{'enzyme_cut_n_or_c_term'}
	);
      @sequence_fragments = @{$sequence_fragments_ref};
      @sequence_fragments_linear_only =
	@{$sequence_fragments_linear_only_ref};
    }

#  	  Alternative digestion methods found in Hekate/Crosslinker yet to be added to Persephone
# 	    elsif (defined $settings{'non_specific_digest'} && $settings{'non_specific_digest'}  == 2) {
#  	      ($sequence_fragments_ref) = amino_peptidase_digest($no_enzyme_min, $no_enzyme_max, $reactive_site, $sequence);
#  	      @sequence_fragments             = @{$sequence_fragments_ref};
#  	  } else {
#  	      ($sequence_fragments_ref) = no_enzyme_digest_proteins($no_enzyme_min, $no_enzyme_max, $reactive_site, $sequence);
#  	      @sequence_fragments             = @{$sequence_fragments_ref};
#  	  }

    print
"Run $results_table: Sequence $count = $sequence_names[$count] \n";
    print "Run $results_table: Digested peptides:",
      scalar(@sequence_fragments), " \n";

    foreach my $fragment (@sequence_fragments)

    {
      add_peptide( $results_dbh, $results_table, $fragment, $count, 0,
	0, '', 0, 0 );
    }

    foreach my $fragment (@sequence_fragments_linear_only) {
      add_peptide( $results_dbh, $results_table, $fragment, $count, 1,
	0, '', 0, 0 );
    }

    $count++;
  }

  $results_dbh->commit;
  $results_dbh->disconnect;

  return $results_table;

}

1;