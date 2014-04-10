use strict;

package Persephone::Import;

use base 'Exporter';
use lib 'lib';
use Persephone::Config;
use Persephone::Data;
use Persephone::Proteins;

our @EXPORT = ('import_fasta', 'import_mgf', 'import_settings', 'import_peptides');

sub import_settings {

my ($settings_ref) = @_;
  my %settings = %{$settings_ref};

  $settings{'protein_sequences'} = $settings{directory} . "/" . $settings{fasta};

  my ($settings_dbh) = connect_db( \%settings );
  my $results_table = save_settings( $settings_dbh, \%settings );
  my ($results_dbh) =
    connect_db_results( $results_table, 0, \%settings );

  print "Results ID:", $results_table, "\n";


  $results_dbh->commit;
  $results_dbh->disconnect;

  return $results_table;


}
sub import_fasta {
     my ($results_table, $settings_ref) = @_;
     

   
     my %settings = %{$settings_ref};
     create_protein_database($results_table, \%settings);
     my ($results_dbh) = connect_db_results( $results_table, 0, \%settings );
      
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

    create_peptide_table($results_dbh);

    
    
    
    my $count = 0;
  
  foreach my $sequence (@sequences) { 
    print "Imported protein $count: $sequence_names[$count]\n";
    add_proteins ($results_dbh, $sequence_names[$count], $sequence, $results_table);
    $count = $count +1;
  }
  
  
#   import decoy?

  if ($settings{'decoy'} == '1') {
       my $count = 0;
    foreach my $sequence (@sequences) { 
    print "Imported decoy $count: $sequence_names[$count]\n";
      my $decoy = reverse $sequence;
        $decoy         =~ tr/KR/RK/;
    add_proteins ($results_dbh, ">decoy" .$sequence_names[$count], $decoy, $results_table);
    $count = $count +1;
  }
  
  }
  
  $results_dbh->commit;
  $results_dbh->disconnect;
}

sub import_peptides {
     my ($results_table, $settings_ref) = @_;
     

   
     my %settings = %{$settings_ref};
     my ($results_dbh) = connect_db_results( $results_table, 0, \%settings );
       
  my $count = 0;  
  
  my $sequences = $results_dbh->prepare("SELECT * from proteins where run_id = ?");
  $sequences->execute($results_table);
  
  while ( my $sequence = $sequences->fetchrow_hashref ) {

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
	$sequence->{'sequence'},
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

    print "Digesting sequence $count = $sequence->{'name'} \n";
    print "Peptides:",
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
     
     
     
     
     return;

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

1;