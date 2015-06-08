#!/usr/bin/perl
use strict;
use warnings;

my %energy_alphabet=(); 
my (@table,@energy_table,@x_axis,@y_axis)=();
my ($x_axe,$y_axe)=();

while(my $line=<>) {
	chomp($line);
	if (my ($letter,$energy) = $line =~ m{"([A-Z]).*\s/\*\s"([0-9.-]+)}) {
		$energy_alphabet{$letter} = $energy;
	} elsif ( my ($x_axis) = $line =~ m{x-axis:(.*)\s\*/}) {
		$x_axe = $x_axis;
	} elsif (my ($y_axis) = $line =~ m{y-axis:(.*)\s\*/}) {
		$y_axe = $y_axis;	
	} elsif (my ($table_line) = $line =~ m{"([A-Z]+)"}) {
		push(@table, $table_line);
	}
}

foreach my $table_line(@table) {
	my (@table_letters) = (split '',$table_line);
	foreach my $table_letter(@table_letters) {
		my $energy_value = letter_to_energy($table_letter,%energy_alphabet);
		push(@energy_table,$energy_value);
	}
	push(@energy_table,"\n");
}

#my @xcoordinates = (split ' ',$x_axe);
#my @ycoordinates = (split ' ',$y_axe);

open(my $fh, '>', 'pes_profile.txt');

	print $fh "$x_axe \n";
	print $fh "$y_axe \n";
	foreach my $energyValue(@energy_table){
		print $fh "$energyValue ";
	}

close $fh;


# subroutine definition
sub letter_to_energy {
	my ($value, %hash) = @_;	# subroutine argument
  	if(exists $hash{$value}) {
		return $hash{$value};
  } else {
		my $fake = 0 ;
      return $fake;
    }
}


