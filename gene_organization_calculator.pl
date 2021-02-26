use strict;

################This script takes TAIR10 annotation, and determine if each gene has convergent or tandem downstream gene and calculate calculate the distance to downstream gene.#################

#Use two annotation files (one sorted by start position, and the other sorted by end position)
my $file="TAIR10_all_nuclear_PCG_TEG_PG32431.txt";
my(@chr, @start, @end, @ID, @annotation, @direction, @downstream_ID, @downstream_kind, @distance_to_next_start, @distance_to_next_end, @covered);
open(IN, $file);
while(my $gene = <IN>){
	chomp($gene);
	my ($chr, $start, $end, $ID, $annotation, $direction) = split(/\t+/, $gene);
	push @chr, $chr;
	push @start, $start;
	push @end, $end;
	push @ID, $ID;
	push @annotation, $annotation;
	push @direction, $direction;
}
close(IN);

my $file2="TAIR10_all_nuclear_PCG_TEG_PG32431_endsort.txt";
my(@chr2, @start2, @end2, @ID2, @annotation2, @direction2, @downstream_ID2, @downstream_kind2, @distance_to_next_start2, @distance_to_next_end2, @covered2);
open(IN, $file2);
while(my $gene = <IN>){
	chomp($gene);
	my ($chr, $start, $end, $ID, $annotation, $direction) = split(/\t+/, $gene);
	push @chr2, $chr;
	push @start2, $start;
	push @end2, $end;
	push @ID2, $ID;
	push @annotation2, $annotation;
	push @direction2, $direction;
}
close(IN);

#Determine if the downstream gene is convergent or tandem, and calculate the distance to the downstream gene.

for (my $i = 0; $i <= 32430; $i++){
	my $j = 0;
	until ($chr[$j] eq $chr[$i]){$j++;}
	while ($chr[$j] eq $chr[$i]){
		if ($start[$i] >= $start[$j] && $end[$i] <= $end[$j] && $i != $j){
			$covered[$i] = "covered";
			$covered[$j] = "covering";
			$downstream_ID[$i] = $ID[$j];
			$downstream_ID[$j] = $ID[$i];
			if ($direction[$j] ne $direction[$i]){
				$downstream_kind[$i] = $downstream_kind[$j] = "convergent";
			}
			elsif ($direction[$j] eq $direction[$i]){
				$downstream_kind[$i] = $downstream_kind[$j] = "tandem";
			}
			$distance_to_next_start[$i] = $start[$j] - $end[$i];
			$distance_to_next_end[$i] = $end[$j] - $end[$i];
			$distance_to_next_start[$j] = $start[$i] - $end[$j];
			$distance_to_next_end[$j] = $end[$i] - $end[$j];
			last;
		}
		else {$j++;}
	}
	#unless ($covered[$i] eq "true"){$covered[$i] = "false";}
}

for (my $i = 0; $i <= 32430; $i++){
	my $j = 0;
	until ($chr2[$j] eq $chr2[$i]){$j++;}
	while ($chr2[$j] eq $chr2[$i]){
		if ($start2[$i] >= $start2[$j] && $end2[$i] <= $end2[$j] && $i != $j){
			$covered2[$i] = "covered";
			$covered2[$j] = "covering";
			$downstream_ID2[$i] = $ID2[$j];
			$downstream_ID2[$j] = $ID2[$i];
			if ($direction2[$j] ne $direction2[$i]){
				$downstream_kind2[$i] = $downstream_kind2[$j] = "convergent";
			}
			elsif ($direction2[$j] eq $direction2[$i]){
				$downstream_kind2[$i] = $downstream_kind2[$j] = "tandem";
			}
			$distance_to_next_start2[$i] = $start2[$j] - $end2[$i];
			$distance_to_next_end2[$i] = $end2[$j] - $end2[$i];
			$distance_to_next_start2[$j] = $start2[$j] - $end2[$i];
			$distance_to_next_end2[$j] = $start2[$j] - $start2[$i];
			last;
		}
		else {$j++;}
	}
	#unless ($covered2[$i] eq "true"){$covered2[$i] = "false";}
}

for (my $i = 0; $i <= 32430; $i++){
	if ($direction[$i] eq "+" && $chr[$i] eq $chr[$i+1] && $covered[$i] eq ""){
		$covered[$i] = "-";
		$downstream_ID[$i] = $ID[$i+1];
		if ($direction[$i+1] eq "-"){
			$downstream_kind[$i] = "convergent";
		}
		elsif ($direction[$i+1] eq "+"){
			$downstream_kind[$i] = "tandem";
		}
		$distance_to_next_start[$i] = $start[$i+1] - $end[$i];
		$distance_to_next_end[$i] = $end[$i+1] - $end[$i];
	}
	elsif ($direction[$i] eq "+" && $covered[$i] eq "false"){$downstream_ID[$i] = "na";$downstream_kind[$i] = "na";$distance_to_next_start[$i] = "na";$distance_to_next_end[$i] = "na";}
}
for (my $i = 0; $i <= 32430; $i++){
	if ($direction2[$i] eq "-" && $chr2[$i] eq $chr2[$i-1] && $covered2[$i] eq ""){
		$covered2[$i] = "-";
		$downstream_ID2[$i] = $ID2[$i-1];
		if ($direction2[$i-1] eq "+"){
			$downstream_kind2[$i] = "convergent";
		}
		elsif ($direction2[$i-1] eq "-"){
			$downstream_kind2[$i] = "tandem";
		}
		$distance_to_next_start2[$i] = $start2[$i] - $end2[$i-1];
		$distance_to_next_end2[$i] = $start2[$i] - $start2[$i-1];
	}
	elsif ($direction[$i] eq "-" && $covered2[$i] eq "false"){$downstream_ID2[$i] = "na";$downstream_kind2[$i] = "na";$distance_to_next_start2[$i] = "na";$distance_to_next_end2[$i] = "na";}
}	


my $output_file="downstream2.txt";
open(OUT, ">$output_file");
print OUT join("\t", "chr","start","end","ID","annotation","direction","downstream_ID","downstream_kind","distance_to_next_start","distance_to_next_end","covered?"), "\n";

for (my $i =0; $i <= 33322; $i++){
	if ($direction[$i] eq "+"){
	print OUT join("\t", $chr[$i],$start[$i],$end[$i],$ID[$i],$annotation[$i],$direction[$i],$downstream_ID[$i],$downstream_kind[$i],$distance_to_next_start[$i],$distance_to_next_end[$i],$covered[$i]), "\n";
	}
}
for (my $i =0; $i <= 33322; $i++){
	if ($direction2[$i] eq "-"){
	print OUT join("\t", $chr2[$i],$start2[$i],$end2[$i],$ID2[$i],$annotation2[$i],$direction2[$i],$downstream_ID2[$i],$downstream_kind2[$i],$distance_to_next_start2[$i],$distance_to_next_end2[$i],$covered2[$i]), "\n";
	}
}
	
close(OUT);