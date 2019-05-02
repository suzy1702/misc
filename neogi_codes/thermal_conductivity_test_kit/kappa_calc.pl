#!/usr/bin/perl

if($#ARGV<3){
   print STDOUT "kappa_calc.pl nruns numfirstrun numcorr temp thickness\n";
   exit;
}
$dirroot="heat";
$nruns=$ARGV[0];
$Tk0=$ARGV[1]; 
$numcorr=$ARGV[2];
$T=$ARGV[3];
$thickness=$ARGV[4];
$corrtime=$numcorr*0.00025;
$eintname="einstein_";
$log=".log";
$diffusionname="diffusion";
$dat=".dat";

for($Tk=$Tk0;$Tk<($nruns+$Tk0);$Tk++){
	chdir "$dirroot.$Tk";
   	print "$dirroot.$Tk \n";
	system "/projects/sane3962/bin/einst2D.x $numcorr $T HEATFLUX 100 > $eintname$corrtime$log \n";
	system "tail $eintname$corrtime$log \n";
	system "mv diffusion.dat $diffusionname$corrtime$dat \n";
#	system "xmgrace -block $diffusionname$corrtime$dat -bxy 1:2\n";

	chdir "../";
					
}

