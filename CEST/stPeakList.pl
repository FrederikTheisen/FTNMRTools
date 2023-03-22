#! /usr/bin/perl

# Creates a peak list for seriesTab given a 
# nmrPipe spektrum and a Sparky peak list
#
# usage: ./stPeakList [pipe spectrum] [Sparky peak list]
#
# Kaare Teilum 061025

$pipeFile = $ARGV[0];
$sparkyList = $ARGV[1];

open HDR, "showhdr $pipeFile |" or die "can't run: showhdr $pipeFile";
while (<HDR>){
        if (/OBS MHz/){
                /([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)/;
                $frqX=$1;
                $frqY=$2;
        };
        if (/DATA SIZE/){
                /([0-9]+)\s+([0-9]+)/;
                $sizeX=$1;
                $sizeY=$2;
        };
        if (/ORIG Hz/){
                /(-*[0-9]+\.[0-9]+)\s+(-*[0-9]+\.[0-9]+)/;
                $origX=$1;
                $origY=$2;
        };
        if (/SW Hz/){
                /([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)/;
                $swX=$1;
                $swY=$2;
        };
};
close(HDR);

$stepX=$swX/$frqX/$sizeX;
$highX=($origX+$swX)/$frqX;
$stepY=$swY/$frqY/$sizeY;
$highY=($origY+$swY)/$frqY;
$i=1;
print "VARS   INDEX X_AXIS Y_AXIS X_PPM Y_PPM VOL ASS\n";
print "FORMAT %5d %9.3f %9.3f %8.3f %8.3f %+e %s\n\n";

open IN, "$sparkyList" or die "Cannot open $sparkyList for read";
        while (<IN>){
                @process = split (/\s+/, $_);
                if ($process[0] eq ""){splice (@process, 0, 1)};
                if ($process[0] ne "Assignment" && $process[0] ne ""){
                        $ptsX=($highX-$process[1])/$stepX;
                        $ptsY=($highY-$process[2])/$stepY;
                        printf "%5d %9.3f %9.3f %8.3f %8.3f %+e %s\n",$i,$ptsX,$ptsY,$process[1],$process[2],700000,$process[0];
                        $i++;
                };
        };
close (IN);

