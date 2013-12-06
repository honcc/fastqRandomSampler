#!/usr/bin/perl/ -w
$|++; #---turn on the auto flush for the progress bar
use strict;
use Math::Random; #---For generating normal distribution random numbers
use File::Path; #---for removing tmp GNUPLOT dat files
use List::Util 'shuffle';#--- for shuffling an array
use POSIX qw(log10);
use Time::HiRes qw( time );
######################################################################################################################################################
#
#	Description
#		This is a perl script to random sample sequences from a fastq file according to user defined number of sequences. Single ane pairend are supporte;
#
#	Input
#		--fastq1Path=		the path of the first mate fastq file if pairend, or the fastq for singleend data;
#		--randomReadNum=	the number of reads to be randomly sampled (in unit of millions);
#		--outDir			the output directory;
#			
#	Output
#
#	Usage
#		perl fastqRandomSampler_v0.1.pl --fastq1Path=/Volumes/BCPNGS/s1.fastq_1 --fastq2Path=/Volumes/BCPNGS/s1.fastq_2 --randomReadNumMil=1 --outDir=/Users/chung/Desktop/NGS/March2011Pilot/ --filterLen=0 --headerTag=HWI
#
#	Assumption 
#
#
#	Version history
#
#		v0.1
#
#		v0.2
#			-it takes and output gz instead of plain txt 
#			-fastq2Path, headerTag and filter length option removed;
#
#####################################################################################################################################################

#==========================================================Main body starts==========================================================================#
#---1 read the parameters
my ($fastq1Path, $randomReadNumMil, $outDir) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

#---generate random index chunks by chunks (use chunk since splice large array is very slow)
my ($randomPickedReadNumAry_ref, $fastq1TotalReadNum, $intervalSize) = generateRandomSamplingIndex($fastq1Path, $randomReadNumMil);

#---pick the preselected seq on the fly
pickPreSelectedSeqOnTheFly($randomPickedReadNumAry_ref, $fastq1TotalReadNum, $intervalSize, $fastq1Path);

printCMDLogOrFinishMessage("finishMessage");

########################################################################## readParameters
sub readParameters {
	
	my $randomReadNumMil = 1; #---default value is one million
	my $fastq1Path = "";
	my $outDir = "./";
	
	foreach my $param (@ARGV) {
		if ($param =~ m/--fastq1Path=/) {$fastq1Path = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--randomReadNumMil=/) {$randomReadNumMil = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);}
	}
	
	#---check fastq pair
	die "Cant open $fastq1Path\n" if not (-s "$fastq1Path");#---single end
	
	chop $outDir if ($outDir =~ m/\/$/); #--remove the last slash
	system ("mkdir -p -m 777 $outDir");
	
	return ($fastq1Path, $randomReadNumMil, $outDir);
}
########################################################################## generateRandomSamplingIndex
sub generateRandomSamplingIndex {

	my ($fastq1Path, $randomReadNumMil) = @_;

	#---estimate the number of lines in the file
    print "Estimating the number of lines.\n";
   	my $fastq1TotalLineNum = `gzip -d -c $fastq1Path | wc -l`;
   	my $fastq1TotalReadNum = $fastq1TotalLineNum/4;
   	print "Estimated to have ".$fastq1TotalReadNum." reads.\n";

	my $intervalSize = int ($fastq1TotalReadNum/100); #---define as 
	$intervalSize = 1000000 if ($intervalSize > 1000000);

	#--convert to millions
	my $randomReadNum = $randomReadNumMil*1000000;
	
	#---quit if randomReadNum > 90% of total Readnum
	die "randomReadNum is more than 90% of the total read num in the fastq fie. Quitting\n" if ($randomReadNum > ($fastq1TotalReadNum*0.9)); #---$randomReadNum is very close to $fastq1TotalReadNum
	
	#---generate a array contains the random read order to pick
	my $chunkSize = 1000; #---0.1 M as a chunk
	my $numberOfChunks = sprintf "%.0f", $fastq1TotalReadNum/$chunkSize;
	$numberOfChunks = 1 if ($numberOfChunks < 1);
	my $samplePerChunk = sprintf "%.0f", ($randomReadNum/$numberOfChunks);
	$samplePerChunk = 1 if ($samplePerChunk < 1);

	#---get the differences of rounded sample number and actual sample number
	#---will be redistributed later to randomly chosen chunk;
	my $roundDifference = $randomReadNum - ($samplePerChunk*$numberOfChunks);
	my @chunkIDAry = (1..$numberOfChunks);
	my %randChunkRedistributeHsh;
	
	#---redistribute the rounded difference into randomly selected chunk ID
	for my $i (0..$roundDifference) {
		my $randomIndex = rand $#chunkIDAry;
		my $randomChunkID = $chunkIDAry[$randomIndex];
		$randChunkRedistributeHsh{$randomChunkID}++;
		splice @chunkIDAry, $randomIndex, 1;
	}
	
	print "Generating an array of $randomReadNum indexes to be sampled\n";
	my @chunkRngAry;

	#---generation the chunk ranges
	for my $i (0..$numberOfChunks) {#----no problem even $chunkSize > $fastq1TotalReadNum  
		my $chunkStart = 1 + $i*$chunkSize;
		my $chunkEnd = ($i+1)*$chunkSize;
		last if ($chunkStart > $fastq1TotalReadNum);
		$chunkEnd = $fastq1TotalReadNum if ($chunkEnd > $fastq1TotalReadNum);
		push @chunkRngAry, $chunkStart.":".$chunkEnd."\n";
	}
	
	#---go through all chunks and sample the index
	my @randomPickedReadNumAry;
	my $chunkProc = 0;
	my $storedNum = 0;
	my $chunkID = 0;
	foreach my $chunkRng (@chunkRngAry) {
		#print "Generating Index for $chunkProc of $numberOfChunks chunks.\n";
		$chunkID++;
		my $roundDiffRedistbtSamplePerChunk = $samplePerChunk;
		$roundDiffRedistbtSamplePerChunk += 1 if (exists $randChunkRedistributeHsh{$chunkID});
		my @chunkRngSplt = split /\:/, $chunkRng;
		my @chunkNumAry = ($chunkRngSplt[0]..$chunkRngSplt[1]);
		my @tmpChunkRandNumAry;
		foreach my $i (1..$roundDiffRedistbtSamplePerChunk) {
			my $randomIndex = rand $#chunkNumAry;
			my $randomNum = $chunkNumAry[$randomIndex];
			push (@tmpChunkRandNumAry, $randomNum);
			splice @chunkNumAry, $randomIndex, 1;
			$storedNum++;
			last if $storedNum == $randomReadNum;
		}
		@tmpChunkRandNumAry = sort {$a <=> $b} @tmpChunkRandNumAry;
		push (@randomPickedReadNumAry, @tmpChunkRandNumAry);
		$chunkProc++; 
	}
	
	print $storedNum." random index stored.\n";
	
	return (\@randomPickedReadNumAry, $fastq1TotalReadNum, $intervalSize);
	
}
########################################################################## pickPreSelectedSeqOnTheFly
sub pickPreSelectedSeqOnTheFly {
	
	my ($randomPickedReadNumAry_ref, $fastq1TotalReadNum, $intervalSize, $fastq1Path) = @_;

	my @randomPickedReadNumAry = @{$randomPickedReadNumAry_ref};
	
	my @fastq1PathSplt = split /\//, $fastq1Path;
	my $fastq1OutName = $fastq1PathSplt[-1];
	$fastq1OutName =~ s/\.(\w+)$/\.resam\.$randomReadNumMil.M\.$1/;
	open (INFQ1, "pigz -c -d $fastq1Path |");
	my $fastq1OutPath = "$outDir/$fastq1OutName";
	open (OUTFQ1, "| pigz -c - >$fastq1OutPath");

	#---define the start time and counters
	my $intervalStart = time();
	my $headerProc = my $progCount = 0;

	print "Sampling process started. End time will be estimated after processing the first $intervalSize reads.\n";

	my (%sampledLenCount1Hsh, %sampledLenCount2Hsh);
	my $randomReadNum = $randomReadNumMil*1000000;
	my $headerCount = my $sampledSeqCount = 0;
	my $totalLineCount = 3;
	while (my $infq1Line = <INFQ1>) {
		$totalLineCount++;
		next if ($totalLineCount % 4 != 0);#---skip all non-header line;
		
		$headerProc++; $progCount++; $headerCount++;
		if ($progCount >= $intervalSize) {
			($progCount, $intervalStart) = reportProgress($progCount, $headerProc, $intervalSize, $fastq1TotalReadNum, $intervalStart);
		}
		
		last if ($sampledSeqCount == $randomReadNum); #---stop if the array is empty
		
		if ($randomPickedReadNumAry[0] == $headerCount) {#---this seq is prepicked
			$sampledSeqCount++;
			splice @randomPickedReadNumAry, 0, 1;
			print OUTFQ1 $infq1Line;
			for my $i (1..3) {
				$infq1Line = <INFQ1>;
				print OUTFQ1 $infq1Line;
				$totalLineCount++;
			}
		}
	}#---end of while (my $infq1Line = <INFQ1>)
	
	print "Totally $sampledSeqCount reads were randomly sampled.\n";

}
########################################################################## reportProgress
sub reportProgress {

	my $progCount = $_[0];
	my $lineProc = $_[1];
	my $intervalSize = $_[2];
	my $fileTotalLineNum = $_[3];
	my $intervalStart = $_[4];

	$progCount=0;
	my $intervalEnd = time();
	my $timeElapsed = $intervalEnd - $intervalStart;
	$timeElapsed = sprintf ("%.2f", $timeElapsed);
	my $estimatedEnd = (($fileTotalLineNum - $lineProc)*$timeElapsed)/$intervalSize;
	$estimatedEnd = sprintf ("%.2f", $estimatedEnd/60);
	print "$lineProc headers processed. Last $intervalSize header:".$timeElapsed." sec. Estimated end: ".$estimatedEnd." mins.                  \r";
	$intervalStart = time();
	
	return ($progCount, $intervalStart);
		
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>./$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {

		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
