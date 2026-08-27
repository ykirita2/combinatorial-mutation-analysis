#!/usr/bin/env perl
#
# Copyright (c) 2024 GeneBay, Inc. All rights reserved.
#
use Bio::DB::Sam;
use FileHandle;

my $version = "1.0";
my $fFAS = shift or usage();
my $fBED = shift or usage();
my $fBAM = shift or usage();

%aacode = (
  TTT => "F", TTC => "F", TTA => "L", TTG => "L",
  TCT => "S", TCC => "S", TCA => "S", TCG => "S",
  TAT => "Y", TAC => "Y", TAA => "*", TAG => "*",
  TGT => "C", TGC => "C", TGA => "*", TGG => "W",
  CTT => "L", CTC => "L", CTA => "L", CTG => "L",
  CCT => "P", CCC => "P", CCA => "P", CCG => "P",
  CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
  CGT => "R", CGC => "R", CGA => "R", CGG => "R",
  ATT => "I", ATC => "I", ATA => "I", ATG => "M",
  ACT => "T", ACC => "T", ACA => "T", ACG => "T",
  AAT => "N", AAC => "N", AAA => "K", AAG => "K",
  AGT => "S", AGC => "S", AGA => "R", AGG => "R",
  GTT => "V", GTC => "V", GTA => "V", GTG => "V",
  GCT => "A", GCC => "A", GCA => "A", GCG => "A",
  GAT => "D", GAC => "D", GAA => "E", GAG => "E",
  GGT => "G", GGC => "G", GGA => "G", GGG => "G",
);

my $sam = Bio::DB::Sam->new(-bam => $fBAM, -fasta => $fFAS);
my @BED = readBED($fBED, $sam);
my $hrR={};
my $full_len=0;

foreach my $bed (@BED) {
  # retrieve read alignments that overlap with the specified region
  my @alignments = $sam->get_features_by_location(-seq_id =>$bed->{chr}, -start=>$bed->{stt}+1, -end=>$bed->{end});
  my @seq=@{$bed->{ar_seq}}; # allowed sequence patterns to be counted at the region
  $full_len += $bed->{end} - $bed->{stt};

  for my $a (@alignments) {
    my $seqid  = $a->seq_id;
    my $start  = $a->start;
    my $end    = $a->end;
    my $strand = $a->strand;
    my $cigar  = $a->cigar_str;
    # where does the alignment start in the query sequence
    my $query_start = $a->query->start;
    my $query_end   = $a->query->end;
    my $ref_dna   = $a->dna;        # reference sequence bases
    my $query_dna = $a->query->dna; # query sequence bases
    my @scores    = $a->qscore;     # per-base quality scores
    my $match_qual= $a->qual;       # quality of the match

    next unless $a->start <= 1; # skip if the alignment is not fully aligned against the reference
    next unless $a->end >= $sam->length($bed->{chr}); # skip if the alignment is not fully aligned against the reference
    next if $a->flag & 2304; # only count primary alignments

    # extract subsequence of the read at the region
    my $subseq = extract_subsequences($bed->{stt}+1, $bed->{end}, $a->start, $a->end, $a->strand, $a->cigar_str, $a->query->start, $a->query->end, $a->query->dna, $a->qname);
    if ($subseq) {
      my $flag=0;
      # check to see if the extracted subsequence does match to the candidate sequence
      foreach my $seq (@seq) {
        $flag++ if uc($subseq) eq uc($seq);
      }
      if ($flag) {
        $hrR->{$a->qname} .= $subseq;
      }
    }
  }
}

print "\t", "Fully aligned reads", "\t", scalar(keys %$hrR), "\n";

my $hrS={};
my $nread=0;
for my $read_id (keys %$hrR) {
  my $seq = $hrR->{$read_id};
  if (length($seq) == $full_len) {
    $hrS->{$seq}++;
    $nread++;
  }
}

my @pattern=();

foreach my $bed (@BED) {
  my @seq=@{$bed->{ar_seq}};
  my @output=();

  if (@pattern) {
    foreach my $pat (@pattern) {
      foreach my $seq (@seq) {
        push @output, $pat . $seq;
      }
    }
    @pattern=@output;
  } else {
    @pattern=@seq;
  }
}

print "\t", "Assigned reads", "\t", $nread, "\n";
foreach my $seq (@pattern) {
  print join("\t", $seq, translate($seq), $hrS->{$seq}||0), "\n";
}

sub extract_subsequences {
    my $from=shift;
    my $to=shift;
    my $rst=shift;
    my $red=shift;
    my $rori=shift;
    my $cigar=shift;
    my $qst=shift;
    my $qed=shift;
    my $qseq=shift;
    my $qname=shift;

    while ($cigar=~m/(\d+)([A-Za-z])/g) {
      my $n=$1;
      my $op=$2;
      if ($op eq 'M') {
        my $ss='';
        while ($rst < $from && $n > 0) {
          $rst++; $n--; $qst++;
          die unless $rst <= $red && $qst <= $qed;
        }
        while ($rst >= $from && $rst <= $to && $n > 0) {
          $ss .= substr($qseq, $qst - 1, 1);
          $rst++; $n--; $qst++;
          die unless $rst <= $red && $qst <= $qed;
        }
        return undef if $rst > $from && $rst <= $to;
        return $ss if $rst > $to;
      } elsif ($op eq 'S') {	# ignore soft clipped sequences
      } elsif ($op eq 'I') {	# skip inserted sequences on the read
        $qst += $n;
        die unless $qst <= $qed;
      } elsif ($op eq 'D') {	# skip deleted sequences on the reference
        $rst += $n;
        die unless $rst <= $red;
      } elsif ($op eq 'H') {	# ignore hard clipped sequences
      } else {
        die "unknown operator $op \n";
      }
      return undef if $rst > $from;
    }
    return undef;
}

sub readBED {
  my $fBED = shift;
  my $sam = shift;

  my $fhBED = FileHandle->new($fBED) or die "$fBED:$!";
  my @region = ();

  while (my $line=<$fhBED>) {
    chomp $line;
    $line=~s/$//;

    my @val=split("\t", $line, -1);
    my $chr=$val[0];
    my $stt=$val[1];
    my $end=$val[2];
    my $name=$val[3];
    my @seq=split(",", $name);

    die "invalid line: $line\n" unless $stt < $end;

    my $len = $sam->length($chr) or die "$chr not found in the reference\n";
    die "invalid line: $line\n" unless 0 <= $stt && $stt <= $len;
    die "invalid line: $line\n" unless 0 <= $end && $end <= $len;

    foreach my $seq (@seq) {
      die "$line\nThe length of the candidate allele ($seq) is different from the specified length of ", $end - $stt, " bp.\n" unless length($seq) == $end - $stt;
    }
    push @region, { chr => $chr, stt => $stt, end => $end, ar_seq => \@seq };
  }
  $fhBED->close;
  return @region;
}

sub translate {
  my $seq = shift;
  my $len = length($seq);
  my $aa = '';

  for (my $i=0; $i < $len; $i+=3) {
    if ($len - $i >= 3) {
      $aa .= $aacode{substr($seq, $i, 3)};
    } else {
      $aa .= '-';
    }
  }
  return $aa;
}

sub usage {
  warn "Usage: $0 fasta bed bam\n";
  exit -1;
}
