#!/usr/bin/perl
package parse_restart_log;
use Moose;
use namespace::autoclean;

# Flags
has 'fq_splitting_started' => ( is => 'rw', isa => 'Bool', default => 0);
has 'alignment_started' => ( is => 'rw', isa => 'Bool', default => 0);
has 'mrsfast_calling_started' => ( is => 'rw', isa => 'Bool', default => 0);
has 'bwa_calling_started' => ( is => 'rw', isa => 'Bool', default => 0);

# Structures
has 'mrsfast_bam_list' => ( is => 'rw', isa => 'HashRef[ArrayRef[Str]', predicate => 'has_mrsfast_bams',);
has 'bwa_bam_list' => ( is => 'rw', isa => 'HashRef[ArrayRef[Str]', predicate => 'has_bwa_bams',);


# TODO: create a parsing program that takes the log file and fills appropriate flags and structures