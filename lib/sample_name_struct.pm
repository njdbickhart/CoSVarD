#!/usr/bin/perl
package sample_name_struct;
use Moose;
use namespace::autoclean;

has 'sample_name' => (
	is => 'rw',
	isa => 'Str',
);

has 'fastq_1' => (
	is => 'rw',
	isa => 'Str',
);

has 'fastq_2' => (
	is => 'rw',
	isa => 'Str',
);

has 'lib_type' => (
	is => 'rw',
	isa => 'Str',
);

has 'insert_size' => (
	is => 'rw',
	isa => 'Int',
);

has 'insert_stdev' => (
	is => 'rw',
	isa => 'Num',
);

has 'split_fq_list' => (
	is => 'rw',
	isa => 'ArrayRef[Str]',
	predicate => 'has_fq',
	clearer => 'clear_mrsfastq',
);

has 'mrsfast_bam_list' => (
	is => 'rw',
	isa => 'ArrayRef[Str]',
	predicate => 'has_bam',
);


__PACKAGE__->meta->make_immutable;
1;