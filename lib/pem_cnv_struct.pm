#!/usr/bin/perl
package pem_cnv_struct;
use Moose;
use namespace::autoclean;

extends 'base_cnv_struct';

has 'outStart' => (
	is => 'rw',
	isa => 'Int'
);

has 'outEnd' => (
	is => 'rw',
	isa => 'Int'
);

has 'sumProb' => (
	is => 'rw',
	isa => 'Num'
);

has 'support' => (
	is => 'rw',
	isa => 'Int'
);

has 'used' => (
	is => 'rw',
	isa => 'Bool',
	default => 0
);

__PACKAGE__->meta->make_immutable;

1;