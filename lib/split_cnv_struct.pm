#!/usr/bin/perl
package split_cnv_struct;
use Moose;
use namespace::autoclean;

extends 'base_cnv_struct';

has 'balSupport' => (
	is => 'rw',
	isa => 'Int'
);

has 'unbalSupport' => (
	is => 'rw',
	isa => 'Int'
);

has 'anchorEdit' => (
	is => 'rw',
	isa => 'Num'
);

has 'splitEdit' => (
	is => 'rw',
	isa => 'Num'
);

has 'used' => (
	is => 'rw',
	isa => 'Bool',
	default => 0
);

__PACKAGE__->meta->make_immutable;

1;