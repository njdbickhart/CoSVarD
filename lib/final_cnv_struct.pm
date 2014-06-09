#!/usr/bin/perl
package final_cnv_struct;
use Moose;
use namespace::autoclean;

extends 'base_cnv_struct';


has 'splitSupport' => (
	is => 'rw',
	isa => 'Num',
	default => '-1'
);

has 'vhSupport' => (
	is => 'rw',
	isa => 'Num',
	default => '-1'
);

has 'docConsist' => (
	is => 'rw',
	isa => 'Int',
	default => '0'
);

# default value is a negative number, meaning no copynumber found
has 'cn' => (
	is => 'rw',
	isa => 'Num',
	default => '-1'
);

has 'outStart' => (
	is => 'rw',
	isa => 'Int'
);

has 'outEnd' => (
	is => 'rw',
	isa => 'Int'
);

has 'eventNum' => (
	is => 'rw',
	isa => 'Int'
);

has 'splitCall' => (
	is => 'rw',
	isa => 'Str',
	lazy => 1,
	default => 'null'
);

has 'vhCall' => (
	is => 'rw',
	isa => 'Str',
	lazy => 1,
	default => 'null'
);

has 'docCall' => (
	is => 'rw',
	isa => 'Str',
	lazy => 1,
	default => 'null'
);

__PACKAGE__->meta->make_immutable;

1;