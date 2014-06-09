#!/usr/bin/perl
package base_cnv_struct;
use Moose;
use Moose::Util::TypeConstraints;
use namespace::autoclean;

has 'chr' => (
	is => 'rw',
	isa => 'Str',
	required => 1
);

has 'insStart' =>(
	is => 'rw',
	isa => 'Int',
	required => 1
);

has 'insEnd' => (
	is => 'rw',
	isa => 'Int',
	required => 1
);

enum 'My:Enum:cnvType', [qw( d v s sv ds dsv )];
# Values: d, v, s, sv or dsv
has 'type' => (
	is => 'rw',
	isa => 'My:Enum:cnvType',
	required => 1
);

enum 'cnvCall' => qw( ins inv del gain mult );
# Values: ins, inv, del, gain
has 'call' => (
	is => 'rw',
	isa => 'cnvCall'
);

__PACKAGE__->meta->make_immutable;

1;