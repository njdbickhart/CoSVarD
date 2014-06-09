#!/usr/bin/perl
package doc_cnv_struct;
use Moose;
use namespace::autoclean;


extends 'base_cnv_struct';

# default value is a negative number, meaning no copynumber found
has 'cn' => (
	is => 'rw',
	isa => 'Num',
	default => '-1'
);

# I am going to construct outer coordinates by default based on the limitation of the 1kb window crunch in the 
# DOC algorithm
has 'outStart' => (
	is => 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_build_out_start'
);

has 'outEnd' => (
	is => 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_build_out_end'
);

has 'used' => (
	is => 'rw',
	isa => 'Bool',
	default => 0
);

sub _build_out_start {
	my $self = shift;
	return ($self->insStart() - 1000 < 0 ? 0 : $self->insStart() - 1000);
}

sub _build_out_end {
	my $self = shift;
	return ($self->insEnd() + 1000);
}

__PACKAGE__->meta->make_immutable;

1;