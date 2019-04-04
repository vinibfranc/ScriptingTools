package Animal::Lion;

use Animal::Cat;
use strict;

our @ISA = qw(Animal::Cat);

sub getSound {
    my ($self) = @_;
    return $self->{name}, " says rowwr"
}