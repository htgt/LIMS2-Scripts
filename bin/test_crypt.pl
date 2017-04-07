#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use Data::Serializer;

my $ser = Data::Serializer->new(
    serializer => 'Data::Dumper',
    digester => 'SHA-256',
    cipher => 'Blowfish',
    secret => 'V7}~q+ztoAk|{KMF9hN?5ZxlsrUNFVI5wJM{d%izT8_~b704zd&%dB5CSwP5(Hh',
    compress => 0,
);

my $test = $ser->freeze({ test => 'pineapple' });

my $one = $ser->thaw('^Data::Dumper|Blowfish|SHA-256|hex|^53616c7465645f5f41c97a20953a07b2cad63b75c8e5b8c5161ac6dde424bd50b6d5cced1022889f57001cb562506bf1852f0349913b43a6b0184c7d5d1dba5a0325372b36a70615874f7a8a447c5238d1aa587d10c4761cf2a7603f2714301bca391df1f3f2dda4e74a5ec6094644f9b9ea352df4dea94d4e705d6b55ab81e095827e31cf2a1cf0357574b448e396f8d7878406ba71c6a56a6145e01be6592588349ac8b781945b197752d503e26f34571800977d739249c515387869e89987');

my $two = $ser->thaw('^Data::Dumper|Blowfish|SHA-256|hex|^53616c7465645f5f691ff6e0530195563e73ab89a6c6e07e11433349fa7c5e055d63b3851767a3dca48c8ba169f586ba1665bb1e4ff4df2c3c2923bbfdc251fb7e62a3d71dae855fc8b69e6d9a5657367c3809ce22049cedd939b5f5079ce568a256cbefbc16ba27816dc3eb7be8311fad9035d237f1abdba1de603284a2e6231acd6990150b0bdfaa2fa9107d00eeea6cdd931d7798ac988aca53376ef44dec9eecbf59a087a0f6709cf3e2d82db33068b71de00ce77933def4c4cb662551a2');
use Data::Dumper;
print Dumper ($one) . "\n";
print "New one - " . Dumper ($two) . "\n";

1;

__END__
