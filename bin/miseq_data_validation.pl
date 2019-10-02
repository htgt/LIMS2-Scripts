#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw ($Bin);
use lib "$Bin/../lib";

use MiseqDataValidation qw(miseq_data_validator);

miseq_data_validator();
