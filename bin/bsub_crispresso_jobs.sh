#!/usr/local/bin/bash
PERL5LIB=/software/perl-5.16.2
BIN="$(readlink -f "$(dirname "$0")")"
echo "$@"
$BIN/bsub_crispresso_jobs.pl "$@"
