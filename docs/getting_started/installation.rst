Installation
============


Install from github::

    pip install -e http://github.com/akhanf/snakedwi


Running the app
===============

Do a dry-run first (``-n``) and simply print (``-p``) what would be run::

    snakedwi /path/to/bids/dir /path/to/output/dir participant -np

Run the app, using all cores::

    snakedwi /path/to/bids/dir /path/to/output/dir participant --cores all

If any workflow rules require containers, then run with the ``--use-singularity`` option.


Generating a report
===================

After your processing is complete, you can use snakemake's ``--report`` feature to generate
an HTML report. This report will include a graph of all the jobs run, with clickable nodes
to inspect the shell command or python code used in each job, along with the config files and
run times for each job. Workflows may also contain append images for quality assurance or to
summarize outputs, by using the ``report(...)`` function on any snakemake output.

To generate a report, run::

    snakedwi /path/to/bids/dir /path/to/output/dir participant --report



