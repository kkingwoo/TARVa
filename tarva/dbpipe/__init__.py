"""
dbpipe: build the TARVa database from the data processing pipeline's output.

Four loaders, each the single writer of its own tables: annotate (gene,
transcript, exon), load-samples (individual), load-calls (site, call) and
load-expression (expression), with `all` to run them in the required order.
schema.py declares every table and holds the helpers they share.
"""
