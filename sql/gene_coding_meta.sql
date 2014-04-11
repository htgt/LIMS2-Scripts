CREATE TABLE gene_coding_meta
(
  id serial NOT NULL,
  marker_symbol text,
  ensembl_gene_id text NOT NULL,
  ensembl_exon_id text NOT NULL,
  cd_start integer NOT NULL,
  cd_end integer NOT NULL,
  chr_start integer NOT NULL,
  chr_end integer NOT NULL,
  chr_strand integer NOT NULL,
  CONSTRAINT gene_coding_meta_pkey PRIMARY KEY (id)
);
ALTER TABLE gene_coding_meta
  OWNER TO lims2;
