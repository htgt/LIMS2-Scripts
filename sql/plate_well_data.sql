WITH RECURSIVE well_hierarchy(root_well_id, process_id, process_type, input_well_id, output_well_id, output_well_name, output_plate_name, output_plate_type, depth, design_id, short_arm_design_id, cassette, cassette_resistance, cassette_promoter, backbone, recombinase, cell_line, gene_id, gene_symbol, crispr_id, nuclease ) AS (
    SELECT pr_out.well_id, pr.id, pr.type, pr_in.well_id, pr_out.well_id, output_well.name, output_plate.name, output_plate.type_id, 0, pr.design_id, pr.short_arm_design_id, pr.cassette, pr.cassette_resistance, pr.cassette_promoter, pr.backbone, pr.recombinase, pr.cell_line, pr.gene_id, pr.gene_symbol, pr.crispr_id, pr.nuclease
    FROM process_data pr
    LEFT OUTER JOIN process_input_well pr_in ON pr_in.process_id = pr.id
    JOIN process_output_well pr_out ON pr_out.process_id = pr.id
    JOIN wells output_well ON output_well.id = pr_out.well_id
    JOIN plates output_plate ON output_plate.id = output_well.plate_id
    WHERE pr_out.well_id IN (
        SELECT starting_well FROM well_list
    ) 
UNION ALL
    SELECT wh.root_well_id, pr.id, pr.type, pr_in.well_id, pr_out.well_id, output_well.name, output_plate.name, output_plate.type_id, wh.depth + 1, pr.design_id, pr.short_arm_design_id, pr.cassette, pr.cassette_resistance, pr.cassette_promoter, pr.backbone, pr.recombinase, pr.cell_line, pr.gene_id, pr.gene_symbol, pr.crispr_id, pr.nuclease
    FROM process_data pr
    LEFT OUTER JOIN process_input_well pr_in ON pr_in.process_id = pr.id
    JOIN process_output_well pr_out ON pr_out.process_id = pr.id
    JOIN wells output_well ON output_well.id = pr_out.well_id
    JOIN plates output_plate ON output_plate.id = output_well.plate_id
    JOIN well_hierarchy wh ON wh.input_well_id = pr_out.well_id
),
process_data ( id, type, cassette, cassette_resistance, cassette_promoter, backbone, design_id, short_arm_design_id, cell_line, recombinase, gene_id, gene_symbol, crispr_id, nuclease ) as (
    SELECT p.id, p.type_id, c.name, c.resistance, c.promoter, b.name, pd.design_id, psd.design_id, cl.name, pr.recombinase_id, gene.design_gene_id, gene.design_gene_symbol, pcr.crispr_id, n.name
    FROM processes p
    LEFT OUTER JOIN process_cassette pc on pc.process_id = p.id
    LEFT OUTER JOIN cassettes c on c.id = pc.cassette_id

    LEFT OUTER JOIN process_backbone pb on pb.process_id = p.id
    LEFT OUTER JOIN backbones b on b.id = pb.backbone_id

    LEFT OUTER JOIN process_design pd on pd.process_id = p.id
    LEFT OUTER JOIN process_global_arm_shortening_design psd on psd.process_id = p.id

    LEFT OUTER JOIN process_crispr pcr on pcr.process_id = p.id

    LEFT OUTER JOIN process_recombinase pr on pr.process_id = p.id

    LEFT OUTER JOIN process_cell_line pcl on pcl.process_id = p.id
    LEFT OUTER JOIN cell_lines cl on cl.id = pcl.cell_line_id

    LEFT OUTER JOIN process_nuclease pn on pn.process_id = p.id
    LEFT OUTER JOIN nucleases n on n.id = pn.nuclease_id

    LEFT OUTER JOIN (
        select distinct on (design_id) design_id, design_gene_id, design_gene_symbol
        FROM summaries
    ) as gene 
    ON pd.design_id = gene.design_id
),
well_list(starting_well) AS (
    SELECT wells.id
    FROM wells inner join plates on plates.id = wells.plate_id
    WHERE plates.id = 6665 
)
SELECT w.*
FROM well_hierarchy w
ORDER BY root_well_id, depth
;
