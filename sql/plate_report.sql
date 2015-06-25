WITH RECURSIVE well_hierarchy(process_id, input_well_id, output_well_id, path) AS (
    SELECT pr.id, pr_in.well_id, pr_out.well_id, ARRAY[pr_out.well_id]
    FROM processes pr
    LEFT OUTER JOIN process_input_well pr_in ON pr_in.process_id = pr.id
    JOIN process_output_well pr_out ON pr_out.process_id = pr.id
    WHERE pr_out.well_id IN (
        SELECT starting_well FROM well_list
    ) 
UNION ALL
    SELECT pr.id, pr_in.well_id, pr_out.well_id, path || pr_out.well_id
    FROM processes pr
    LEFT OUTER JOIN process_input_well pr_in ON pr_in.process_id = pr.id
    JOIN process_output_well pr_out ON pr_out.process_id = pr.id
    JOIN well_hierarchy ON well_hierarchy.input_well_id = pr_out.well_id
),
well_list(starting_well) AS (
    SELECT platewells.id FROM wells platewells, plates
    WHERE plates.name = 'HEPD0996_3'
    AND platewells.plate_id = plates.id
)
SELECT w.process_id, w.input_well_id, w.output_well_id, pd.design_id, w.path[1] "original_well", w.path
FROM well_hierarchy w, process_design pd
WHERE w.process_id = pd.process_id
GROUP BY w.process_id, w.input_well_id, w.output_well_id, pd.design_id,"original_well", w.path;
