# Analysis Containers

This directory holds analysis outputs grouped by analysis type.

## Primary result (current)

`triplet_contact_story_100ns/`

- `plots/`: final triplet contact-story figures (mean trace across replicates + pooled occupancy heatmap)
- `tables/`: summary, occupancy, timing, and mean-trace CSVs used by the plots
- `config/`: run configuration (triplet definitions)

Start here:

- `triplet_contact_story_100ns/plots/`
- `triplet_contact_story_100ns/tables/selection_summary.csv`

## Convention

For future analyses, create one sibling folder per analysis type:

`results/analysis/<analysis_name>/{plots,tables,config}/`

This keeps exploratory outputs separate from primary deliverables.
