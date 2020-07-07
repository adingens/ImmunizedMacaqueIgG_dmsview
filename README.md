## Mutational antigenic profiling of purified IgG from immunized macaques

Adam S. Dingens, Jesse D. Bloom, in collaboration with Chris Cottrell and Andrew Ward

To view the data from XX using [`dms-view`](dms-view.github.io) you need three files: a [data file](process_ImmunoMacaque.csv), a [protein structure file](5fyl_trimer_renumber.pdb), and a [metadata file](ImmunizedMacaqueIgG.md).

### files for [`dms-view`](dms-view.github.io)
UPDATE ALL OF THIS
- data file: [`process_ImmunoMacaque.csv`](process_ImmunoMacaque.csv), made by [`process_ImmunoMacaque.py`](process_IDC561andmAbs.py)
- protein structure file: [`5fyl_trimer_renumber.pdb`](5fyl_trimer_renumber.pdb). This is the [5fyl pdb file](https://www.rcsb.org/structure/5FYL), with the antibody chains removed. The 321a insertion moved to 322a such that numbering matches the DMS data.  
- metadata file: [`ImmunizedMacaqueIgG.md`](ImmunizedMacaqueIgG.md)

### processing script
- [`process_ImmunoMacaque.py`](process_ImmunoMacaque.py)): Python script to create the [`process_ImmunoMacaque.csv`](process_ImmunoMacaque.csv) data file. This script retrieves the site and mutation data from the Hutch server.
