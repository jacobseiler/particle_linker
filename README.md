Given a list of particle IDs, this script searches through a particle snapshot to find matches. The properties of the matched particles are stored and written out to create a matched *pseudo snapshot*.

## Notes 

The input match particle IDs are split up over multiple files. This program is parallelised and hence each rank will open up one snapshot subfile, then search all particle ID files for matches.  
