nucmer  -c 1000  -l 50 refgenome qrygenome       # Whole genome alignment. Any other alignment can also be used.
delta-filter -m -i 90 -l 100 out.delta >out.filtered.delta     # Remove small and lower quality alignments
show-coords -THrd out.filtered.delta >out.filtered.coords      # Convert alignment information to a .TSV format as required by SyRI
python3 /data/lix/software/syri/syri/bin/syri -c out.filtered.coords -d out.filtered.delta -r refgenome -q qrygenome
python3 /data/lix/software/syri/syri/bin/plotsr syri.out refgenome qrygenome -H 8 -W 5
