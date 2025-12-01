## SimSpace conda env
python bench.py \
  --methods simspace \
  --grid-sizes 10 20 30 50 75 100 150 200 \
  --repeats 5 \
  --threads 1 \
  --outfile results_simspace.csv

## sccube conda env
python bench.py \
  --methods sccube \
  --grid-sizes 10 20 30 50 75 100 150 200 \
  --repeats 5 \
  --threads 1 \
  --outfile results_sccube.csv

python bench.py \
  --methods scmultisim \
  --grid-sizes 10 15 20 25 \
  --repeats 5 \
  --threads 1 \
  --outfile results_scmultisim.csv

## Make plots
python plot_results.py --csv results.csv --outdir plots