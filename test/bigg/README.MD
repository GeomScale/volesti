# BiGG Models
This directory contains a single instance `e_coli_core.json`, a small *E. coli* core metabolic model used by the test and the benchmark suite. It is included directly in the repository for convenience.
## Additional Models
Additional models can be downloaded from the [BiGG database](http://bigg.ucsd.edu/). Download the JSON files of interest and place them in this directory. The benchmark will automatically pick up all `.json` files present in the directory.

You can also download a single model with id <model_id>, using the following command:
```bash
curl -O http://bigg.ucsd.edu/static/models/<model_id>.json

```