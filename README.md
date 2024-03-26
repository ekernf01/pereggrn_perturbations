This is a companion package to our [collection of perturbation data](https://github.com/ekernf01/perturbation_data).

```python
import pereggrn_perturbations
# Set this to point to the "perturbations" folder in the perturbation data collection. 
pereggrn_perturbations.set_data_path("path/to/perturbation_data/perturbations")
# What datasets are available?
pereggrn_perturbations.load_perturbation_metadata()
# Grab one
nakatake_et_al = pereggrn_perturbations.load_perturbation("nakatake") 
```

