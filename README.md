## **Usage**

### **Command-line Execution**

Run the pipeline with a YAML configuration file:

```bash
python run_pipeline.py /path/to/config.yaml
```

Run in experiment mode to test multiple feature combinations:

```bash
python run_pipeline.py /path/to/config.yaml --mode experiment
```

### **Command-line Help**

To see available options:

```bash
python run.py -h
```

Output:

```
usage: run_pipeline.py [-h] config

Run the analysis pipeline with the specified YAML configuration.

positional arguments:
  config      Path to the YAML configuration file.

optional arguments:
  -h, --help  show this help message and exit
```

---

## **Configuration**

The pipeline is configured via a YAML file. Below is an example configuration:

```yaml
# Input settings
input_type: pepxml
input_files:
  - /path/to/file1.pep.xml
  - /path/to/file2.pep.xml

# Output settings
output_dir: /path/to/output
allele: [HLA-A*02:06]
test_fdr: 0.01
visualization: true

# Global parameters
global_parameters:
  remove_modification: true
  remove_pre_nxt_aa: false
  n_processes: 32
  show_progress: true

# Feature generator configurations
feature_generators:
  - type: LadderPeptide
    min_overlap_length: 8
    min_entropy: 0
  - type: Basic
  - type: PWM
    mhc_class: I
  - type: MHCflurry
  - type: NetMHCpan
    mode: all
  - type: NetMHCIIpan
    mode: best
  - type: DeepLC
    calibration_criteria_column: spscore
    lower_score_is_better: false
    calibration_set_size: 0.1
  - type: SpectraSimilarity
    mzML_dir: /path/to/mzML/files
    spectrum_id_pattern: (.+?)\.\d+\.\d+\.\d+
    fragmentation_method: CID
    tolerance_ppm: 20

# Rescore settings
rescore:
  test_fdr: 0.01
  model: percolator
  n_jobs: 4
```

## **Experiment Mode**

Experiment mode allows you to run multiple experiments with different feature combinations in a single execution. Each experiment uses a different set of features from the generated feature set.

### **Experiment Configuration**

To configure experiment mode, add an `experiments` section to your YAML file:

```yaml
# Standard configuration as above
# ...

# Experiments configuration
experiments:
  - name: "Baseline"
    source: ["Original"]
    model: "percolator"
  
  - name: "Basic+RT"
    source: ["Original", "Basic", "DeepLC"]
    model: "percolator"
  
  - name: "Complete"
    source: ["Original", "Basic", "DeepLC", "MHCflurry", "NetMHCpan", "PWM", "LadderPeptide"]
    model: "xgboost"
```

### **Experiment Parameters**

Each experiment in the `experiments` array can have the following parameters:

- `name`: Name of the experiment (used for output directory)
- `source`: List of feature sources to use in this experiment
- `model`: ML model to use for rescoring ("percolator", "xgboost", or "random_forest")

### **Example Full Experiment Configuration**

```yaml
# Input settings
input_type: pepxml
input_files:
  - /path/to/file1.pep.xml
  - /path/to/file2.pep.xml

# Output settings
output_dir: /path/to/experiment_results
allele: [HLA-A*02:06]
visualization: true
score: spscore

# Global parameters
global_parameters:
  remove_modification: true
  remove_pre_nxt_aa: false
  n_processes: 32
  show_progress: true

# Feature generator configurations
feature_generators:
  - type: Basic
  - type: DeepLC
    calibration_criteria_column: spscore
    lower_score_is_better: false
    calibration_set_size: 0.1
  - type: MHCflurry
  - type: NetMHCpan
    mode: all
  - type: LadderPeptide
    min_overlap_length: 7
    min_entropy: 0
  - type: PWM
    mhc_class: I

# Experiment definitions
experiments:
  - name: Baseline
    source: ["Original"]
    model: percolator
  
  - name: Basic+RT
    source: ["Original", "Basic", "DeepLC"]
    model: percolator
  
  - name: Basic+RT+Affinity
    source: ["Original", "Basic", "DeepLC", "MHCflurry", "NetMHCpan"]
    model: percolator
  
  - name: Basic+RT+Affinity+PWM
    source: ["Original", "Basic", "DeepLC", "MHCflurry", "NetMHCpan", "PWM"]
    model: percolator
  
  - name: Basic+RT+Affinity+Ladder
    source: ["Original", "Basic", "DeepLC", "MHCflurry", "NetMHCpan", "LadderPeptide", "LadderGroupFeatures"]
    model: percolator
  
  - name: Full_Model
    source: ["Original", "Basic", "DeepLC", "MHCflurry", "NetMHCpan", "PWM", "LadderPeptide", "LadderGroupFeatures"]
    model: xgboost
```

---
