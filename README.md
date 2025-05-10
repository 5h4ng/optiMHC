# README

## Usage

### Command-line Execution (Modern CLI)

You can run optiMHC using either a YAML configuration file or by specifying parameters directly via the command line.

#### Using a YAML Configuration File

```bash
optimhc pipeline --config /path/to/config.yaml
```

Run in experiment mode:

```bash
optimhc experiment --config /path/to/config.yaml
```

#### Using Direct Command-Line Parameters

```bash
optimhc pipeline \
  --input-type pepxml \
  --input-file ./data/YE_20180428_SK_HLA_A0202_3Ips_a50mio_R1_01.pep.xml \
  --decoy-prefix DECOY_ \
  --output-dir ./results \
  --visualization \
  --num-processes 32 \
  --allele HLA-A*02:02 \
  --feature-generator 'name: Basic' \
  --feature-generator 'name: DeepLC' \
  --test-fdr 0.01 \
  --model Percolator
```

You can mix and override YAML config values with CLI parameters as needed.

#### Full CLI Help

```bash
optimhc --help
optimhc pipeline --help
optimhc experiment --help
```

---

## YAML Configuration File

The pipeline is configured using a YAML file. This file defines the input settings, the list of feature generators, rescore parameters, and (optionally) experiment configurations. Below you will find a table summarizing the main configuration parameters along with examples and descriptions.

### Configuration Parameters

| Parameter           | Type               | Example                                                          | Description                                                             |
|---------------------|--------------------|------------------------------------------------------------------|-------------------------------------------------------------------------|
| `experimentName`    | String             | `classI_example`                                                 | Name of the experiment and output subdirectory name.                  |
| `inputType`         | String             | `pepxml`                                                         | Type of input file. Supported values: `pepxml`, `pin`.                  |
| `inputFile`         | String or List     | `./data/YE_20180428_SK_HLA_A0202_3Ips_a50mio_R1_01.pep.xml`        | Path(s) to the input PSM file(s).                                       |
| `decoyPrefix`       | String             | `DECOY_`                                                         | Prefix used to identify decoy sequences.                              |
| `outputDir`         | String             | `./results`                                                      | Base directory where output files, logs and figures are stored.         |
| `visualization`     | Boolean            | `True`                                                           | Enable or disable generation of visualization plots.                  |
| `removePreNxtAA`    | Boolean            | `False`                                                          | Remove pre/post neighboring amino acids in sequence processing.         |
| `numProcesses`      | Integer            | `32`                                                             | Number of parallel processes to use.                               |
| `showProgress`      | Boolean            | `True`                                                           | Show progress information during execution.                           |
| `modificationMap`   | Dictionary         | `{ '147.035385': 'UNIMOD:35' }`                                    | Maps modification masses to their 'UNIMOD' identifiers. See https://www.unimod.org/ for details           |
| `allele`            | List               | `[HLA-A*02:02]`                                                  | List of alleles for which predictions will be computed.                 |
| `featureGenerator`  | List of Dictionaries | See table below                                                  | List of feature generator configurations (each with a `name` and optional `params`). |
| `rescore`           | Dictionary         | See table below                                                 | Rescore settings including FDR threshold, model and number of jobs.     |

---

### Feature Generator Configurations

Each feature generator is specified with its `name` and an optional `params` subsection. Some common generators include:

| Generator Name      | Example Parameters                                                                                                           | Description                                                                              |
|---------------------|------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------|
| `Basic`             | N/A                                                                                                                          | Generates basic sequence features.                                                       |
| `SpectraSimilarity` | `mzmlDir: ./data`<br>`spectrumIdPattern: (.+?)\.\d+\.\d+\.\d+`<br>`model: AlphaPeptDeep_ms2_generic`<br>`collisionEnergy: 28`<br>`instrument: LUMOS`<br>`tolerance: 20`<br>`numTopPeaks: 36`<br>`url: koina.wilhelmlab.org:443` | Computes features based on the similarity between experimental spectra and predicted spectra. See more options on https://koina.proteomicsdb.org/ |
| `DeepLC`            | `calibrationCriteria: expect`<br>`lowerIsBetter: True`<br>`calibrationSize: 0.1`                                             | Creates retention time predictions by calibrating using DeepLC.                           |
| `OverlappingPeptide`     | `minOverlapLength: 7`<br>`minLength: 7`<br>`maxLength: 20`<br>`overlappingScore: expect`                                           | Generates overlapping peptide features for grouping similar peptides.                         |
| `PWM`               | `class: I`                                                                                                                   | Generates position weight matrix features for MHC class I and class II peptides.                        |
| `MHCflurry`         | N/A                                                                                                                          | Predicts class I binding affinities using the MHCflurry model.                                     |
| `NetMHCpan`         | N/A                                                                                                                          | Predicts class I peptide-MHC binding affinity using NetMHCpan.                                     |
| `NetMHCIIpnan` | N/A                                                                                                                          | Predicts class II peptide-MHC binding affinity using NetMHCIIpan.                                   |
---

### Rescore Settings

Rescore parameters control how the rescoring step is executed and include:

| Parameter  | Type    | Example     | Description                                                                |
|------------|---------|-------------|----------------------------------------------------------------------------|
| `testFDR`  | Float   | `0.01`      | The false-discovery rate threshold at which to evaluate the learned models.       |
| `model`    | String  | `Percolator`| Model to use for rescoring (valid options include `Percolator`, `XGBoost`, or `RandomForest`). |
| `numJobs`  | Integer | `4`         |The number of parallel jobs to run. This value is passed to Scikit-learn's n_jobs parameter to control parallelism for model training or scoring. Set to -1 to use all available CPU cores.                   |

---

### Experiment Mode (Optional)

When running in experiment mode, the YAML file can include an `experiments` section that defines multiple experiments with different feature combinations and models. Each experiment entry uses the following parameters:

| Parameter | Type   | Example                                                                                       | Description                                                              |
|-----------|--------|-----------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|
| `name`    | String | `"Baseline"`                                                                                  | Name of the experiment; also used to create a dedicated output directory. |
| `source`  | List   | `["Original", "Basic", "DeepLC"]`                                                             | List of feature source names to include in the experiment.               |
| `model`   | String | `"Percolator"` or `"XGBoost"`                                                                 | Rescoring model to use for the experiment.                                |

---

### Example YAML Configuration

Below is an example YAML configuration for class I based on the latest pipeline version:

```yaml
experimentName: classI_example
inputType: pepxml
inputFile:
  - ./data/YE_20180428_SK_HLA_A0202_3Ips_a50mio_R1_01.pep.xml
decoyPrefix: DECOY_
outputDir: ./results
visualization: True
removePreNxtAA: False
numProcesses: 32
showProgress: True
modificationMap:
  '147.035385': 'UNIMOD:35'

# Allele settings
allele:
  - HLA-A*02:02

# Feature generator configurations
featureGenerator:
  - name: Basic
  - name: SpectraSimilarity
    params:
      mzmlDir: ./data
      spectrumIdPattern: (.+?)\.\d+\.\d+\.\d+
      model: AlphaPeptDeep_ms2_generic
      collisionEnergy: 28
      instrument: LUMOS
      tolerance: 20  
      numTopPeaks: 36
      url: koina.wilhelmlab.org:443
  - name: DeepLC
    params:
      calibrationCriteria: expect
      lowerIsBetter: True
      calibrationSize: 0.1
  - name: OverlappingPeptide
    params:
      minOverlapLength: 7
      minLength: 7
      maxLength: 20
      overlappingScore: expect 
  - name: PWM
    params:
      class: I
  - name: MHCflurry
  - name: NetMHCpan

# Rescore settings
rescore:
  testFDR: 0.01
  model: Percolator
  numJobs: 4
```

---

## Try It Yourself

You can also run the test suite as a full example:

```bash
pytest tests/
```

