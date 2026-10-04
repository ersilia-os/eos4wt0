# Morgan fingerprints in binary form (radius 3, 2048 dimensions)

Represents a molecule as a 2,048-bit binary fingerprint recording which circular substructures are present within a radius of three bonds. Rogers and Hahn introduced extended-connectivity fingerprints specifically for structure-activity modelling rather than substructure search, iteratively updating atom identifiers with information from their neighbourhoods and hashing the result. Bit collisions mean a set bit can arise from more than one substructure, and the binary form records presence while discarding how many times a feature occurs.

This model was incorporated on 2023-12-01.Last packaged on 2026-08-31.

## Information
### Identifiers
- **Ersilia Identifier:** `eos4wt0`
- **Slug:** `morgan-binary-fps`

### Domain
- **Task:** `Representation`
- **Subtask:** `Featurization`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Descriptor`, `Fingerprint`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `2048`
- **Output Consistency:** `Fixed`
- **Interpretation:** 2048-bit binary fingerprint where each bit flags a circular substructure within radius three.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| feat_0000 | integer |  | Morgan fingerprint bit index 0 |
| feat_0001 | integer |  | Morgan fingerprint bit index 1 |
| feat_0002 | integer |  | Morgan fingerprint bit index 2 |
| feat_0003 | integer |  | Morgan fingerprint bit index 3 |
| feat_0004 | integer |  | Morgan fingerprint bit index 4 |
| feat_0005 | integer |  | Morgan fingerprint bit index 5 |
| feat_0006 | integer |  | Morgan fingerprint bit index 6 |
| feat_0007 | integer |  | Morgan fingerprint bit index 7 |
| feat_0008 | integer |  | Morgan fingerprint bit index 8 |
| feat_0009 | integer |  | Morgan fingerprint bit index 9 |

_10 of 2048 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos4wt0](https://hub.docker.com/r/ersiliaos/eos4wt0)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos4wt0.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos4wt0.zip)

### Resource Consumption
- **Model Size (Mb):** `1`
- **Environment Size (Mb):** `443`
- **Image Size (Mb):** `441.79`

**Computational Performance (seconds):**
- 10 inputs: `29.82`
- 100 inputs: `20.68`
- 10000 inputs: `46.67`

### References
- **Source Code**: [https://www.rdkit.org/docs](https://www.rdkit.org/docs)
- **Publication**: [https://doi.org/10.1021/ci100050t](https://doi.org/10.1021/ci100050t)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2010`
- **Ersilia Contributor:** [GemmaTuron](https://github.com/GemmaTuron)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [BSD-3-Clause](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos4wt0
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos4wt0
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
