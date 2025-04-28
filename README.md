# Morgan fingerprints in binary form (radius 3, 2048 dimensions)

The Morgan Fingerprints are one of the most widely used molecular representations. They are circular representations (from an atom,search the atoms around with a radius n) and can have thousands of features. This implementation uses the RDKit package and is done with radius 3 and 2048 dimensions, providing a binary vector as output. For Morgan counts, see model eos5axz.

This model was incorporated on 2023-12-01.

## Information
### Identifiers
- **Ersilia Identifier:** `eos4wt0`
- **Slug:** `morgan-binary-fps`

### Domain
- **Task:** `Representation`
- **Subtask:** `Featurization`
- **Biomedical Area:** `Any`
- **Target Organism:** `Not Applicable`
- **Tags:** `Descriptor`, `Fingerprint`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `2048`
- **Output Consistency:** `Fixed`
- **Interpretation:** Binary vector representing the SMILES

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| dim_0000 | integer |  | Morgan fingerprint bit index 0 |
| dim_0001 | integer |  | Morgan fingerprint bit index 1 |
| dim_0002 | integer |  | Morgan fingerprint bit index 2 |
| dim_0003 | integer |  | Morgan fingerprint bit index 3 |
| dim_0004 | integer |  | Morgan fingerprint bit index 4 |
| dim_0005 | integer |  | Morgan fingerprint bit index 5 |
| dim_0006 | integer |  | Morgan fingerprint bit index 6 |
| dim_0007 | integer |  | Morgan fingerprint bit index 7 |
| dim_0008 | integer |  | Morgan fingerprint bit index 8 |
| dim_0009 | integer |  | Morgan fingerprint bit index 9 |

_10 of 2048 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos4wt0](https://hub.docker.com/r/ersiliaos/eos4wt0)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos4wt0.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos4wt0.zip)

### Resource Consumption
- **Model Size (Mb):** `1`
- **Environment Size (Mb):** `435`


### References
- **Source Code**: [https://www.rdkit.org/docs](https://www.rdkit.org/docs)
- **Publication**: [https://pubmed.ncbi.nlm.nih.gov/20426451/](https://pubmed.ncbi.nlm.nih.gov/20426451/)
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
