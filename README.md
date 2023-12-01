# Morgan Fingerprints

The Morgan Fingerprints are one of the most widely used molecular representations. They are circular representations (from an atom,search the atoms around with a radius n) and can have thousands of features. This implementation uses the RDKit package and is done with radius 3 and 2048 dimensions,providing a binary vector as output. For Morgan counts, see eos5axz.

## Identifiers

* EOS model ID: `eos4wt0`
* Slug: `morgan-fps`

## Characteristics

* Input: `Compound`
* Input Shape: `Single`
* Task: `Representation`
* Output: `Descriptor`
* Output Type: `Integer`
* Output Shape: `List`
* Interpretation: Binary vector representing the SMILES

## References

* [Publication](https://pubmed.ncbi.nlm.nih.gov/20426451/)
* [Source Code](https://www.rdkit.org/docs)
* Ersilia contributor: [GemmaTuron](https://github.com/GemmaTuron)

## Ersilia model URLs
* [GitHub](https://github.com/ersilia-os/eos4wt0)
* [AWS S3](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos4wt0.zip)
* [DockerHub](https://hub.docker.com/r/ersiliaos/eos4wt0) (AMD64, ARM64)

## Citation

If you use this model, please cite the [original authors](https://pubmed.ncbi.nlm.nih.gov/20426451/) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).

## License

This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a BSD-3.0 license.

Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission!