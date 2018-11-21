# Statistics
Version: 0.1

## Short Description
A collection of scripts for performing downstream analysis

## Description

Moderated t-statistics and regression based batch/covarite effect adjustment

## Key features

- Downstream analysis

## Approaches

- Metabolomics / Untargeted
- Metabolomics / Targeted

## Instrument Data Types

- MS / LC-MS
- MS / GC-MS

## Tool Authors

Smyth GK

## Container Contributors

- [Payam Emami](https://github.com/PayamEmami) (Stockholm U. - NBIS)

## Website

- http://msbi.ipb-halle.de/msbi/CAMERA/

## Git Repository

- https://github.com/sneumann/CAMERA

## Installation

No installation is needed on PhenoMeNal Cloud Research Environments.

For advanced Docker usage:

- Go to the directory where the dockerfile is.
- Create container from dockerfile:

```
docker build -t container-statistics .
```

Alternatively, pull from repo:

```
docker pull container-registry.phenomenal-h2020.eu/phnmnl/container-statistics
```

## Usage Instructions

On a PhenoMeNal Cloud Research Environment, go to "Statistical Analysis" category. 

Alternatively, to use locally through the docker image:

```
docker run --entrypoint=<script-name> container-registry.phenomenal-h2020.eu/phnmnl/container-statistics <arguments for script>
```

Where script can be any of: 

```
correctBatchEffect.r
univeriateLimma.r
```

## Publications

- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), e47.



