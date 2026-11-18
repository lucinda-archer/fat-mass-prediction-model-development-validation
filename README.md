# Development and validation of a prediction model for fat mass in adults

## Overview

This repository contains the Stata code used for the development and validation of a prediction model for fat mass in adults.

The code covers model development, internal validation, internal-external cross-validation and external validation in five countries.

## Repository structure

The `code/` folder contains the following Stata scripts:

| File | Analysis |
|---|---|
| `01_model_development.do` | Model development |
| `02_bootstrap_internal_validation.do` | Bootstrap internal validation |
| `03_internal_external_cross_validation.do` | Internal-external cross-validation |
| `04_external_validation_namibia.do` | External validation in Namibia |
| `05_external_validation_kuwait.do` | External validation in Kuwait |
| `06_external_validation_usa.do` | External validation in the USA |
| `07_external_validation_singapore.do` | External validation in Singapore |
| `08_external_validation_south_korea.do` | External validation in South Korea |

Each script relates to a separate analysis and can be run independently.

## Data

The datasets used in these analyses are not included in this repository. Access to the data is subject to the relevant data access, governance and ethical requirements.

## Software

All analyses were conducted using Stata 19.5. The scripts may require access to the relevant datasets and any user-written Stata commands specified within them.
