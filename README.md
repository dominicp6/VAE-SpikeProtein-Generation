# VAE-SpikeProtein-Generation
Implementation of "Generating Sentences from a Continuous Space" paper using pytorch and torchtext. Generates SARS-COV2 Spike proteins based on the GISAID dataset.

## Prerequisites
* Python 2.7
* Pytorch 0.4.0

## Data
GISAID data not included; must be stored locally.

## Training
To train the model from the beginning, run:
```
python main.py
```
To resume training from a saved model, run:
```
python main.py --resume_training=True
```
## Testing
To generate samples from a saved model, run:
```
python main.py --to_train=False
```
