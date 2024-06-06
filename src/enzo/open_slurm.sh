#!/bin/bash

srun --job-name=compile --pty --nodes=1  --cpus-per-task=16 -p gpu --gpus=1 bash
