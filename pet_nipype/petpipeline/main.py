#!/usr/bin/env python3

import sys
import argparse
import yaml
import os
from PETPipeline import PETPipeline
from config import _EnvConfig, \
                   _MotionCorrectionConfig, \
                   _PartialVolumeCorrectionConfig, \
                   _ReconAllConfig, \
                   _CoregistrationConfig \

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", default=os.path.realpath('config.yaml'))
    parser.add_argument("-e", "--experiment_dir", default=None,
                        help="The experiment directory")
    parser.add_argument("-o", "--output_dir", default=None,
                        help="The output directory (relative to the experiment directory)")
    parser.add_argument("-w", "--working_dir", default=None,
                        help="The working directory (relative to the experiment directory)")
    parser.add_argument("-d", "--data_dir", default=None,
                        help="The data directory containing the bids dataset")                        
    return parser.parse_args()

def parse_yaml(args):
    config = None
    with open(args.config,"r") as stream:
        try:
            config = yaml.load(stream, yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc)
    # override vars from config file with command line params 
    if args.experiment_dir is not None:
        config['environment']['experiment_dir'] = args.experiment_dir
    if args.output_dir is not None:
        config['environment']['output_dir'] = args.output_dir
    if args.working_dir is not None:
        config['environment']['working_dir'] = args.working_dir
    if args.data_dir is not None:
        config['environment']['data_dir'] = args.data_dir
    return config
   
def main(argv):
    args = parse_args(argv)

    if args.config:
        config = parse_yaml(args)
        env_config = _EnvConfig(**config['environment'])
        motion_correction_config = _MotionCorrectionConfig(**config['motion_correction']) 
        coregistration_config = _CoregistrationConfig(**config['coregistration'])
        reconall_config = _ReconAllConfig(**config['reconall'])
        pvc_config = _PartialVolumeCorrectionConfig(**config['partial_volume_correction'])

    else:
        env_config = _EnvConfig(experiment_dir=args.experiment_dir, \
                        output_dir=args.output_dir, \
                        working_dir=args.working_dir, \
                        data_dir=args.data_dir)
    

    pipeline = PETPipeline(env_config=env_config,
                           motion_correction_config = motion_correction_config, \
                           coregistration_config = coregistration_config, \
                           reconall_config = reconall_config, \
                           pvc_config = pvc_config)
    pipeline.PETWorkflow()
    pipeline.run()
    
if __name__ == "__main__":
    sys.exit(main(sys.argv))
    print(os.getcwd())
