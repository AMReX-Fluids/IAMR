import os
import numpy as np
import yt
import argparse, os
#from utils.etl_utils import etl
from utils.config import simulations_path, experiments_path
import matplotlib.pyplot as plt

class FlexibleData:
    def __init__(self, file_path):
      self.file_path = file_path
      self.member_variable = {}

    def read_file_inputc(self, count=0):
        with open(self.file_path, 'r') as file:
            for line in file:
                items = line.strip().split()
                if items:
                    if items[0] != '#':
                        self.member_variable[items[0]] = items[2]

def read_file_plt(simulation_path, variables, level_num, visualize, verbose):
    if not os.path.exists('output'):
        os.makedirs('output')
    bads = []
    count = 0
    level_num = int(level_num)
    for entry in os.scandir(simulation_path):
        if entry.is_dir() and 'plt' in entry.name and entry.name != 'plt.':
            plt_name = entry.name
            ds = yt.load(simulation_path + plt_name)
            ds.force_periodicity()
            for i in range(level_num):
                if verbose:
                    print("current level is")
                    print(i)
                dim = ds.domain_dimensions * 2**i
                dim[2] = 1
                ad = ds.covering_grid(level=i, left_edge=ds.domain_left_edge, dims=dim)
                level_path = os.path.join('output', f'level{i}')
                if not os.path.exists(level_path):
                    os.makedirs(level_path)
                for variable in variables:
                    full_path = os.path.join(level_path, variable)
                    if not os.path.exists(full_path):
                        os.makedirs(full_path)
                    try:
                        P_array = ad[variable].to_ndarray()
                        if verbose: 
                            print("here is the shape of numpy array")
                            print(P_array.shape)
                        
                        fig = plt.figure()
                        ax = fig.add_subplot(111, projection='3d')

                        # Create a 3D scatter plot
                        x, y, z = np.meshgrid(range(P_array.shape[0]), range(P_array.shape[1]), range(P_array.shape[2]))
                        ax.scatter(x, y, z, c=P_array, cmap='viridis')
                        ax.set_xlabel('X')
                        ax.set_ylabel('Y')
                        ax.set_zlabel('Z')
                        plt.title(f'3D Scatter Plot_{count} level{i} {variable} {plt_name}')
                        figure_name = f"fig{count}__level{i}__{variable}__{plt_name}.png"
                        figure_directory = full_path
                        os.makedirs(figure_directory, exist_ok=True)  # Create the "figure" directory if it doesn't exist
                        figure_path = os.path.join(figure_directory, figure_name)
                        if visualize == 'yes':
                            plt.savefig(figure_path)
                        full_path = os.path.join(level_path, variable)
                        np.save(full_path, P_array)
                    except ValueError:
                        print('ValueError in ', simulation_path, plt_name,f'{variable}. Skipping.')
                        bads.append(plt_name)
                        continue
                    count += 1


def main():
    parser = argparse.ArgumentParser(
                    prog='etl_experiment',
                    description='Extract-transform-load simulations in an experiment: provide experiment file with list of simulations in experiment ',
                    epilog='')
    
    parser.add_argument('experiment_id', help='Experiment ID')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('variables', nargs='*', help='Physical variables to process and visualize')
    parser.add_argument('level_num', type=int, help='Number of levels in the simulation grid to process')
    parser.add_argument('visualize', choices=['yes', 'no'], help='Enable or disable visualization')


    args = parser.parse_args()
    experiment_id = args.experiment_id
    verbose = args.verbose
    variables = args.variables
    level_num = args.level_num
    visualize = args.visualize
    experiment_path = experiments_path + experiment_id
    
    simulation_paths = []
    with open(experiment_path) as experiment:
        for i, line in enumerate(experiment.readlines()):
            if i == 0:
                print(line)
            else:
                simulation_paths.append(simulations_path + line.split()[0] + '/')

    ###############Open through plt files###################
    for simulation_path in simulation_paths:
        read_file_plt(simulation_path, variables, level_num, visualize, verbose)

if __name__ == "__main__":
    main()
