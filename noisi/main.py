from __future__ import print_function
import sys
import os
import io
import click
import json
import time

from noisi.scripts.source_grid import setup_sourcegrid as setup_sgrid
from noisi.scripts.run_correlation import run_corr
from noisi.util.prepare_sem_input import prepare_specfem_input
from noisi.scripts.run_measurement import run_measurement
from noisi.scripts.run_adjointsrcs import run_adjointsrcs
from noisi.scripts.run_kernel import run_kern
from noisi.scripts.run_preprocessing import run_preprocessing
from noisi.scripts.run_preprocessing_data import run_preprocess_data
from noisi.scripts.assemble_gradient import assemble_ascent_dir
from noisi.util.setup_new import setup_proj
@click.group()
def run():
    """
    Main routine for noise correlation modeling and noise source inversion.
    """
    pass
    

###########################################################################
### Setting up a new project
###########################################################################
@run.command(help='Initialize a new project.')
@click.argument('project_name')
def setup_project(project_name):
    if os.path.exists(project_name):
        click.echo('Project exists already, must give it a new name.')
        exit()
    else:
        setup_proj(project_name)
    
    click.secho("Copied default config.json to project directory, please edit.")


###########################################################################
### Setting up the discretized source grid
###########################################################################  
@run.command(help='Determine the source grid and get specfem STATIONS file.')
@click.argument('project_path')
def setup_sourcegrid(project_path):
    setup_sgrid(os.path.join(project_path,'config.json'))


###########################################################################
### Initialize a source model
###########################################################################   
@run.command(help='Initialize a new source model.')
@click.argument('source_model')
def setup_source(source_model):

    if os.path.exists(source_model):
        click.echo('Source exists already, must give it a new name.')
        exit()

    if not os.path.exists('config.json'):
        click.echo('No config file for project found \
        (detailing e.g. source grid). Run setup_project first.')
        exit()

    os.makedirs(os.path.join(source_model,'step_0'))
    os.mkdir(os.path.join(source_model,'observed_correlations'))
  
    for d in ['adjt','grad','corr','kern']:
        os.mkdir(os.path.join(source_model,'step_0',d))

    from . import _ROOT
    
    with io.open(os.path.join(_ROOT,'config','source_config.json'),'r') as fh:
        conf = json.loads(fh.read())
        conf['date_created'] = str(time.strftime("%Y.%m.%d"))
        conf['project_name'] = os.path.basename(os.getcwd())
        conf['project_path'] = os.getcwd()
        conf['source_name'] = source_model
        conf['source_path'] = os.path.abspath(source_model)
        
    with io.open(os.path.join(source_model,'source_config.json'),'w') as fh:
        cf = json.dumps(conf,sort_keys=True, indent=4, separators=(",", ": "))
        fh.write(cf)
    
    with io.open(os.path.join(_ROOT,'config','measr_config.json'),'r') as fh:
        conf = json.loads(fh.read())
        conf['date_created'] = str(time.strftime("%Y.%m.%d"))
       
    with io.open(os.path.join(source_model,'measr_config.json'),'w') as fh:
        cf = json.dumps(conf,sort_keys=True, indent=4, separators=(",", ": "))
        fh.write(cf)
    
    from . import _ROOT
    os.system('cp {} {}'.format(os.path.join(_ROOT,'jnotebks/\
setup_noisesource.ipynb'),
    source_model))  
    os.system('cp {} {}'.format(os.path.join(_ROOT,'util/setup_noisesource.py'),
    source_model))
    os.system('cp {} {}'.format(os.path.join(_ROOT,
        'util/wavefield_from_instaseis.py'),source_model))
    click.secho("Copied default source_config.json and measr_config.json \
to source model directory, please edit. \
Please run setup_noisesource.ipynb or setup_noisesource.py after editing to \
create starting model.")
    

###########################################################################
### Preprocess the sytnthetic wavefields
########################################################################### 
@run.command(help='Filter & truncate synthetics before correlation.')
@click.argument('source_model')
def preprocess_synthetics(source_model):
    source_model = os.path.join(source_model,'source_config.json')
    source_config = json.load(open(source_model))
    if source_config['preprocess_do']:
        
        dir = os.path.join(source_config['source_path'],'wavefield_processed')
        
        try:
            os.mkdir(dir)
        except:
            pass    
            
        run_preprocessing(source_config)


###########################################################################
### Correlations <3
###########################################################################
@run.command(help='Calculate correlations for selected source model.')
@click.argument('source_model')
@click.argument('step')
def correlation(source_model,step):
    source_model = os.path.join(source_model,'source_config.json')
    run_corr(source_model,step)
    

###########################################################################
### Measure and get adjoint sources
###########################################################################
@run.command(help='Run measurement and adjoint sources.')
@click.argument('source_model')
# To do: Include a --test option that produces only plots 
# To do: include a get_parameters_options or something, so that there is no 
# extra step necessary in run_measurement
@click.argument('step')
@click.option('--ignore_network',is_flag=True)
@click.option('--step_test',is_flag=True)
def measurement(source_model,step,ignore_network,step_test):
    
    measr_config = os.path.join(source_model,'measr_config.json')
    source_model = os.path.join(source_model,'source_config.json')
    
    run_measurement(source_model,measr_config,int(step),ignore_network,
        step_test)
    if not step_test:
        run_adjointsrcs(source_model,measr_config,int(step),ignore_network)


###########################################################################
### Get kernels (without residuals multiplied)
###########################################################################
@run.command(help='Calculate preliminary kernels.')
@click.argument('source_model')
@click.argument('step')
@click.option('--ignore_network',is_flag=True)
def kernel(source_model,step,ignore_network):
    source_model = os.path.join(source_model,'source_config.json')
    run_kern(source_model,step,ignore_network=ignore_network)


###########################################################################
### Step length test forward model
###########################################################################
@run.command(help='Calculate fewer correlations for step length test.')
@click.argument('source_model')
@click.argument('step')
def step_test(source_model,step):
    source_model = os.path.join(source_model,'source_config.json')
    run_corr(source_model,step,steplengthrun=True)


###########################################################################
### Assemble the gradient by multplying kernels by residuals and summing
###########################################################################
@run.command(help='Assemble ascent direction from spatial kernels and \
measurements')
@click.argument('source_model')
@click.argument('step')
@click.option('--snr_min',default=0.0)
@click.option('--n_min',default=0)
@click.option('--normalize',default=False)
def gradient(source_model,step,snr_min,n_min,normalize):
    snr_min = float(snr_min)
    source_model = os.path.join(source_model,'source_config.json')
    assemble_ascent_dir(source_model,step,snr_min,
        n_min,normalize_gradient=normalize)
    


###########################################################################
### Older stuff, might be useful again but maybe not
###########################################################################

###########################################################################
### Old: prepare input for specfem
###########################################################################
@run.command(help='Prepare specfem input.')
@click.argument('project_path')
def specfem_input(project_path):
    prepare_specfem_input(os.path.join(project_path,'config.json'))


###########################################################################
### Old: Preprocess data (filtering is done anyway by measurement, if asked!)
########################################################################### 
@run.command(help='Preprocess observed correlations')
@click.argument('source_model')
@click.option('--bandpass',help='Bandpass filter, format: freq1 freq2 corners.',
    default=None)
@click.option('--decimator',help='Decimation factor. Default obspy antialias \
filter will be run before decimating.',default=None)
@click.option('--fs_new',help='New sampling rate. Ensure that filtering is \
performed before!',default=None)
def preprocess_data(source_model,bandpass,decimator,fs_new):

    if bandpass is not None:
        bandpass = [float(bp) for bp in bandpass.split()]

    if fs_new is not None:
        fs_new = float(fs_new)

    if decimator is not None:
        decimator = int(decimator)

    
    run_preprocess_data(source_model,bandpass=bandpass,
        decimator=decimator,Fs_new=fs_new)
