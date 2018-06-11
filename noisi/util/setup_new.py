import os
import io
import time
import json
from noisi import _ROOT

def setup_proj(project_name):

    os.makedirs(os.path.join(project_name))
    
    with io.open(os.path.join(_ROOT,'config','config.json'),'r+') as fh:
        conf = json.loads(fh.read())
    
    # For gaussian grid
    with io.open(os.path.join(_ROOT,'config','config_gauss.json'),'r+') as fh:
        conf_gauss = json.loads(fh.read())
        
    conf['date_created'] = time.strftime("%Y.%m.%d")
    conf['project_name'] = project_name
    conf['project_path'] = os.path.abspath(project_name)
    
    conf_gauss['date_created'] = time.strftime("%Y.%m.%d")
    conf_gauss['project_name'] = project_name
    conf_gauss['project_path'] = os.path.abspath(project_name)
    
    with io.open(os.path.join(project_name,'config.json'),'w') as fh:
        cf = json.dumps(conf,sort_keys=True, indent=4, separators=(",", ": "))
        fh.write(cf)
        
    # For gaussian grid
    with io.open(os.path.join(project_name,'config_gauss.json'),'w') as fh:
        cf = json.dumps(conf_gauss,sort_keys=True, indent=4, separators=(",", ": "))
        fh.write(cf)
        
    # Copy gaussian grid notebook
    os.system('cp {} {}'.format(os.path.join(_ROOT,'jnotebks/setup_gaussian_grid.ipynb'),
    project_name))