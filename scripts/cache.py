#!/usr/bin/env python
import os,sys
import subprocess
import yaml

def container(config):
    print('Checking container cache status...')
    cachedir = config['cachedir']
    containerdir = '%s/container'%(cachedir)
    if not os.path.exists(containerdir):
        subprocess.check_call(['mkdir', '-p', containerdir])
    
    # Download builded images
    for img in config['simg']:
        if not os.path.exists('%s/%s.simg.ok'%(containerdir, img)):
            pass
        else:
            with open('%s/%s.simg.ok'%(containerdir, img), 'r') as infile:
                if infile.readline().strip() == config['simg'][img]:
                    print('Container', img, 'exists, skip...')
                    continue
                else:
                    pass
        try:
            subprocess.check_call(['wget', config['simg'][img], '-O', '%s/%s.simg'%(containerdir, img)])
        except subprocess.CalledProcessError:
            raise Exception("Command exited with non-zero code.")
        with open('%s/%s.simg.ok'%(containerdir, img), 'w') as savefile:
            savefile.write(config['simg'][img]+'\n')
            print('Container', img, 'downloaded...')
    
    # Pull from docker registries
    for registry in ['docker','singularity']:
        for img in config[registry]:
            if not os.path.exists('%s/%s.simg.ok'%(containerdir, img)):
                pass
            else:
                with open('%s/%s.simg.ok'%(containerdir, img), 'r') as infile:
                    if infile.readline().strip() == config[registry][img]:
                        print('Container', img, 'exists, skip...')
                        continue
                    else:
                        pass
            try:
                cmd = 'module load singularity && singularity pull %s/%s.simg %s'%(containerdir, img, config[registry][img])
                print(cmd)
                subprocess.check_call([cmd], shell=True, executable="/bin/bash")
            except subprocess.CalledProcessError:
                raise Exception("Command exited with non-zero code.")
            with open('%s/%s.simg.ok'%(containerdir, img), 'w') as savefile:
                savefile.write(config[registry][img]+'\n')
                print('Container', img, 'pulled...')
    
    print('Container cache finished!')
    print('#'*80)
    return True

def reference(config):
    print('Checking reference cache status...')
    cachedir = config['cachedir']
    refdir = '%s/reference'%(cachedir)
    
    if not os.path.exists(refdir):
        subprocess.check_call(['mkdir', '-p', refdir])
    elif not os.path.exists('%s/ref.ok'%(refdir)):
        pass
    else:
        with open('%s/ref.ok'%(refdir)) as infile:
            refexists = [i.strip() for i in infile.readlines()]
        for item in config['references']:
            if config['references'][item]['link'] in refexists:
                print('Reference', item, 'exists, skip...')
                continue
            else:
                print(item, 'in refdir does not match current config:', config['references'][item]['link'])
                print('Redo the reference cache...')
                subprocess.check_call(['rm', '-rf', refdir])
                subprocess.check_call(['mkdir', '-p', refdir])
                
    # Download reference packages
    if os.path.exists('%s/download.ok'%(refdir)):
        with open('%s/download.ok'%(refdir)) as infile:
            downloaded = [i.strip() for i in infile.readlines()]
    else:
        downloaded = []
        
    for item in config['references']:
        if config['references'][item]['link'] in downloaded:
            continue
        try:
            if 'pack' in config['references'][item]:
                cmd = 'cd %s && wget %s'%(refdir, config['references'][item]['link'])
            else:
                if 'http' in config['references'][item]['link']:
                    cmd = 'cd %s && wget %s -O %s/%s'%(refdir, config['references'][item]['link'], cachedir, config['references'][item]['path'])
                else:
                    cmd = 'cd %s && ln -s %s %s/%s'%(refdir, config['references'][item]['link'], cachedir, config['references'][item]['path'])
            subprocess.check_call(cmd, shell=True)
            with open('%s/download.ok'%(refdir), 'a') as savefile:
                savefile.write(config['references'][item]['link'] + '\n')
        except subprocess.CalledProcessError:
            raise Exception("Command exited with non-zero code.")
    print('Reference download finished!')
    
    # Unpack and organize files
    if os.path.exists('%s/unpack.ok'%(refdir)):
        with open('%s/unpack.ok'%(refdir)) as infile:
            unpacked = [i.strip() for i in infile.readlines()]
    else:
        unpacked = []
    
    ## hmftools
    item = 'hmftools'
    if config['references'][item]['link'] not in unpacked:
        cmd  = 'cd %s && '%(refdir)
        cmd += 'tar zxvf %s && '%(config['references'][item]['link'].split('/')[-1])
        cmd += 'mv %s/%s %s/%s '%(cachedir, config['references'][item]['pack'], cachedir, config['references'][item]['path'])
        for clean in config['references'][item]['clean']:
            cmd += '&& rm -rf %s/%s'%(cachedir, clean)
        try:
            subprocess.check_call([cmd], shell=True)
            with open('%s/unpack.ok'%(refdir), 'a') as savefile:
                savefile.write(config['references'][item]['link']+'\n')
        except subprocess.CalledProcessError:
            raise Exception("Command exited with non-zero code.")
    print('Reference unpack finished!')
    
    with open('%s/ref.ok'%(refdir), 'w') as savefile:
        for item in config['references']:
            savefile.write(config['references'][item]['link'] + '\n')
    print('Reference cache finished!')
    print('#'*80)
    return True
    
def test_container(configpath):
    with open(configpath, 'r') as infile:
        config = yaml.safe_load(infile)
    config['pipelinedir'] = '/'.join(os.path.dirname(__file__).split('/')[:-1])
    if config['cachedir'][0] != '/':
        config['cachedir'] = '%s/%s'%(config['pipelinedir'], config['cachedir'])
    container(config)
    return True
    
def test_reference(configpath):
    with open(configpath, 'r') as infile:
        config = yaml.safe_load(infile)
    config['pipelinedir'] = '/'.join(os.path.dirname(__file__).split('/')[:-1])
    if config['cachedir'][0] != '/':
        config['cachedir'] = '%s/%s'%(config['pipelinedir'], config['cachedir'])
    reference(config)
    return True

def cleancache(configpath):
    with open(configpath, 'r') as infile:
        config = yaml.safe_load(infile)
    config['pipelinedir'] = '/'.join(os.path.dirname(__file__).split('/')[:-1])
    if config['cachedir'][0] != '/':
        config['cachedir'] = '%s/%s'%(config['pipelinedir'], config['cachedir'])
    subprocess.check_call(['rm -rf', config['cachedir']], shell=True, executable="/bin/bash")
    
    
if __name__ == '__main__':
    test_c = test_container('%s/config/container.yaml'%('/'.join(os.path.dirname(__file__).split('/')[:-1])))
    if not test_c:
        raise Exception("Container function test failed.")
    test_r = test_reference('%s/config/reference.yaml'%('/'.join(os.path.dirname(__file__).split('/')[:-1])))
    if not test_r:
        raise Exception("Reference function test failed.")
    cleancache('%s/config/container.yaml'%('/'.join(os.path.dirname(__file__).split('/')[:-1])))
        
    