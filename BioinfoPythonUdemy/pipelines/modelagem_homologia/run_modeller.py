#run_modeller.py

#-*- coding:utf-8 -*-
from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()

env.io.atom_files_directory = ["/home/vinibfranc/Documentos/Scripting/BioinfoPythonUdemy/pipelines/modelagem_homologia"]

env.io.hetatm = True
env.io.water = True

a = automodel(env, alnfile = 'alinha.pir', knowns = '4hpg', sequence = 'bgl')
a.starting_model = 1
a.ending_model = 5
a.make()