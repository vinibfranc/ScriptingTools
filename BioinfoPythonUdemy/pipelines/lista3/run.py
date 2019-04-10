#-*- coding:utf-8 -*-
from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()

env.io.atom_files_directory = ['/home/vinibfranc/Documentos/Scripting/BioinfoPythonUdemy/pipelines/lista3']

env.io.hetatm = True
env.io.water = True

a = automodel(env, alnfile = 'new_alinha.pir', knowns = '3UR8', sequence = 'beta_glicosidase.fasta')
a.starting_model = 1
a.ending_model = 2
a.make() 
