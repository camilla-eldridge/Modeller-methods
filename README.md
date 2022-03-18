# Modeller-methods


Step by step of methods to model and evaluate a protein structure based on homologus template (later automated ==> see Modeller-automation repo) 

'Modeller is used for homology or comparative modeling of protein three-dimensional structures.
The user provides an alignment of a sequence to be modeled with known related structures and MODELLER automatically calculates a model containing all non-hydrogen atoms'

for manual see: https://salilab.org/modeller/9.19/manual.pdf





Example work through with tsp23 using the example scripts and methods provided  at https://salilab.org/modeller/tutorial/basic.html 


1. Find homologous structures for templates - I usually use Phyre2 for this...



2. Download pdb files from pdb database.


 
3. Make sure mod9.21 is in $PATH.



4. Put pdb template files and sequence in working dir with modeller scripts found here https://salilab.org/modeller/tutorial/basic.html or copy from below.




5. Run compare.py (if you dont know which the closest template is, use this:




        /usr/bin/python2.7 compare.py > compare.log
       
       
        compare.py: 
      
        from modeller import *

        env = environ()
        aln = alignment(env)
        for (pdb, chain) in (('1b8p', 'A'), ('1bdm', 'A'), ('1civ', 'A'),
                             ('5mdh', 'A'), ('7mdh', 'A'), ('1smk', 'A')):
            m = model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
            aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)
        aln.malign()
        aln.malign3d()
        aln.compare_structures()
        aln.id_table(matrix_file='family.mat')
        env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)
      

        grep -A 50 "Weighted pair-group" compare.log 
        
        
 Look for structure in dendrogram with: 
- Greatest sequence similarity
- Crystallographic resolution (@n at tip of branch on tree), the larger the better.



For tsp23 we already knew that the only close templates at the time were 5TCX and 2M7Z so this step was skipped.




7.Align the sequence to template.


         /usr/bin/python2.7 align_2D_single.py > align_2D_single.log


         align_2D_single.py:

         from modeller import *

         env = environ()
         aln = alignment(env)
         mdl = model(env, file='5TCX', model_segment=('FIRST:A','LAST:A'))
         aln.append_model(mdl, align_codes='5TCX', atom_files='5TCX.pdb')
         aln.append(file='sttsp23.ali', align_codes='sttsp23')
         aln.align2d()
         aln.write(file='sttsp23-5tcx.ali', alignment_format='PIR')
         aln.write(file='sttsp23-5tcx.pap', alignment_format='PAP')
      
      
      
      
8.Create a model based on the alignment.


           /usr/bin/python2.7 model-single.py > model-single.log


           model-single.py:

           from modeller import *
           from modeller.automodel import *
           #from modeller import soap_protein_od

           env = environ()
           a = automodel(env, alnfile='sttsp23-5tcx.ali',
                         knowns='5TCX', sequence='sttsp23',
                         assess_methods=(assess.DOPE,
                                         #soap_protein_od.Scorer(),
                                         assess.GA341))
           a.starting_model = 1
           a.ending_model = 5
           a.make()

               
         
         
9.Check dope scores for models (lowest == best)


          grep -A 20 ">> Summary"  model-single.log
          

          Summary of successfully produced models:
           Filename                      molpdf     DOPE score    GA341 score
          ----------------------------------------------------------------------
          sttsp23.B99990001.pdb          980.14014   -22154.52930        0.19803
          sttsp23.B99990002.pdb          967.76953   -21628.17188        0.16171
          sttsp23.B99990003.pdb          894.84656   -21737.41211        0.09677
          sttsp23.B99990004.pdb          889.46710   -21908.58203        0.04608
          sttsp23.B99990005.pdb          948.67993   -21889.07617        0.09870

        
             
        
        
10. Evaluate the template by viewing the energy profile and looking at the DOPE scores (model and template should align).


          /usr/bin/python2.7 evaluate_template.py > evaluate_template.log


          evaluate_template.py:

          from modeller import *
          from modeller.scripts import complete_pdb

          log.verbose()    # request verbose output
          env = environ()
          env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
          env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

          # directories for input atom files
          env.io.atom_files_directory = './:../atom_files'

          # read model file
          mdl = complete_pdb(env, '5TCX.pdb', model_segment=('FIRST:A', 'LAST:A'))

          s = selection(mdl)
          s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='5tcx.profile',
                        normalize_profile=True, smoothing_window=15)

        
   
   
   ![image](https://user-images.githubusercontent.com/12966869/159078417-3e78a32f-85bd-4c53-a633-610ed59616ec.png)

   
   
   
 
        
        
**For a multiple model template:**

1.Align 3D structures to each other using salign()


           /usr/bin/python2.7 salign_tsp.py > salign_tsp.log


           salign_tsp.py:

           from modeller import *

           log.verbose()
           env = environ()
           env.io.atom_files_directory = './:../atom_files/'

           aln = alignment(env)
           for (code, chain) in (('2M7Z', 'A'), ('5TCX', 'A')):
               mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain))
               aln.append_model(mdl, atom_files=code, align_codes=code+chain)

           for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                               ((1., 0.5, 1., 1., 1., 0.), False, True),
                                               ((1., 1., 1., 1., 1., 0.), True, False)):
               aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                          rr_file='$(LIB)/as1.sim.mat', overhang=30,
                          gap_penalties_1d=(-450, -50),
                          gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
                          dendrogram_file='tsp.tree',
                          alignment_type='tree', # If 'progresive', the tree is not
                                                 # computed and all structues will be
                                                 # aligned sequentially to the first
                          feature_weights=weights, # For a multiple sequence alignment only
                                                   # the first feature needs to be non-zero
                          improve_alignment=True, fit=True, write_fit=write_fit,
                          write_whole_pdb=whole, output='ALIGNMENT QUALITY')

           aln.write(file='tsp.pap', alignment_format='PAP')
           aln.write(file='tsp.ali', alignment_format='PIR')

           aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
                      rr_file='$(LIB)/as1.sim.mat', overhang=30,
                      gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
                      gap_gap_score=0, gap_residue_score=0, dendrogram_file='tsp23.tree',
                      alignment_type='progressive', feature_weights=[0]*6,
                      improve_alignment=False, fit=False, write_fit=True,
                      write_whole_pdb=False, output='QUALITY')        
        
        
        
        
        
2.Align sequence to 3D structure alignment 


            /usr/bin/python2.7 align_2dmult.py > multalign_tsp.23.log


             align_2dmult.py:

             from modeller import *

             log.verbose()
             env = environ()

             env.libs.topology.read(file='$(LIB)/top_heav.lib')

             # Read aligned structure(s):
             aln = alignment(env)
             aln.append(file='tsp.ali', align_codes='all')
             aln_block = len(aln)

             # Read aligned sequence(s):
             aln.append(file='sttsp23.ali', align_codes='sttsp23')

             # Structure sensitive variable gap penalty sequence-sequence alignment:
             aln.salign(output='', max_gap_length=20,
                        gap_function=True,   # to use structure-dependent gap penalty
                        alignment_type='PAIRWISE', align_block=aln_block,
                        feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
                        gap_penalties_1d=(-450, 0),
                        gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
                        similarity_flag=True)

             aln.write(file='sttsp23-mult.ali', alignment_format='PIR')
             aln.write(file='sttsp23-mult.pap', alignment_format='PAP')         
         
         
         
    
         
         
3. Model based on alignment to template.


            /usr/bin/python2.7 model-mult_dope.py > model_tsp23.log


            model-mult_dope.py:

            from modeller import *
            from modeller.automodel import *

            env = environ()
            a = automodel(env, alnfile='sttsp23-mult.ali',
                          knowns=('2M7ZA','5TCXA'), sequence='sttsp23',
                    assess_methods=(assess.DOPE,
                    #soap_protein_od.Scorer(),
                                assess.GA341))
            a.starting_model = 1
            a.ending_model = 5
            a.make()

        
4. Choose best model based on lowest DOPE score.    

           grep -A 20 ">> Summary"  model_tsp23.log 


           Summary of successfully produced models:
           Filename                      molpdf          DOPE score       GA341 score
           ----------------------------------------------------------------------
           sttsp23.B99990001.pdb         3085.58496   -22093.03906        0.11341
           sttsp23.B99990002.pdb         3216.13330   -21675.56445        0.09355
           sttsp23.B99990003.pdb         3107.46289   -21665.86719        0.11687
           sttsp23.B99990004.pdb         3018.95239   -21871.83789        0.14612
           sttsp23.B99990005.pdb         3223.13135   -21736.32422        0.09807


From modeller manual:
"The best model was chosen based on DOPE score, you could pick the model with the lowest value of the MODELLER objective function or the DOPE or SOAP assessment scores, or with the highest GA341 assessment score, which are reported at the end of the log file, above. GA341 scores always range from 0.0 (worst) to 1.0 (native-like); however GA341 is not as good as DOPE or SOAP at distinguishing 'good' models from 'bad' models."


5. Rename best model
           mv sttsp23.B99990001.pdb STTSP23_mult.pdb       

6. Evaluate model
           /usr/bin/python2.7 evaluate_model.py > evaluate_model.log


7.Take profile from single template (sttsp23.profile and plot with new profile from multiple templates sttsp23_mult.profile)


             /usr/bin/python2.7 plot_profiles_multi.py  


             plot_profiles_multi.py:
             import pylab
             import modeller

             def r_enumerate(seq):
                 """Enumerate a sequence in reverse order"""
                 # Note that we don't use reversed() since Python 2.3 doesn't have it
                 num = len(seq) - 1
                 while num >= 0:
                     yield num, seq[num]
                     num -= 1

             def get_profile(profile_file, seq):
                 """Read `profile_file` into a Python array, and add gaps corresponding to
                    the alignment sequence `seq`."""
                 # Read all non-comment and non-blank lines from the file:
                 f = open(profile_file)
                 vals = []
                 for line in f:
                     if not line.startswith('#') and len(line) > 10:
                         spl = line.split()
                         vals.append(float(spl[-1]))
                 # Insert gaps into the profile corresponding to those in seq:
                 for n, res in r_enumerate(seq.residues):
                     for gap in range(res.get_leading_gaps()):
                         vals.insert(n, None)
                 # Add a gap at position '0', so that we effectively count from 1:
                 vals.insert(0, None)
                 return vals

             e = modeller.environ()
             a = modeller.alignment(e, file='sttsp23-mult.ali')


             model = get_profile('sttsp23_mult.profile', a['sttsp23'])
             template = get_profile('sttsp23.profile', a['sttsp23'])


             # Plot the template and model profiles in the same plot for comparison:
             pylab.figure(1, figsize=(10,6))
             pylab.xlabel('Alignment position')
             pylab.ylabel('DOPE per-residue score')
             pylab.plot(model, color='red', linewidth=2, label='Multiple templates')
             pylab.plot(template, color='green', linewidth=2, label='Basic model')
             pylab.legend()
             pylab.savefig('dope_profile.png', dpi=65)     

        
      
  ![image](https://user-images.githubusercontent.com/12966869/159078470-51e2bfe7-f443-4d16-8d95-6349f8b9b574.png)
      
          
        
        
8. Loop refinement

Attempted loop refinement between 120:190 

           /usr/bin/python2.7 loop_refine.py > loop_refine.log


           loop_refine.py:

           from modeller import *
           from modeller.automodel import *

           log.verbose()
           env = environ()

           # directories for input atom files
           env.io.atom_files_directory = './:../atom_files'

           # Create a new class based on 'loopmodel' so that we can redefine
           # select_loop_atoms (necessary)
           class MyLoop(loopmodel):
               # This routine picks the residues to be refined by loop modeling
               def select_loop_atoms(self):
                   # 10 residue insertion
                   return selection(self.residue_range('120', '190'))

           m = MyLoop(env,
                      inimodel='STTSP23_mult.pdb', # initial model of the target
                      sequence='sttsp23')          # code of the target

           m.loop.starting_model= 1           # index of the first loop model 
           m.loop.ending_model  = 10          # index of the last loop model
           m.loop.md_level = refine.very_fast # loop refinement method; this yields
                                              # models quickly but of low quality;
                                              # use refine.slow for better models



9.Model energies

        
         from modeller import *
         from modeller.scripts import complete_pdb

         log.verbose()    # request verbose output
         env = environ()
         env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
         env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

         for i in range(1, 5):
             # read model file
             code = "sttsp23.BL%04d0001.pdb" % i
             mdl = complete_pdb(env, code)
             s = selection(mdl)
             s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='sttsp23_mult.profile',
                           normalize_profile=True, smoothing_window=15)       

        
        
10. Plot all dope scores



          Plotdope_scores.py

          import modeller

          k = modeller.environ()
          y = modeller.alignment(k, file='sttsp23-mult.ali')

          model = get_profile('sttsp23_multi_model.profile', y['sttsp23'])  #   multi model 
          template = get_profile('sttsp23_single_model.profile', y['sttsp23'])    #  single model
          loop_refined = get_profile('sttsp23_multi_model_loop_refined.profile', y['sttsp23'])


          #* Plot the template and model profiles in the same plot for comparison *#

          pylab.figure(1, figsize=(10,6))
          pylab.xlabel('Alignment position')
          pylab.ylabel('DOPE per-residue score')

          pylab.plot(model, color='red', linewidth=2, label='Multiple templates:' + '5TCX' + " " +  '2M7Z')
          pylab.plot(template, color='green', linewidth=2, label='Basic model:' + '5TCX')
          pylab.plot(loop_refined, color='blue', linewidth=2, label='Loop refinement(DOPE)')

          pylab.legend()
          pylab.savefig(sttsp23 + 'dope_profile.png', dpi=65)

     
        
        
       ![image](https://user-images.githubusercontent.com/12966869/159079481-383752ca-c785-4d7b-b460-8390de6586f4.png)
 
        
        
        
        I interpreted this as the model with loop refinement being no better than the basic template model....
