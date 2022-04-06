#! /usr/bin/env python3
import glob,re
from pymol import cmd

def align_allfiles_to_allfiles(files=None,selection='name ca',cutoff=2,cycles=5,debug=0,full_matrix=0,method='align'):
  """
  Aligns all models in a list to all other models in the list

  usage:
    align_allfiles_to_allfiles [file_list][selection][cutoff=2][cycles=5][debug=0][full_matrix=0][method='align']

        where method can be align, super, cealign or rmscur

        where selection, cutoff and cycles are options passed to the align command.
        (e.g. one chain of each protein),

        cutoff and cycles are options passed to the align or super command.

    Setting debug=1 prints more information to the terminal or external GUI.
    Setting full_matrix=1 prints out the full symmetric matrix, rather than
    simply the top-half matrix

    Example:
      align_allfiles_to_allfiles files=name1 name2 name3 name4, selection=c. a & n. ca, full_matrix=1
  """
  file_list = glob.glob(files)
  file_list.sort()
  extension = re.compile( '(^.*[\/]|\.(pdb|ent|brk))' )

  cutoff = int(cutoff)
  cycles = int(cycles)
  full_matrix = int(full_matrix)
  debug=int(debug)

  object_list = []
  for filename in file_list:
    object_list.append(extension.sub('',filename))
  print(object_list)
  rmsd = {}
  rmsd_list = []
  for i in range(len(file_list)):
    obj_name1 = extension.sub('',file_list[i])
    cmd.load(file_list[i],obj_name1)
    for j in range(i+1,len(file_list)):
      obj_name2 = extension.sub('',file_list[j])
      cmd.load(file_list[j],obj_name2)

      if method == 'align':
        rms = cmd.align('%s & %s' % (object_list[j],selection),'%s & %s' % (object_list[i],selection),cutoff=cutoff,cycles=cycles)
      elif method == 'super':
        rms = cmd.super('%s & %s' % (object_list[j],selection),'%s & %s' % (object_list[i],selection),cutoff=cutoff,cycles=cycles)
      elif method == 'cealign':
        rmsdict = cmd.cealign('%s & %s' % (object_list[i],selection),'%s & %s' % (object_list[j],selection))
        rms = [rmsdict['RMSD'],rmsdict['alignment_length'],1,0,0]
      elif method == 'rms_cur':
#        print'mobile: %s & %s' % (object_list[j],selection),'target: %s & %s' % (object_list[i],selection)
        num_atoms = cmd.select('junkselection','%s & %s' % (object_list[j],selection))
        cmd.delete('junkselection')
        rms = [cmd.rms_cur('%s & %s' % (object_list[j],selection),'%s & %s' % (object_list[i],selection)),num_atoms]

      rmsd.setdefault(obj_name1,{})[obj_name2] = rms[0]
      rmsd_list.append((rms[0],obj_name1,obj_name2))
      if debug:
        print("Alignment of %s to %s:" % (obj_name2,obj_name1))
        print("     Initial RMS: %6.3f for %d atoms" % (rms[3],rms[4]))
        print("     Final RMS: %6.3f for %d atoms after %d cycles\n" % (rms[0],rms[1],rms[2]))

      cmd.delete(obj_name2)
    cmd.delete(obj_name1)

  rmsd_list.sort()
# loop over dictionary and print out matrix of final rms values
  if debug:
    for object_name in object_list[:-1]:
      print("%s: %s" % (object_name,str(rmsd[object_name])))

    for r in rmsd_list:
      print("%6.3f  %s  %s" % r)

  print("%6s" % " ", end=' ')
  if full_matrix:
# fill in other half of matrix
    for i in range(len(object_list)):
      for j in range(i+1,len(object_list)):
        rmsd.setdefault(object_list[j],{})[object_list[i]] = rmsd[object_list[i]][object_list[j]]
      rmsd[object_list[i]][object_list[i]] = 0

    for i in range(len(rmsd)):
      print("%6s" % object_list[i], end=' ')
    print("")
    for i in range(len(object_list)):
      print("%6s" % object_list[i], end=' ')
      for j in range(len(object_list)):
        print("%6.3f" % (rmsd[object_list[i]][object_list[j]]), end=' ')
      print("")
  else:
    for i in range(len(rmsd)):
      print("%6s" % object_list[i+1], end=' ')
    print("")
    for i in range(len(object_list)):
      print("%6s" % object_list[i], end=' ')
      for k in range(i):
        print("%6s" % " ", end=' ')
      for j in range(i+1,len(object_list)):
        print("%6.3f" % (rmsd[object_list[i]][object_list[j]]), end=' ')
      print("")

cmd.extend('align_allfiles_to_allfiles',align_allfiles_to_allfiles)
