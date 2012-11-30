## @package sbml2rdme
#
#
#  The sbml2rdme class definitions
#
#  @author V.Gerdin 

import sys
import math

## Model class
# The complete model built by the parser. A model's collection of entities are connected to each other
#
class model():
  ## Constructor
  # @param name Optional name of the model. Default is "model"
  def __init__(self, name="model"):
    self.setName(name)
    self.globals = {}        #empty dictionary
    self.reactions = {}      #empty dictionary
    self.reactions_order = [] #empty list
    self.species = {}        #empty dictionary
    self.speciesTypes = {}   #empty dictionary
    self.speciesTypes_order = [] #empty list
    self.compartmentTypes = {} #empty dictionary
    self.compartmentTypes_order = [] #empty list
    self.compartments = {}     #empty dictionary
    self.data = {} #empty list
  
  ## getName
  # @return Returns name of the model
  def getName(self):
    return self.name

  ## setName
  # @param name New name of model
  # Will check that name is propper and if so set it. Otherwise set "model"
  def setName(self, name):
    if name and str(name).strip():
      self.name = str(name).strip()
    else:
      model_warning("Attempt to put model name as empty string or None. Will use \"model\" instead.")
      self.name = "model"

##
# COMPARTMENT TYPES
##

  ## addCompartmentType
  # @param ctype CompartmentType object to add
  # Add CompartmentType object to model. If there already is a CompartmentType with the same id it will be replace
  # with a warning. Inpropper input will raise SystemExit.
  def addCompartmentType(self, ctype):
    if not isinstance(ctype, CompartmentType):
      model_error("Attempting to add non compartment type (" + ctype.getId() + ")  to model.")
      #will exit(1)
    if self.getCompartmentType(ctype.getId()):
      model_warning("Reassigning compartment type " + ctype.getId())
    else:
      self.compartmentTypes_order.append(ctype.getId())
    self.compartmentTypes[ctype.getId()] = ctype

  ## getCompartmentType
  # @param sid id of CompartmentType
  # @return CompartmentType object with sid or None if not found.
  def getCompartmentType(self, sid):
    if sid:
      return self.compartmentTypes.get(sid)
    else:
      return None
  
  ## getNextSD 
  # Fetches what the next sub domain should be labled. These are assigned in
  # cronological order of compartment type creation.
  # @return integer of next sub domain index
  def getNextSD(self):
    return len(self.compartmentTypes)+1

##
#  SPECIES TYPE
##
 
  ## getSpeciesType
  # @param sid id of speciesType
  # @return speciesType object with sid or None if not found.
  def getSpeciesType(self, sid):
    if sid:
      return self.speciesTypes.get(sid)
    else:
      return None

  ## addSpeciesType
  # @param stype SpeciesType object to add
  # Add SpeciesType object to model. If there already is a SpeciesType with the same id it will be replace
  # with a warning. Inpropper input will raise SystemExit.
  def addSpeciesType(self, stype):
    if not isinstance(stype, SpeciesType):
      model_error("Attempting to add non species type (" + stype.sid + ")  to model.")
      #will exit(1)
    if self.getSpeciesType(stype.sid):
      model_warning("Reassigning species type " + stype.sid)
    else:
      self.speciesTypes_order.append(stype.sid)
    self.speciesTypes[stype.sid] = stype
  
##
#  COMPARTMENT
##

  ## addCompartment
  # @param comp Compartment object to add
  # Add Compartment object to model. If there already is a Compartment with the same id it will be replace
  # with a warning. Inpropper input will raise SystemExit.
  def addCompartment(self, comp):
    if not isinstance(comp, Compartment):
      model_error("Attempting to add non compartment (" + comp.getId() + ")  to model.")
      #will exit(1)
    if self.getCompartment(comp.getId()):
      model_warning("Reassigning compartment " + comp.getId())
    self.compartments[comp.getId()] = comp
    comp.getType().addMember(comp)

  ## getCompartment
  # @param cid id of Compartment
  # @return Compartment object with cid or None if not found.
  def getCompartment(self, cid):
    if cid:
      return self.compartments.get(cid)
    else:
      return None

##
#  SPECIES
##

  def addSpecies(self, spec):
    if not isinstance(spec, Species):
      model_error("Attempting to add non speices (" + spec.getId() + ")  to model.")
      #will exit(1)
    if self.getSpecies(spec.getId()):
      model_warning("Reassigning species " + spec.getId())
    self.species[spec.getId()] = spec
    spec.getType().addMember(spec)

  def getSpecies(self, sid):
    if sid and sid in self.species:
      return self.species[sid]
    else:
      return None

##
#  PARAMETERS
##
  def addGlobal(self, param):
    if not isinstance(param, Parameter):
      model_error("Attempting to add non Parameter (global) (" + param.getId() + ") to model.")
      #will exit(1)
    if self.getGlobal(param.getId()):
      model.warning("Reassigning global parameter " + param.getId())
    self.globals[param.getId()] = param

  def getGlobal(self, pid):
    if pid:
      return self.globals.get(pid)
    else:
      return None

  def addData(self, param):
    if not isinstance(param, Parameter):
      model_error("Attempting to add non Parameter (data) (" + param.getId() + ") to model.")
      #will exit(1)
    if not self.getData(param.getId()):
      param.data = len(self.data)
      self.data[param.getId()] = param

  def getData(self, pid):
    if pid and pid in self.data:
      return self.data.get(pid)
    else:
      return None

##
#  REACTION
##

  def addReaction(self, reaction):
    if not isinstance(reaction, Reaction):
      model_error("Attempting to add non Reaction (" + reaction.getId() + ") to model.")
      #will exit(1)
    if self.getReaction(reaction.getId()):
      model.warning("Reassigning Reaction " + reaction.getId())
    else:
      self.reactions_order.append(reaction.getId())
    self.reactions[reaction.getId()] = reaction

  def getReaction(self, rid):
    if rid and rid in self.reactions:
      return self.reactions.get(rid)
    else:
      return None

##
#  GENERATE FILES
##

  ## Generate c-file
  # @param outdir Directory for output 
  # Creates the model.c file from the current model object. Output file is named after the model.
  def c(self, outdir="./"):
    # it is not certain that outdir will end with a '/' which is neccessary for concat  
    if outdir[len(outdir)-1] != '/':
      outdir = outdir + '/'

    filename = self.getName() + ".c"
    cfile = model_openFile(outdir + filename, 'w')
    sys.stderr.write("Creating model c-file " + filename + "\n") 

    # HEADER INFORMATION
    cfile.write("/**\n")
    cfile.write(" * " + filename + "\n")
    cfile.write(" * This is an autogenerated model file, generated by sbml2rdme with intended use for\n")
    cfile.write(" * URDME 1.2 or later using Comsol mutliphysics 4.2 or later.\n")
    cfile.write(" * Please be sure to validate the generated code, especially if sbml2rdme produced any warnings.\n")
    cfile.write("**/\n\n")

    cfile.write("#include <stdlib.h>\n")
    cfile.write("#include \"propensities.h\"\n")
    cfile.write("\n")

    species = self.speciesTypes_order

    cfile.write("\n/* SPECIES */\n")

    spec_i = 0 #species are array indexes thus start at 0
    cfile.write("enum {")
    for spec in species:
      cfile.write(str(spec))
      if spec_i == 0:
        cfile.write('=' + str(spec_i))
      if spec_i < len(species)-1: #we start at 0 so the last element is at len-1
        cfile.write(',')
      spec_i = spec_i + 1
    cfile.write("};\n\n")

    compartments = self.compartmentTypes_order
    nc = len(compartments)

    cfile.write("\n/* COMPARTMENTS */\n")

    if nc > 1 :
      cfile.write("/*\n  Compartments or subdomains are distinguished by their sd value defined below.\n")
      cfile.write("  Make sure that these values correspond to the values and domains in your mesh model.\n*/\n")
      sd_i = 1 #sd is a consequential list starting at 1
      cfile.write("enum {")
      for comp in compartments:
        cfile.write(str(comp))
        if sd_i == 1:
          cfile.write("="+str(sd_i))
          if sd_i < len(compartments): #we start at 1 so the last element is at len
            cfile.write(",")
        sd_i = sd_i+1
      cfile.write("};\n\n")
    else:
      cfile.write("/*\n  The given model does not contains enough different compartments (" + str(nc) + " are detected)\n")
      cfile.write("  for it to be neccessary to add subdomain support.\n*/\n")

    reactions = self.reactions_order
    nr = len(reactions)

    cfile.write("const double NR = " + str(nr) + "; /* Number of reactions */\n")

    #for i in range(0,nr):
    #  re_params = reactions.get(i).getKineticLaw().getListOfParameters();
    #  for j in range(0,re_params.size()):
    #    param = re_params.get(i)
    #    if param.isSetValue():
    #      global_parameters.append(param)

    for param in self.globals:
      cfile.write("const double " + str(param) + ' = ' + str(self.globals[param].getValue()) + ";\n")

    cfile.write("\n")

    if nr > 0:
      # print propensity functions if there are any reactions. If no reactions excist it is a pointless model and warnings
      # are placed here instead.
      for i in range(0,nr):
        pfun = self.getReaction(reactions[i])
        cfile.write("/* Propensity function " + str(i+1) + "*/\n")
        cfile.write("double " + str(pfun.getId()) + "(const int *x, double t, double vol, const double *data, int sd)\n{\n")
        for param_key in pfun.locals:
          param = pfun.getLocal(param_key)
          if param and param.isSet():
            cfile.write("  const double " + str(param.getId()) + " = " + str(param.getValue()) + ";\n" )

        comp = pfun.getCompartment()
        if comp and nc > 1:
          cfile.write("  if(sd == " + str(comp.getId()) + "){\n  ")
        cfile.write("  return " + pfun.getKinetic() + ";\n")
        if comp and nc > 1:
          cfile.write("  }\n  return 0.0;\n")
        cfile.write("}\n\n")

      cfile.write("\n/* COLLECTION OF PROPENSITIES (STATIC CODE) */\n")

      cfile.write("PropensityFun *ALLOC_propensities(void)\n")
      cfile.write("{\n")
      cfile.write("  PropensityFun *ptr = malloc(sizeof(PropensityFun)*NR);\n")
      cfile.write("\n")
      for i in range(0,nr):
        cfile.write("  ptr[" + str(i) + "]=" + str(reactions[i]) + ";\n")
      cfile.write("\n")
      cfile.write("  return ptr;\n")
      cfile.write("}\n")
      cfile.write("\n")
      cfile.write("void FREE_propensities(PropensityFun* ptr)\n")
      cfile.write("{\n")
      cfile.write("  free(ptr);\n")
      cfile.write("}\n")
    else:
      cfile.write("/* WARNING:\n   The model contains no reactions and thus has no propensities. This is a useless model and will not work with URDME.\n   Please ad reactions/propensities to your model before running it.\n*/")



  ## Generate m-file
  # @param outdir Directory for output 
  # Creates the model.m file from the current model object. Output file is named after the model.
  def m(self, outdir="./"):
    # it is not certain that output will end with a '/' which is neccessary for concat  
    if outdir[len(outdir)-1] != '/':
      outdir = outdir + '/'

    filename = self.getName() + ".m"
    mfile = model_openFile(outdir + filename, 'w')
    sys.stderr.write("Creating model m-file " + filename + "\n")
    
    # HEADER INFORMATION
    mfile.write("function umod = " + self.getName() + "(umod)\n")
    mfile.write("%%\n")
    mfile.write("% " + filename + "\n")
    mfile.write("% This is an autogenerated model file, generated by sbml2rdme with intended use for\n")
    mfile.write("% URDME 1.2 or later using Comsol mutliphysics 4.2 or later.\n")
    mfile.write("% Please be sure to validate the generated code, especially if sbml2rdme produced any warnings.\n")
    mfile.write("%\n")
    mfile.write("% After this file has been run the following fields must be assigned in the umod structure\n")
    mfile.write("% for the urdme core to function. Some are assigned by comsol2urdme and the rest are assigned by\n")
    mfile.write("% this file. You are free to alter any preassigned Matrixes if your model requires it but do so\n")
    mfile.write("% only if you know what you are doing.\n")
    mfile.write("%\n")
    mfile.write("% umod.D - Diffusion Matrix*\n")
    mfile.write("% umod.N - Stoichiometry Matrix\n")
    mfile.write("% umod.G - Dependancy Matrix\n")
    mfile.write("% umod.u0 - Initial Concentrations Vector\n")
    mfile.write("% umod.sd - Sub Domain Divison Vector*\n")
    mfile.write("% umod.tspan - Time Vector\n")
    mfile.write("% umod.data - Data Vector\n")
    mfile.write("% umod.vol - Voxel Volume Vector*\n")

    mfile.write("%%\n\n")

    # DIMENSIONS
    mfile.write("\n% DIMENSIONS\n")
    
    stypes = self.speciesTypes_order
    reactions = self.reactions_order
    ctypes = self.compartmentTypes_order
    mspecies = len(stypes)
    nreactions = len(reactions)
    ndomains = len(ctypes)

    mfile.write("Mspecies = " + str(mspecies) + ";\n")
    mfile.write("Nreactions = " + str(nreactions) +";\n")
    mfile.write("Ndomains = " + str(ndomains) + ";\n")
    mfile.write("Ncells = umod.Ncells;\n")
    mfile.write("Ndofs = Mspecies*Ncells;\n")

    mfile.write("species_order = {")

    for i in range(0,mspecies):
      spec = self.getSpeciesType(stypes[i])
      if spec and spec.name:
        mfile.write("\'" + spec.name + "\'")
      else:
        mfile.write("\'" + stypes[i] + "\'")
      if i < mspecies-1:
        mfile.write(",")
    mfile.write("};\n")

    mfile.write("reactions_order = {")

    for i in range(0,nreactions):
      mfile.write("\'" + reactions[i] + "\'")
      if i < nreactions-1:
        mfile.write(",")
    mfile.write("};\n")

    mfile.write("domains_order = {")
    for i in range(0,ndomains):
      mfile.write("\'" + ctypes[i] + "\'")
      if i < ndomains-1:
        mfile.write(",")
    mfile.write("};\n")

    mfile.write("\n% There is a risk that our species' order does not match the one given in Comsol.\n")
    mfile.write("perm = [1:Mspecies];\n")
    mfile.write("if isfield(umod, 'species')\n")
    mfile.write("  [foo,perm] = ismember(species_order, umod.species);\n")
    mfile.write("end\n")

    # If the model is so large that the Stoichiometry or Dependancy matrix gets
    # dimensions larger than 20, then we will output them in vector form instead
    # of a full matrix.
    if mspecies > 20 or nreactions > 20 or mspecies+nreactions > 20:
      sparse_format_full = False
    else:
      sparse_format_full = True

    # STOICHIOMETRY MATRIX
    mfile.write("\n% STOICHIOMETRY MATRIX\n")
    if sparse_format_full:
      mfile.write("% Using full format\n")
      mfile.write("% Every column corresponds to a reaction.\n\n")
      mfile.write("umod.N=sparse([")    
      for spec_key in stypes:
        st_str = ""
        for reaction_key in reactions:
          st_str = st_str + " "
          st = self.reactions[reaction_key].getStoichiometry(spec_key)
          if st >= 0:
            st_str = st_str + " "
          st_str = st_str + str(st)
        mfile.write(st_str + ";\n                    ")
      mfile.write("]);\n")
    else:
      mfile.write("% Using vector format\n")
      mfile.write("% st_i is a vector of species referenses (rows)\n")
      mfile.write("% st_j is a vector of reaction referenses (columns)\n")
      mfile.write("% st_s is a vector of stoichiometries (values)\n\n")
      i_str = "st_i = ["
      j_str = "st_j = ["
      s_str = "st_s = [" 
      for i in range(0,mspecies):
        for j in range(0,nreactions):
          st = self.reactions[reactions[j]].getStoichiometry(stypes[i])
          if st != 0:
            i_str = i_str + str(i+1) + ','
            j_str = j_str + str(j+1) + ','
            s_str = s_str + str(st) + ','
      mfile.write(i_str[:-1] + "];\n")
      mfile.write(j_str[:-1] + "];\n")
      mfile.write(s_str[:-1] + "];\n")
      mfile.write("st_i = perm(st_i); % permute the species references to match the comsol model.\n")
      mfile.write("umod.N = sparse(st_i, st_j, st_s);\n")

    # DEPENDENCY MATRIX
    mfile.write("\n% DEPENDENCY MATRIX\n")

    if sparse_format_full:
      mfile.write("% Using full format:\n")
      mfile.write("% The first Mspecies columns tells which propensities\n")
      mfile.write("% needs to be updated that species diffuses. The following Mreactions\n")
      mfile.write("% columns does the same thing but for reaction events.\n\n")
      mfile.write("umod.G=sparse([")
      for row in reactions:
        dep_str = ""
        for spec_key in stypes:
          if self.speciesTypes[spec_key].isDependant(self.reactions[row]):
            dep_str = dep_str + " 1"
          else:
            dep_str = dep_str + " 0"
        for reaction_key in reactions:
          if self.reactions[reaction_key].isDependant(self.reactions[row]):
            dep_str = dep_str + " 1"
          else:
            dep_str = dep_str + " 0"
        mfile.write(dep_str + ";\n                    ")
      mfile.write("]);\n\n") 
    else:
      mfile.write("% Using vector format:\n")
      mfile.write("% de_i is a vector of reaction referenses (rows)\n")
      mfile.write("% de_j is a vector of both species and reaction referenses (columns)\n")
      mfile.write("%   The first Mspecies (" + str(mspecies) + ") columns tells which propensities\n")
      mfile.write("%   needs to be updated when the amount of a species changes.\n")
      mfile.write("%   The following Nreactions (" + str(nreactions) + ") columns does the same thing\n")
      mfile.write("%   but for reaction events.\n\n")
      i_str = "de_i = ["
      j_str = "de_j = ["
      for j in range(0,mspecies):
        for i in range(0,nreactions):
          if self.speciesTypes[stypes[j]].isDependant(self.reactions[reactions[i]]):
            i_str = i_str + str(i+1) + ','
            j_str = j_str + str(j+1) + ','
      mfile.write(j_str[:-1] + "];\n")
      mfile.write("de_j = perm(de_j); \n")

      j_str = "de_j = [de_j,"
      for j in range(0,nreactions):
        for i in range(0,nreactions):
          if self.reactions[reactions[j]].isDependant(self.reactions[reactions[i]]):
            i_str = i_str + str(i+1) + ','
            j_str = j_str + str(j+mspecies+1) + ','
      mfile.write(i_str[:-1] + "];\n")
      mfile.write(j_str[:-1] + "];\n")
      mfile.write("umod.G = sparse(de_i, de_j, 1);\n")

    # INITIAL VALUES
    # When dealing with initial values we will add together all values per type per compartment.
    mfile.write("\n% INITIAL VALUES\n")
    mfile.write("% Ndomains x Mspecies with initial values\n")
    mfile.write("% init_values is a #domains x #species matrix with the initial\n")
    mfile.write("% distribution of species per compartment. Default value is zero.\n")

    mfile.write("init_values = [")
    for stype_key in stypes:
      row = []
      for i in range(0,ndomains):
        row.append(0)
      stype = self.getSpeciesType(stype_key)
      for spec in stype.members:
        if not spec:
          print spec_key
        col = spec.getCompartment().getSd()-1
        row[col] = row[col] + spec.getInitialValue()
      for i in range(0,ndomains):
        mfile.write(str(row[i]))
        if i < ndomains-1:
          mfile.write(", ")
      mfile.write(";\n               ")
    mfile.write("]';\n")

    mfile.write("\n% The following snippet will generate a randomized distribution of the initial values,\n")
    mfile.write("% scattered evenly (standard uniform distribution) over the voxels of \n")
    mfile.write("% given compartment with account taken for the voxel volume.\n")
    mfile.write("%\n")
    mfile.write("% If you want some other initial distribution, please replace the code\n")
    mfile.write("% below.\n\n")

    mfile.write("u0 = zeros(Mspecies,Ncells);\n\n")

    mfile.write("for sd = 1:" + str(ndomains) + "\n")
    mfile.write("    voxels = find(umod.sd == sd);\n")
   # mfile.write("    distro = zeros(Mspecies,size(voxels,2));\n")
    mfile.write("    vol_sum = [0; cumsum(umod.vol(voxels))];\n")
    mfile.write("    for spec = 1:Mspecies\n")
    mfile.write("        rnd = vol_sum(end)*rand(1,init_values(sd,spec));\n")
    mfile.write("        cnt = histc(rnd,vol_sum);\n")
    mfile.write("        u0(spec,voxels) = cnt(1:end-1);\n")
    mfile.write("    end\n")
    mfile.write("end\n")
    mfile.write("umod.u0 = u0;\n")

    # DATA ARRAY
    mfile.write("\n% DATA ARRAY\n")
    mfile.write("% The data array in urdme is a collection of evaluated expressions used to send\n")
    mfile.write("% data into the propensity functions. When using sbml2rdme any parameter that does not\n")
    mfile.write("% have a set value is expected to be evaluated and put in the data array. The order, index\n")
    mfile.write("% of the array is very important as they are matched with the propensity functions in\n")
    mfile.write("% " + self.getName() + ".c\n\n")
    
    if len(self.data) > 0:

      mfile.write("% This default interpretation assumes that the variables can be evaluated by the corresponding\n")
      mfile.write("% model in Comsol Multiphysics by the same name (read ID) give in the sbml model. Make sure\n")
      mfile.write("% you look through these definitions and that they are correctly interpreted by Comsol.\n\n")

      mfile.write("% The top alternative uses the API used by Comsol 4.2a and the lower one uses the older\n")
      mfile.write("% Comsol 3.5a interface.\n\n")

      mfile.write("if iscmp4x(umod.comsol) %Comsol 3.5 or 4.x?\n")
      mfile.write("  xmi = mphxmeshinfo(umod.comsol);\n")
      for data_key in self.data:
        param = self.getData(data_key)
        #mfile.write("  % data[" + str(param.data) + "] = " + str(param.getId()) +"\n")
        mfile.write("  umod.data(" + str(param.data+1) + ",:) = mphinterp(umod.comsol,'" + str(param.getId()) + "','coord', xmi.dofs.coords(:,1:Mspecies:end), 'solnum', 1);\n")
      mfile.write("  umod.data = umod.data(xmi.dofs.nodes(1:Mspecies:end)+1);\n")
      mfile.write("else\n")
      mfile.write("  dofs = xmeshinfo(umod.comsol,'Out','dofs');\n")
      for data_key in self.data:
        param = self.getData(data_key)
        #mfile.write("  % data(" + str(param.data+1) + ") = " + str(param.getId()) +"\n")
        mfile.write("  umod.data = postinterp(umod.comsol,'" + str(param.getId()) + "',dofs.coords(:,1:Mspecies:end));\n")
      mfile.write("  umod.data = umod.data(dofs.nodes(1:Mspecies:end));\n")
      mfile.write("end\n")
    else:
      mfile.write("% In the given sbml model there are no such expressions defined. So we leave it as an empty vektor.\n\n")
      mfile.write("umod.data = zeros(0,Ncells);\n")

    # TIME INTERVAL
    mfile.write("\n% TIME INTERVAL\n")
    mfile.write("% The number of time steps and their sizes that urdme will use for evaluation is not a feature present in\n")
    mfile.write("% SBML lvl 2 so it is here set at a default value of 1:10 (1 to 10). Simply alter the tspan vektor for other\n")
    mfile.write("% settings.\n\n")
    
    mfile.write("umod.tspan = [1:10];\n")
    
    mfile.close()

# END MODEL CLASS
###

## CompartmentType class
# Represents the types of compartments in case there are several compartments sharing the 
# same properties (i.e same sd)
class CompartmentType():
  def __init__(self, cid, sd = None, name = None):
    self.cid = propperId(cid)
    self.sd = int(sd)
    self.name = str(name)
    self.members = [] # empty list

  def getSd(self):
    return self.sd
  
  def getId(self):
    return self.cid

  def addMember(self, comp):
    if isinstance(comp, Compartment):
      self.members.append(comp)
    else:
      model_error("Attempting to add a non compartment (" + comp.getId() + ") as member in compartment type " + str(self.cid))
      # will exit(1)


## Species Type
#  Species type 
#
class SpeciesType():
  def __init__(self, sid, member=None, auto=False,  name=None):
    self.sid = propperId(sid)
    self.members = []    # empty list
    self.dependancy = [] # empty list
    self.auto = auto
    if name and str(name).strip():
      self.name = str(name).strip()
    else:
      self.name = None
    if member:
      self.addMember(member)
      
  def getId(self):
    return self.sid

  def getName(self):
    return self.name

  def addMember(self, member):
    if not isinstance(member, Species):
      model_warning("Ignoring attempt to add a non species as member in species type " + self.sid)  
    else:
      self.members.append(member)

  def isDependant(self, reaction):
    return reaction in self.dependancy

  def addDependancy(self, reaction):
    if not isinstance(reaction, Reaction):
      model_warning("Ignoring attempt to add non reaction as dependency in " + str(self.getId()))
    elif not self.isDependant(reaction):
      self.dependancy.append(reaction)

## Compartment class
#
class Compartment():
  def __init__(self, cid, sd = None, name = None, ctype = None, vol = 0):
    self.cid = propperId(cid)

    if isinstance(ctype, CompartmentType):
      self.ctype = ctype
    else:
      model_warning("Attempting to set a non compartment type in compartment " + self.cid)
      self.ctype = None

    if not sd: 
      if self.ctype and self.ctype.getSd():
        self.sd = self.ctype.getSd()
      else:
        self.sd = None
    else:
      self.sd = int(sd)

    self.name = str(name)
    self.volume = vol
    
  def getId(self):
    return self.cid

  def getType(self):
    return self.ctype

  def getVolume(self):
    return self.volume

  def getSd(self):
    if self.sd:
      return self.sd
    elif self.ctype:
      return self.ctype.getSd()
    else:
      return None

## Species class
#  Representation of the pairing of a SpeciesType within a Compartment.
#  Not to be comfused with SpeciesType
#  @see SpeciesType
#  @see Compartment
class Species():
  def __init__(self, sid, comp = None, stype = None, init=0):
    self.sid = propperId(sid)
    if not comp:
      model_error("Attempting to set None as compartment in species " + self.sid)
      #will exit(1)
    if isinstance(comp, Compartment):
      self.compartment = comp
    else:
      model_error("Attempting to set a non compartment in species " + self.sid)
      #will exit(1)
    if isinstance(stype, SpeciesType):
      self.stype = stype
    else:
      model_warning("Attempting to set a non species type in species " + self.sid)
      self.stype = None
    self.initialValue = int(init)

  def getId(self):
    return self.sid

  def getType(self):
    return self.stype

  def getCompartment(self):
    return self.compartment

  def getInitialValue(self):
    return self.initialValue

  def setInitialValue(self, init):
    self.initialValue = self.initialValue + int(init)

class Parameter():
  keywords = ['x', 't', 'vol', 'data', 'sd']

  def __init__(self, pid, value=None, globalparam=False, data=-1):
    self.pid = propperId(pid)
    if value:
      self.value = float(value)
    else:
      self.value = None

    if globalparam == True:
      self.globalparam = True
    else:
      self.globalparam = False

    self.data = data

  def getId(self):
    return self.pid

  def isGlobal(self):
    return self.globalparam

  def isSet(self):
    return self.value != None

  def getValue(self):
    return self.value

  def isData(self):
    return int(self.data) >=0 

  def isKeyword(self):
    return self.pid in self.keywords

class Reaction():
  def __init__(self, rid):
    self.rid = propperId(rid)
    self.reactants = {}
    self.products = {}
    self.locals = {}
    self.kinetic = ''
    self.compartment = None

  def getId(self):
    return self.rid

  def addReactant(self, spec, st=1):
    if not isinstance(spec, SpeciesType):
      model_error("Attempted to add not Species type (" + spec.getId() + ") as reactant in reaction " + self.getId())
      # will exit(1)
    if math.isnan(st):
      st = 1.0
    if not float(st).is_integer():
      model_warning("Stoichiometry is not an integer and has been floored. Truncation error may occur")
    st = math.trunc(st)
    react_tuple = self.getReactant(spec.getId())
    if react_tuple:
      model_warning("Multiple addition of reactant " + spec.getId() + ". Adding stochiometry of all entries.")
      react_tuple = (spec, st + react_tuple[1])
    else:
      react_tuple = (spec, st)
    self.reactants[spec.getId()] = react_tuple

  def addProduct(self, spec, st=1):
    if not isinstance(spec, SpeciesType):
      model_error("Attempted to add not Species type (" + spec.getId() + ") as product in reaction " + self.getId())
      # will exit(1)
    if math.isnan(st):
      st = 1.0
    if not float(st).is_integer():
      model_warning("Stoichiometry is not an integer and has been floored. Truncation error may occur")
    st = math.trunc(st)
    prod_tuple = self.getProduct(spec.getId())
    if prod_tuple:
      model_warning("Multiple addition of product " + spec.getId() + ". Adding stochiometry of all entries.")
      prod_tuple = (spec, st + prod_tuple[1])
    else:
      prod_tuple = (spec, st)
    self.products[spec.getId()] = prod_tuple

  def getReactant(self, sid):
    if sid and sid in self.reactants:
      return self.reactants[sid]
    else:
      return None

  def getProduct(self, sid):
    if sid and sid in self.products:
      return self.products[sid]
    else:
      return None

  def addLocal(self, parameter):
    if not isinstance(parameter, Parameter):
      model_error("Attempted to add non parameter type (" + str(parameter.getId()) + ") as local in reaction " + str(self.getId()))
    if self.locals.get(parameter.getId()):
      model_warning("Reassigning local parameter " + str(parameter.getId()) + "in reaction " + str(self.getId()))
    self.locals[parameter.getId()] = parameter

  def getLocal(self, pid):
    if pid and pid in self.locals:
      return self.locals.get(pid)
    else:
      return None

  def setKinetic(self, math):
    self.kinetic = str(math)

  def getKinetic(self):
    return self.kinetic

  ## getStoichiometry
  # Calculate the stochiometry for a reaction given a species.
  # @param sid Species Type ID for which the stochiometry 
  # @return The stoichiometry value for the species
  def getStoichiometry(self, sid):
    stoichiometry = 0
    reactant = self.getReactant(sid)
    if reactant:
      stoichiometry = stoichiometry - reactant[1]
    product = self.getProduct(sid)
    if product:
      stoichiometry = stoichiometry + product[1]

    return stoichiometry

  def isDependant(self, reaction):
    for reactant in self.reactants:
      if self.reactants[reactant][0].isDependant(reaction):
        return True
    for product in self.products:
      if self.products[product][0].isDependant(reaction):
        return True
    return False
  
  def getCompartment(self):
    return self.compartment

  def setCompartment(self, comp):
    if not isinstance(comp, CompartmentType):
      model_error("Attempting to add non compartment type (" + str(comp.getId()) + ") in reaction " + str(self.getId()))
      #will exit(1)
    if self.compartment:
      model_warning("Reassinging compartment type in " + str(self.getId()))
    self.compartment = comp

##
#  UTILITY FUNCTIONS
##

def propperId(sid = None):
  if sid and str(sid).strip():
    return str(sid).strip()
  else:
    model_error(str(sid) + " is not a proper Id")
    #will exit(1)

def model_error(e_string):
  if e_string:
    sys.stderr.write("Error: " + str(e_string) + "\n")
    exit(1)

def model_warning(w_string):
  if w_string:
    sys.stderr.write("Warning: " + str(w_string) + "\n")

def model_openFile(filename, opts = None):
  retfile = None
  try:
    retfile = open(filename, 'r')
    retfile.close()
    model_error(str(filename) + " does already exist. Please remove previous version before creating a new one.")
  except IOError:
    if retfile:
      retfile.close()
    return open(filename, opts)
