#! /usr/bin/python

## @package sbml2rdme
#
#  SBML2RDME is a tool for creating the neccessary model files (model.m & model.c) for URDME from
#  a model described in SBML (www.sbml.org). SBML2RDME was first released with URDME 1.2 
#
#  SBML2RDME is uses the official sbml library for Python available 
# 
#  @author V.Gerdin


import sys            # System environmental access
from libsbml import * # Official SBML library from sbml.org
import urdme_model    # Model classes, found in src/urdme_model.py

model = urdme_model

## Version
# Returns the version number of SBML2RDME
# @return string with version label
def version():
  return "SBML2RDME 0.1 alpha"

## Options parser
# Uses a standard module to handle command line arguments. The module depends on the running version of Python.
# As of Python version 2.7 the module optparse was depreciated and replaced by argparse.
# @return Returns an object with the parsed arguments available by name.   
def getOptions():
  if sys.version_info >= (2,7):
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Turn SBML model into URDME 1.2 model files.')
    parser.add_argument('-f','--file', action='store',help="Input file name (SBML format)", default = "model.xml", dest="input")
    parser.add_argument('-o','--output',action='store', help="Output destination (directory)", default=".", dest="output")
    return parser.parse_args()
  else:
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f","--file", help="Input file name (SBML format)", type="string", default = "model.xml", dest="input")
    parser.add_option("-o","--output", help="Output destination (directory)", type="string", default=".", dest="output")
    (opts, args) = parser.parse_args()
    return opts

## sbml2rdme parser
# The parser function traverses the sbml model building the rdme object as it goes. 
# The parsing is done in the order the parts are expected by sbml level 2 standard but
# should not matter if out of order.
# @param sbml An sbml model object to parse.
# @return rdme An rdme model object, ready to print
def sbml2rdme(sbml):
  if not isinstance(sbml, Model):
    parser_error("Inproper SBML model used for translating")
    # will exit(1)
  
  # beginning of model definition
  # create urdme model object 
  rdme = model.model(name = sbml.getName())

  ###
  #  START LIST OF FUNCTION DEFINITIONS (not supported)
  sbml_functions = sbml.getListOfFunctionDefinitions()

  if sbml_functions.size() > 0:
    parser_warning(version() +" does not support the Function Definintion component. Ignoring it.")

  del sbml_functions
  #  END LIST OF FUNCTION DEFINITIONS
  ###

  ###
  #  START LIST OF UNIT DEFENITIONS (not supported)
  sbml_units = sbml.getListOfUnitDefinitions()

  if sbml_units.size() > 0:
    parser_warning(version() +" does not support the Unit Definintion component. Ignoring it.")

  del sbml_units
  #  END LIST OF UNIT DEFENITIONS
  ###

  ###
  # START LIST OF COMPARTMENT TYPES (optional)
  # this list contains the types of compartments, where every member compartment has the same
  # sd value, so urdme can't tell them apart.
  ctypes = sbml.getListOfCompartmentTypes()
  for i in range(0,ctypes.size()):
    ctype = ctypes.get(i)
    rdme.addCompartmentType(model.CompartmentType(ctype.getId(),sd = i+1,name = ctype.getName()))
    # sd = i+1 in order to make unit offset
  del ctypes
  # END LIST OF COMPARTMENT TYPES
  ###

  ###
  # START LIST OF SPECIES TYPES (optional)
  stypes = sbml.getListOfSpeciesTypes()
  for i in range(0,stypes.size()):
    stype = stypes.get(i)
    rdme.addSpeciesType(model.SpeciesType(stype.getId(),name=stype.getName()))
  del stypes
  # END LIST OF SPECIES TYPES
  ###
  
  ###
  # START LIST OF COMPARTMENTS (compulsory)
  # this list is the actual compartments. Though urdme can't tell them appart except by
  # compartment type (with sd) multiple compartments of the same type will be redundant and
  # added together when calculating initial concentrations.  
  comps = sbml.getListOfCompartments()

  if comps.size() == 0:
    parser_warning("The given model does not contain any list of compartments. To urdme this makes no sense.")

  for i in range(0,comps.size()):
    w_string = ""
    sbml_comp = comps.get(i)
    comp_type_attr = sbml_comp.getCompartmentType()
    comp_type = None
    if not comp_type_attr:
      #comp_type_attr not set in sbml, will use id as type
      comp_type_attr = sbml_comp.getId()
      w_string = w_string + "Compartment type attribute missing for compartment " + str(sbml_comp.getId()) + "\n\t"
      w_string = w_string + "Using " + str(sbml_comp.getId()) + "instead.\n"

    comp_type = rdme.getCompartmentType(comp_type_attr)

    if not comp_type:
      # if comp_type is None (not present or registered in the model) we create a new compartment type
      w_string = w_string + "Unknown compartment type (" + str(comp_type_attr) + ") for compartment " + str(sbml_comp.getId()) + "\n\t"
      w_string = w_string + "Creating a new compartment type with id " + str(comp_type_attr) + "\n"

      comp_type = model.CompartmentType(comp_type_attr,rdme.getNextSD(),sbml_comp.getName())
      rdme.addCompartmentType(comp_type)

    if sbml_comp.isSetVolume():
      vol = sbml_comp.getVolume()
    else:
      vol = 0.0

    rdme.addCompartment(model.Compartment(sbml_comp.getId(),name=sbml_comp.getName(),ctype=comp_type, vol=vol))
  del comps
  # END LIST OF COMPARTMENTS
  ###

  ###
  # START LIST OF SPECIES (compulsory)
  # an sbml species is a connection between a species type and a compartment
  specs = sbml.getListOfSpecies()

  if specs.size() == 0:
    parser_warning("The given model does not contain any list of species. To urdme this makes no sense.")

  for i in range(0,specs.size()):
    w_string = ""
    sbml_spec = specs.get(i)
    spec_type_attr = sbml_spec.getSpeciesType()
    rdme_stype = None
    if not spec_type_attr:
      # species type has NOT been defined for a species in the xml file. using species id as species type id
      w_string = w_string + "Species type attribute missing for " + str(sbml_spec.getId()) + "\n\t"
      w_string = w_string + "Using species type id " + sbml_spec.getId() + " instead.\n"
      spec_type_attr = sbml_spec.getId()
    
    rdme_stype = rdme.getSpeciesType(spec_type_attr)

    if not rdme_stype:
      w_string = w_string + "Unknown species type (" + str(spec_type_attr) + ") for species " + str(sbml_spec.getId()) + "\n\t"
      w_string = w_string + "Creating new species type with id " + str(spec_type_attr) + "\n"       
      rdme_stype = model.SpeciesType(spec_type_attr,name=sbml_spec.getName())
      rdme.addSpeciesType(rdme_stype)
    
    if w_string:
      parser_warning(w_string)

    rdme_comp = rdme.getCompartment(sbml_spec.getCompartment())

    initVal = 0
    if sbml_spec.isSetInitialAmount():
      initVal = sbml_spec.getInitialAmount()
    elif sbml_spec.isSetInitialConcentration():
      # not ammount but concentration is set. check if compartment volume is proper. if so, multiply, else set 0.
      if rdme_comp.getVolume() > 0:
        initVal = int(rdme_comp.getVolume() * sbml_spec.getInitialConcentration())
      else:
        initVal = 0
    rdme.addSpecies(model.Species(sbml_spec.getId(), rdme_comp, rdme_stype, initVal))
  del specs
  # END LIST OF SPECIES
  ###

  ###
  # START LIST OF PARAMETERS (optional)
  sbml_params = sbml.getListOfParameters()
  
  for i in range(0, sbml_params.size()):
    sbml_param = sbml_params.get(i)
    if sbml_param.isSetValue():
      rdme.addGlobal(model.Parameter(pid=sbml_param.getId(),value=sbml_param.getValue(), globalparam=True))
    else:
      rdme_param = model.Parameter(pid=sbml_param.getId(), globalparam=True)
      rdme.addGlobal(rdme_param)
      if not param_obj.isKeyword():
        rdme.addData(rdme_param)

  del sbml_params
  # END LIST OF PARAMETERS
  ###

  ###
  # START LIST OF INITIAL ASSIGNMENTS (not supported)  
  sbml_assignments = sbml.getListOfInitialAssignments()

  if sbml_assignments.size() > 0:
    parser_warning(version() + " does not support the Initial Assignment component. Ignoring it.")

  del sbml_assignments
  # END LIST OF INITIAL ASSIGNMENTS
  ###

  ###
  # START LIST OF RULES (not supported)  
  sbml_rules = sbml.getListOfRules()

  if sbml_rules.size() > 0:
    parser_warning(version() + " does not support the Rules component. Ignoring it.")

  del sbml_rules
  # END LIST OF RULES
  ###

  ###
  # START LIST OF CONSTRAINTS (not supported)  
  sbml_constraints = sbml.getListOfConstraints()

  if sbml_constraints.size() > 0:
    parser_warning(version() + " does not support the Constraints component. Ignoring it.")

  del sbml_constraints
  # END LIST OF CONSTRAINTS
  ###

  ###
  # START LIST OF REACTIONS (compulsory)
  reactions = sbml.getListOfReactions()
  if reactions.size() == 0:
    parser_warning("The given model does not contain any list of reactions. To urdme this makes no sense.")

  for i in range(0, reactions.size()):
    sbml_reaction = reactions.get(i)
    rdme_reaction = model.Reaction(sbml_reaction.getId())

    comp_sbml = sbml_reaction.getCompartment()
    if comp_sbml:
      comp = rdme.getCompartment(comp_sbml)
      if not comp:
        parser_warning("Reaction occurs in non defined compartment type. Ignoring it.")
      else:
        rdme_reaction.setCompartment(comp)

    reactants = sbml_reaction.getListOfReactants()
    if reactants:
      for ii in range(0,reactants.size()):
        w_string = ""
        reference = reactants.get(ii)
        spec = rdme.getSpecies(reference.getSpecies())
        if not spec:
          # Reactant id not found in species, ignore.
          w_string = w_string + "Reactant id (" + reference.getSpecies() + ") is not a registered Species.\n\t"
          w_string = w_string + "Ignoring reactant."
        else:
          # Reactant id found in species
          rdme_reaction.addReactant(spec.getType(), reference.getStoichiometry())
          if not rdme_reaction.getCompartment():
            rdme_reaction.setCompartment(spec.getCompartment().getType())
          elif not rdme_reaction.getCompartment() is spec.getCompartment().getType():
            parser_warning("Multiple compartments used in reaction " + rdme_reaction.getId() + ". Will use " + rdme_reaction.getCompartment().getId())
        if w_string:
          parser_warning(w_string)
      del reactants

    products = sbml_reaction.getListOfProducts()
    if products:
      for ii in range(0,products.size()):
        w_string = ""
        reference = products.get(ii)
        spec = rdme.getSpecies(reference.getSpecies())
        if not spec:
          # Product id not found in species, ignore
          w_string = w_string + "Product id (" + reference.getSpecies() + ") is not a registered Species.\n\t"
          w_string = w_string + "Ignoring reactant."
        else:
          # Product id found in species
          rdme_reaction.addProduct(spec.getType(), reference.getStoichiometry())
          if not rdme_reaction.getCompartment():
            rdme_reaction.setCompartment(spec.getCompartment().getType())
          elif not rdme_reaction.getCompartment() is spec.getCompartment().getType():
            parser_warning("Multiple compartments used in reaction " + rdme_reaction.getId() + ". Will use " + rdme_reaction.getCompartment().getId())
        if w_string:
          parser_warning(w_string)
      del products

    kinLaw = sbml_reaction.getKineticLaw()
    params = kinLaw.getListOfParameters()
  
    for i in range(0, params.size()):
      param = params.get(i)
      if param.isSetValue():
        rdme_reaction.addLocal(model.Parameter(pid=param.getId(),value=param.getValue(), globalparam=False))
      else:
        param_obj = model.Parameter(pid=param.getId(), globalparam=False)
        rdme_reaction.addLocal(param_obj)
        if not param_obj.isKeyword():
          rdme.addData(param_obj)

    del params
    
    # math is and ASTnode tree which SBML uses to represent the kinetic formula for a reaction
    # a simple way to traverse the tree is by either depth-first och breadth-first. Either of which
    # will be equally fast as we are traversing the entire tree.
    math = kinLaw.getMath()
    queue = [] # empty list

    while math:
      if math.isName():
        varname = str(math.getName())
        spec = rdme.getSpecies(varname)
        if spec:
          stype = spec.getType()
          if stype:
            stype_id = stype.getId()
            if not (stype_id in rdme_reaction.reactants or stype_id in rdme_reaction.products):
              parser_warning("Species refered to in reaction " + str(rdme_reaction.getId()) + " is not defined in its reactants or products.")
            # if varname is a Species: add dependancy (as the reactions kinetics are dependant on the species)
            # and: rewrite the variable to work for the model.c propensity format.
            stype.addDependancy(rdme_reaction)
            math.setName('x[' + stype.getId() +']')
        else:
          comp = rdme.getCompartment(varname)
          if comp:
            if comp.getType() is rdme_reaction.getCompartment():
              # compartment references means their volume in SBML so we use this as voxel volume.
              math.setName("vol")
            else:
              parser_warning("Compartment reference (" + str(comp.getId()) + ") in reaction " + str(rdme_reaction.getId()) + " does not match the reaction.") 
          else:
            # if varname is a parameter in (local or global) we check if it is in the data array.
            # if in the data array we change the entry to data[#]  
            param = rdme_reaction.getLocal(varname)
            if not param:
              param = rdme.getGlobal(varname)
            if not param:
              parser_warning("Found non-identified parameter (" + str(varname) + ") in reaction " + str(rdme_reaction.getId()))
            elif param.isData():
              math.setName('data[' + str(param.data) + ']')
            
                
      # below is neccessary for the BFS to continue
      if math.getNumChildren() > 0:
        for i in range(0,math.getNumChildren()):
          queue.append(math.getChild(i)) # for DFS, change to queue.insert(0,math.getChild(i))
      if len(queue) > 0:
        math = queue.pop()
      else:
        math = None
    #end of 'while math:'
      
    rdme_reaction.setKinetic(kinLaw.getFormula())

    rdme.addReaction(rdme_reaction)

  del reactions
  # END LIST OF REACTIONS
  ###
  
  return rdme

## Main function
# The main workflow: read options, read sbml-file, create sbml_model object from said file, 
# parse sbml_model and create rdme_model, rdme_model creates model.c and model.m files.
def main():

  # Read input arguments
  opts = getOptions()

  if opts:
    infile = opts.input
    outdir = opts.output
  else:
    infile = 'model.xml'
    outdir = './'

  # Create SBMLdocument object with SBML func 
  reader = SBMLReader()
  document = reader.readSBML(infile)

  # PROCESS ANY SBML ERRORS
  num_errors = document.getNumErrors()

  if (num_errors):
    w_string = str(document.getNumErrors()) + " error(s) were found when reading " + infile+ ".\n\t"
    w_string = w_string + "Read below for details.\n"
    for i in range(0,document.getNumErrors()):
      w_string = w_string + "\t" + str(document.getError(i).getMessage()) + "\n"
   
  # GET MODEL FROM DOCUMENT
  sbml_model = document.getModel()

  # REWORK MODEL
  rdme_model = sbml2rdme(sbml_model)

  # CREATE model.c
  rdme_model.c(outdir)

  # CREATE model.m
  rdme_model.m(outdir)

## Error handler
# parser_error() will write an error string to stderr and raise a system exit exception.
# @param e_string A string that describes the error
def parser_error(e_string):
  if e_string:
    sys.stderr.write("Error: " + str(e_string) + "\n")
    exit(1)

## Warning handler
# parser_warning() will write a warning string to stderr. Program will continue afterwards.
# @param w_string A string that describes the warning. 
def parser_warning(w_string):
  if w_string:
    sys.stderr.write("Warning: " + str(w_string) + "\n")


if __name__ == '__main__':
  main()
